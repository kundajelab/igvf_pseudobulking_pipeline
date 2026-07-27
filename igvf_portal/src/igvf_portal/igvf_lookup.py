import dataclasses
from collections.abc import Iterable
from typing import cast

from igvf_portal.connection import PConnection
from igvf_portal.constants import VERSION
from igvf_portal.enums import (
    AnalysisStep,
    IgvfMode,
)
from igvf_portal.types import (
    AccessionId,
    Alias,
    IgvfRecord,
)


@dataclasses.dataclass(kw_only=True, slots=True)
class IgvfLookup:
    connection: PConnection
    igvf_mode: IgvfMode
    record_lookups: dict[AccessionId | Alias, IgvfRecord] = dataclasses.field(
        default_factory=dict
    )

    @classmethod
    def new(cls, igvf_mode: IgvfMode | str) -> IgvfLookup:
        _igvf_mode = (
            igvf_mode if isinstance(igvf_mode, IgvfMode) else IgvfMode[igvf_mode]
        )
        return cls(connection=PConnection(igvf_mode=_igvf_mode), igvf_mode=_igvf_mode)

    def lookup_record(self, key: AccessionId | Alias) -> IgvfRecord:
        """Lookup record, using table of previous lookups if it's present, otherwise adding to table."""
        record = self.record_lookups.get(key, None)
        if record is None:
            record = self.connection.get(rec_ids=key, ignore404=False)
            if record is None:
                raise ValueError(f"Could not find record for '{key}'")
            self.record_lookups[key] = record
        return record

    def infer_principal_accessions(
        self, intermediate_accessions: Iterable[AccessionId]
    ) -> set[AccessionId]:
        """Check intermediate accessions to find the principal accessions they derive from."""
        # check all the supplied intermediate accessions
        to_check = set(intermediate_accessions)
        principal_accessions: set[AccessionId] = set()
        while len(to_check) > 0:
            # pop off one of the intermediate accessions and get its record
            intermediate_accession = to_check.pop()
            intermediate_record = self.lookup_record(intermediate_accession)
            if intermediate_record["status"] == "deleted":
                continue
            principal_ids = intermediate_record.get("input_for", None)
            if principal_ids is None:
                # it's not input for anything, so it must be a principal accession
                principal_accessions.add(AccessionId(intermediate_accession))
            else:
                # it's input for these principal ids.
                for principal_id in principal_ids:
                    # get the record for this principal ID
                    principal_record = self.lookup_record(principal_id)
                    if principal_record["status"] == "deleted":
                        continue
                    # get its accession and add it to the output set
                    principal_accessions.add(principal_record["accession"])
                    # to decrease lookups, remove everything that was input to it from the IDs to check
                    # (this is VERY effective for data sets with many intermediate accessions)
                    intermediate_inputs = (
                        rec["accession"]
                        for rec in principal_record.get("input_file_sets", [])
                    )
                    to_check.difference_update(intermediate_inputs)

        return principal_accessions

    def lookup_analysis_step_version(self, analysis_step: AnalysisStep) -> list[Alias]:
        step_record = cast(dict[str, object], self.lookup_record(analysis_step.value))
        for version_dict in cast(
            list[dict[str, object]], step_record["analysis_step_versions"]
        ):
            for software_versions in cast(
                list[dict[str, str]], version_dict["software_versions"]
            ):
                if (
                    software_versions["name"]
                    == f"igvf_pseudobulking_pipeline-v{VERSION}"
                ):
                    version_id = cast(Alias, version_dict["@id"])
                    return self.lookup_record(version_id)["aliases"]
        raise ValueError(
            f"Unable to find version of analysis step {analysis_step} for 'igvf_pseudobulking_pipeline-v{VERSION}'"
        )
