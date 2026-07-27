import dataclasses
import json
import logging
from collections.abc import Collection
from typing import Callable

from igvf_portal import utils
from igvf_portal.connection import PConnection
from igvf_portal.enums import IgvfMode


@dataclasses.dataclass(slots=True, kw_only=True)
class RegisterConfig:
    igvf_mode: IgvfMode
    dry_run: bool
    profile_id: str
    num_tries: int
    delay: float
    backoff: float
    overwrite_array_values: bool
    remove_properties: list[str]
    upload_file: bool
    upload_duplicate: bool

    @property
    def rm_patch(self) -> bool:
        return len(self.remove_properties) > 0

    @property
    def cleaned_profile_id(self) -> str:
        profile_id = self.profile_id.strip("/").split("/", 1)[0].lower()
        # Multi-word profile names are hypen-separated, i.e. genetic-modifications.
        profile_id = profile_id.replace("-", "_")
        return profile_id

    @property
    def connection(self) -> PConnection:
        return PConnection.new(
            igvf_mode=self.igvf_mode, submission=True, dry_run=self.dry_run
        )

    def retry[**P, R](
        self,
        no_retry_exceptions: Collection[type] = (json.decoder.JSONDecodeError,),
        logger: logging.Logger | None = None,
    ) -> Callable[[Callable[P, R]], Callable[P, R]]:
        return utils.retry(
            num_tries=self.num_tries,
            delay=self.delay,
            backoff=self.backoff,
            no_retry_exceptions=no_retry_exceptions,
            logger=logger,
        )
