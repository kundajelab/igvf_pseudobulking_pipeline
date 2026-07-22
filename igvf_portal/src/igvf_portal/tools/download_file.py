from pathlib import Path

import boto3
from botocore import UNSIGNED
from botocore.config import Config

from igvf_portal.utils import check_keys
from igvf_portal.enums import IgvfMode
from igvf_portal.igvf_lookup import IgvfLookup
from igvf_portal.types import Alias


def download_file(
    key: str, igvf_mode: IgvfMode = IgvfMode.prod, output: Path | None = None
) -> None:
    """Infer the principal analysis accession IDs from the input metadata file. Display to stdout

    Args:
        metadata_file: Path to annotations file.
        igvf_mode: Mode for accessing the IGVF Portal.
    """
    check_keys()
    igvf_lookup = IgvfLookup.new(igvf_mode=igvf_mode)
    record = igvf_lookup.lookup_record(Alias(key))
    s3_uri = record["s3_uri"]
    _, _, bucket_name, object_key = s3_uri.split("/", 3)
    # note for later:
    # public_url = f"https://{bucket_name}.s3.{region}.amazonaws.com/{object_key}"
    if output is None:
        output = Path(Path(record["submitted_file_name"]).name)
    s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))
    with output.open("wb") as f:
        s3.download_fileobj(bucket_name, object_key, f)
