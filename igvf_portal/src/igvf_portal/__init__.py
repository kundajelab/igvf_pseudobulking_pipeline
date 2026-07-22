import logging
from collections.abc import Callable

import defopt

from igvf_portal import utils
from igvf_portal.constants import VERSION as VERSION
from igvf_portal.tools.download_file import download_file
from igvf_portal.tools.gen_upload_script import (
    gen_upload_script,
)
from igvf_portal.tools.get_references import get_references
from igvf_portal.tools.get_species import get_species
from igvf_portal.tools.get_url import get_url
from igvf_portal.tools.lookup_record import lookup_record
from igvf_portal.tools.make_pseudobulk_tracker import make_pseudobulk_tracker
from igvf_portal.tools.register import register

tools: list[Callable] = [
    download_file,
    gen_upload_script,
    get_references,
    get_species,
    get_url,
    lookup_record,
    make_pseudobulk_tracker,
    register,
]


def main() -> None:
    utils.fix_igvf_logging(level=logging.INFO)
    defopt.run(tools)
