import logging
from collections.abc import Callable

import defopt
import igvf_utils

from igvf_portal import utils
from igvf_portal.constants import VERSION as VERSION
from igvf_portal.tools.download_file import download_file
from igvf_portal.tools.gen_upload_script import (
    gen_upload_script,
)
from igvf_portal.tools.get_url import get_url
from igvf_portal.tools.infer_principal_analysis import infer_principal_analysis
from igvf_portal.tools.make_pseudobulk_tracker import make_pseudobulk_tracker

tools: list[Callable] = [
    gen_upload_script,
    get_url,
    infer_principal_analysis,
    download_file,
    make_pseudobulk_tracker,
]


def fix_igvf_logging(level: int = logging.WARNING):
    utils.setup_logger(igvf_utils.debug_logger, level)
    utils.setup_logger(logging.getLogger(name="root"), level)


def main() -> None:
    fix_igvf_logging(level=logging.WARNING)
    defopt.run(tools)
