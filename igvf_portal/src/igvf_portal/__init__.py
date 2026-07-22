import logging
import sys
from typing import Callable

import defopt

import igvf_utils

from igvf_portal.constants import VERSION as VERSION
from igvf_portal.tools.gen_upload_script import (
    gen_upload_script,
)
from igvf_portal.tools.infer_principal_analysis import infer_principal_analysis
from igvf_portal.tools.download_file import download_file

tools: list[Callable] = [gen_upload_script, infer_principal_analysis, download_file]


def fix_igvf_logging():
    debug_logger = igvf_utils.debug_logger
    debug_logger.handlers.clear()
    debug_logger.setLevel(logging.INFO)
    ch = logging.StreamHandler(stream=sys.stdout)
    f_formatter = logging.Formatter("%(asctime)s:%(name)s:\t%(message)s")
    ch.setLevel(logging.INFO)
    ch.setFormatter(f_formatter)
    debug_logger.addHandler(ch)


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    fix_igvf_logging()
    defopt.run(tools)
