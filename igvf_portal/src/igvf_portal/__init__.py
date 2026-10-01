import logging
from collections.abc import Callable, Iterator
from pathlib import Path

import defopt

from igvf_portal import utils
from igvf_portal.constants import VERSION as VERSION


def _iter_tools() -> Iterator[Callable]:
    """Find all the top-level commands in tools"""
    tools = __import__(f"{__package__}.tools", globals(), locals(), "tools")
    for _path in sorted(Path(f"{tools.__file__}").parent.glob("*.py")):
        _name = _path.stem
        if _name.startswith("_"):
            continue
        _mod = __import__(
            f"{__package__}.tools.{_name}",
            globals(),
            locals(),
            _name,
            0,
        )
        _tool = getattr(_mod, _name, None)
        if isinstance(_tool, Callable):
            yield _tool


commands: list[Callable] = list(_iter_tools())
"""Make top-level list of commands to execute."""


def main() -> None:
    """Use defopt to display help or launch the requested tool."""
    utils.fix_igvf_logging(level=logging.INFO)
    defopt.run(commands, version=VERSION)
