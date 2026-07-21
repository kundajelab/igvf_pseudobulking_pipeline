from importlib.metadata import version
from typing import Final

VERSION: Final[str] = version(f"{__package__}".split(".", 1)[0])
