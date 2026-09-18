"""Shared library setup and checks for the paper's computation scripts."""
import os
from pathlib import Path
import sys

if sys.version_info < (3, 10):
    raise RuntimeError("Use Python 3.10 or newer (see README.md).")

LIBRARY = Path(os.environ.get(
    "HOMOGENEOUS_VARIETIES_PATH",
    Path(__file__).resolve().parent.parent / "homogeneous-varieties",
)).expanduser().resolve()
if not (LIBRARY / "grassmannians.py").is_file():
    raise RuntimeError(
        f"Cannot find homogeneous-varieties at:\n  {LIBRARY}\n"
        "Download https://github.com/marcorampazzo/homogeneous-varieties "
        "and place its extracted folder, named 'homogeneous-varieties', "
        "beside 'the-last-CY-pair'.\n"
        f"The file {LIBRARY / 'grassmannians.py'} should exist.\n"
        "See README.md for the folder layout."
    )
sys.path.insert(0, str(LIBRARY))

from grassmannians import Grassmannian, RHom
from hodge import dim_dict


def nonzero(cohomology):
    """Return the total dimension in each degree with nonzero cohomology."""
    return {q: n for q, n in dim_dict(cohomology).items() if n}


def require(condition, message):
    """Report a failed mathematical check, even when Python uses optimization."""
    if not condition:
        raise AssertionError(message)
