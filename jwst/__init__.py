"""Process JWST data with Python."""

import re
from importlib.metadata import version

__all__ = ["__version__", "__version_commit__"]

__version__ = version(__name__)

_regex_git_hash = re.compile(r".*\+g(\w+)")
__version_commit__ = ""
if "+" in __version__:
    commit = _regex_git_hash.match(__version__).groups()  # type: ignore[union-attr]
    if commit:
        __version_commit__ = commit[0]
