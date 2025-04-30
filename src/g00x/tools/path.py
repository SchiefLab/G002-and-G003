import os
from pathlib import Path
from typing import Any

from pandera.typing import DataFrame, Index, Series


class cd:
    """Context manager for changing the current working directory"""

    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def pathing(path: Path | str, new: bool = False, overwrite: bool = True) -> Path:
    """Guarantees correct expansion rules for pathing.
    :param Union[str, Path] path: path of folder or file you wish to expand.
    :param bool new: will check if distination exists if new  (will check parent path regardless).
    :return: A pathlib.Path object.
    >>> pathing('~/Desktop/folderofgoodstuffs/')
    /home/user/Desktop/folderofgoodstuffs
    """
    path = Path(path)
    # Expand shortened path
    if str(path)[0] == "~":
        path = path.expanduser()
    # Exand local path
    if str(path)[0] == ".":
        path = path.resolve()
    else:
        path = path.absolute()
    # Making sure new paths don't exist while also making sure existing paths actually exist.
    if new:
        if not path.parent.exists():
            raise ValueError(f"ERROR ::: Parent directory of {path} does not exist.")
        if path.exists() and not overwrite:
            raise ValueError(f"ERROR ::: {path} already exists!")
    else:
        if not path.exists():
            raise ValueError(f"ERROR ::: Path {path} does not exist.")
    return path


def pd_replace_home_with_tilde(x: Series[Any]) -> Series[Any]:
    """Replaces home path with tilde for dynamic pathing per user for every cell.

    Parameters
    ----------
    x : Series[Any]
        A pandas series of Any.
    Returns
    -------
    A pandas series of strings containing the tilde if the path is in the home directory.

    >>> df.applymap(pd_replace_home_with_tilde)
    """

    if not isinstance(x, (Path, str)):
        return x
    if len(x) > 256:
        return x
    try:
        if pathing(x, new=True).exists():
            return str("~" / pathing(x).relative_to(Path.home()))
    except ValueError:
        return x
    return x


def pd_expand_path(x: Series[Any]) -> Series[Any]:
    """Expands home path with tilde for dynamic pathing per user for every cell.

    Parameters
    ----------
    x : Series[Any]
        A pandas series of Any.
    Returns
    -------
    A pandas series of expanded paths if path contains '.' or '~'

    >>> df.applymap(pd_expand_path)
    """

    if not isinstance(x, (Path, str)):
        return x
    if len(x) > 256:
        return x
    try:
        if pathing(x, new=True).exists():
            return str(pathing(x))
    except ValueError:
        return x
    return x
