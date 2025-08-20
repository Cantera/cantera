# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from pathlib import Path as _Path
from ._utils import get_data_directories, add_directory


def list_data_files(ext=".yaml"):
    """
    Lists input data files. Includes files in subdirectories, except for subdirectories
    of the current working directory.

    :param ext:
        Extension of files to be displayed.
    :return:
        List of input data files.
    """
    data_files = set()
    for folder in get_data_directories():
        here = _Path(folder)
        if folder == ".":
            data_files.update(f.name for f in here.glob(f"*{ext}"))
        elif here.is_dir():
            data_files.update(str(f.relative_to(here)) for f in here.glob(f"**/*{ext}"))
    return sorted(data_files)
