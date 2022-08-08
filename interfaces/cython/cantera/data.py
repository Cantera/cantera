# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

from pathlib import Path
from ._utils import get_data_directories, add_directory


def list_data_files(ext=".yaml"):
    """
    Lists input data files.

    :param ext:
        Extension of files to be displayed.
    :return:
        List of input data files.
    """
    data_files = []
    for folder in get_data_directories():
        here = Path(folder)
        if here.is_dir():
            files = [x.name for x in here.iterdir() if x.name.endswith(ext)]
        data_files.extend(files)
    data_files = list(set(data_files))
    data_files.sort()
    return data_files
