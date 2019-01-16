from pathlib import Path


def list_data_files(ext=''):
    """
    Lists input data files.

    :param ext:
        Extension of files to be displayed.
    :return:
        List of input data files.
    """
    here = Path(__file__).parent
    data_files = [x.name for x in here.iterdir() if not x.name.startswith('__') and x.name.endswith(ext)]
    return data_files
