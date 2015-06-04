import sys
cdef int _pythonMajorVersion = sys.version_info[0]

cdef CxxPythonLogger* _logger = new CxxPythonLogger()
CxxSetLogger(_logger)

cdef string stringify(x):
    """ Converts Python strings to std::string. """
    # This method works with both Python 2.x and 3.x.
    if isinstance(x, bytes):
        return string(<bytes>x)
    else:
        tmp = bytes(x.encode())
        return string(tmp)

cdef pystr(string x):
    cdef bytes s = x.c_str()
    if _pythonMajorVersion == 2:
        # Python 2.x
        return s
    else:
        # Python 3.x
        return s.decode()

def add_directory(directory):
    """ Add a directory to search for Cantera data files. """
    CxxAddDirectory(stringify(directory))

if get_sundials_version():
    __sundials_version__ = '.'.join(str(get_sundials_version()))
else:
    __sundials_version__ = None

__version__ = pystr(get_cantera_version())

def appdelete():
    """ Delete all global Cantera C++ objects """
    CxxAppdelete()

cdef Composition comp_map(X) except *:
    if isinstance(X, (str, unicode, bytes)):
        return parseCompString(stringify(X))

    # assume X is dict-like
    cdef Composition m
    for species,value in X.items():
        m[stringify(species)] = value
    return m

cdef comp_map_to_dict(Composition m):
    return {pystr(species):value for species,value in m.items()}
