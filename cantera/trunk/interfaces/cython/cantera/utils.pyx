cdef string stringify(x):
    """ Converts Python strings to std::string. """
    # This method works with both Python 2.x and 3.x.
    tmp = bytes(x.encode())
    return string(tmp)
