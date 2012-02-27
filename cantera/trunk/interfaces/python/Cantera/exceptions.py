"""
Cantera exceptions
"""

import _cantera

def getCanteraError():
    """
    Get an error message generated when Cantera throws an exception.
    """
    return _cantera.get_Cantera_Error()

class CanteraError(Exception):
    def __init__(self, msg = ""):
        if msg == "":
            msg = _cantera.get_Cantera_Error()
        self.msg = msg
    def __str__(self):
        print '\n\n\n#######################  CANTERA ERROR  ######################\n'
        print '    ',self.msg
        print '\n##############################################################\n'


class OptionError(CanteraError):
    def __init__(self, msg):
        self.msg = 'Unknown option: '+msg
