from Cantera import CanteraError
import _cantera
from Numeric import *
import types

class Boundary1D:
    def __init__(self, type):
        self.domainType = 2
        self.__bndry_id = _cantera.bndry_new(type)
        self.x = zeros((1,2),'d')
        self.x[0,1] = 300.0
        
    def __del__(self):
        _cantera.bndry_del(self.__bndry_id)

    def shape(self):
        return (1,2)

    def show(self):
        pass

    def restore(self, file='', solution=''):
        pass
    
    def bndry_id(self):
        return self.__bndry_id

    def set(self, mdot = -999.0, V = -999.0, T = -999.0, X = None, Y = None):
        if mdot > 0.0:
            self.mdot = mdot
            self.x[0,0] = mdot
            _cantera.bndry_setmdot(self.__bndry_id, mdot)
        if T > 0.0:
            self.T = T
            self.x[0,1] = T
            _cantera.bndry_settemperature(self.__bndry_id, T)
        if X:
            self.X = X
            if type(X) == types.StringType:
                _cantera.bndry_setxinbyname(self.__bndry_id, X)
            else:
                _cantera.bndry_setxin(self.__bndry_id, X)

def Inlet1D():
    return Boundary1D(1)

def Symm1D():
    return Boundary1D(2)

def Surf1D():
    return Boundary1D(3)


