from Cantera import CanteraError
import _cantera
from Numeric import *

class FlowBoundary:
    def __init__(self, type, phase, kin=None):
        self.phase = phase
        self.kin = kin
        if kin:
            self.__bdry_id = _cantera.bdry_new(type, phase.phase_id(), kin.kin_id())
        else:
            self.__bdry_id = _cantera.bdry_new(type, phase.phase_id(), 0)
            
        FlowBoundary.set(self, mdot = 0.0, V = 0.0, T = phase.temperature(),
                 X = phase.moleFractions())

    def __del__(self):
        _cantera.bdry_del(self.__bdry_id)
    def bdry_id(self):
        return self.__bdry_id
    def set(self, mdot = -999.0, V = -999.0, T = -999.0, X = None, Y = None):
        if mdot > 0.0:
            self.mdot = mdot
            _cantera.bdry_set(self.__bdry_id, 1, mdot, None)
        if V > 0.0:
            self.V = V
            _cantera.bdry_set(self.__bdry_id, 2, V, None)
        if T > 0.0:
            self.T = T
            _cantera.bdry_set(self.__bdry_id, 3, T, None)
        if X:
            self.phase.setMoleFractions(X)
            yy = self.phase.massFractions()
            self.X = self.phase.moleFractions()            
            _cantera.bdry_set(self.__bdry_id, 4, 0.0, yy)            
        if Y:
            self.phase.setMassFractions(Y)
            yy = self.phase.massFractions()
            self.X = self.phase.moleFractions()                        
            _cantera.bdry_set(self.__bdry_id, 4, 0.0, yy)


class Inlet(FlowBoundary):
    def __init__(self, phase):
        FlowBoundary.__init__(self, 0, phase)


class Outlet(FlowBoundary):
    def __init__(self, phase):
        FlowBoundary.__init__(self, 1, phase)
    def set(self, mdot = -999.0, V = -999.0, T = -999.0, X = None, Y = None):
        raise CanteraError("outlet properties cannot be set.")

    
class Surface(FlowBoundary):
    def __init__(self, phase, surf):
        self._surf = surf
        self.domainType = 1
        FlowBoundary.__init__(self, 2, phase, surf)
        self.__surf1d_id = _cantera.surf1d_new(self.kin.kin_id())
        self.x = zeros((1,self.kin.nSpecies()),'d')
        self.x[0,:] = self.kin.coverages()
        self.species_on = []
        self.species_off = []        
    def surf_id(self):
        return self.__surf1d_id
    def setCoverages(self, theta):
        self.kin.setCoverages(theta)
        self.x[0,:] = self.kin.coverages()
    def coverages(self):
        return self.x[0] # kin.coverages()
    def show(self):
        self.kin.setCoverages(self.x)
        self.kin.show()
    def shape(self):
        return (1, self.kin.nSpecies())
    def integrate(self, dt):
        self.kin.setCoverages(self.x)
        self.kin.integrate(dt)
        c = self.kin.coverages()
        self.x[0,:] = c
        for k in range(len(c)):
            #print k, c[k]
            self.fixSpecies(k,c[k])
        self.setSpeciesEqn(on = self.species_on, off = self.species_off)
    def fixSpecies(self, k, c):
        _cantera.surf1d_fixspecies(self.__surf1d_id, k, c)
    def setTemperature(self, t):
        _cantera.surf1d_settemperature(self.__surf1d_id, t)
    def temperature(self):
        return _cantera.surf1d_temperature(self.__surf1d_id)
    def setMultiplier(self,f):
        for k in range(self.kin.nSpecies()):
            _cantera.surf1d_setmultiplier(self.__surf1d_id,k,f)
        
    def setSpeciesEqn(self, on=None, off=None, loglevel=0):
        """Enable or disable surface species equations."""
        
        self.species_on = on
        self.species_off = off

        on_msg = `on`
        off_msg = `off`
        if on == 'all':
            on = self.kin.speciesNames()
            on_msg = 'all species'
        if off == 'all':
            off = self.kin.speciesNames()
            off_msg = 'all species'
            
        a = zeros(self.kin.nSpecies(),'d')
        
        if loglevel > 0:
            if on:
                print '\n\n--------------------------------------------------\n'
                print 'Enabling species equation for',on_msg
                print '\n--------------------------------------------------\n\n'                
            if off:
                print '\n\n--------------------------------------------------\n'
                print 'Disabling species equation for',off_msg
                print '\n--------------------------------------------------\n\n'                                                
        if on:
            for sp in on:
                k = self.kin.speciesIndex(sp)
                a[k] = 1.0
        if off:
            for sp in off:
                k = self.kin.speciesIndex(sp)
                a[k] = 0.0
                
        _cantera.surf1d_solvespecies(self.__surf1d_id, len(a), a)        



class SymmPlane(FlowBoundary):
    def __init__(self, phase):
        FlowBoundary.__init__(self, 3, phase)
    def set(self, mdot = -999.0, V = -999.0, T = -999.0, X = None, Y = None):
        raise CanteraError("outlet properties cannot be set.")

