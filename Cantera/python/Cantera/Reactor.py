"""
Zero-dimensional reactors. More text.

"""

import _cantera
from Numeric import array, zeros
import types

class ReactorBase:
    """Base class for reactors.""" 

    def __init__(self, name = '', contents = None,
                 volume = 1.0, energy = 'on',
                 type = -1, verbose = 0):
        """
        Create a new ReactorBase instance. If 'contents' is specified,
        method 'insert' is invoked. The 'type' parameter determines
        the type of C++ Reactor object that is instantiated (1 = Reactor,
        2 = Reservoir).
        """
        self.__reactor_id = _cantera.reactor_new(type)
        self._inlets = []
        self._outlets = []
        self._walls = []
        self._name = name
        self._verbose = verbose
        if contents:
            self.insert(contents)
        self.setInitialVolume(volume)
        self.setEnergy(energy)

    def __del__(self):
        """Delete the reactor instance."""
        if self._verbose:
            print 'Deleting '+self._name
        _cantera.reactor_del(self.__reactor_id)

    def __str__(self):
        return self._name

    def name(self):
        return self._name
    
    def reactor_id(self):
        """The integer index used to access the kernel reactor
        object. For internal use.  """
        return self.__reactor_id
    
    def insert(self, contents):
        """
        Insert 'contents' into the reactor. Sets the objects used to compute
        thermodynamic properties and kinetic rates.
        """
        self.contents = contents
        _cantera.reactor_setThermoMgr(self.__reactor_id, contents._phase_id)
        _cantera.reactor_setKineticsMgr(self.__reactor_id, contents.ckin)

        
    def setInitialTime(self, t0):
        """Set the initial time. Restarts integration from this time
        using the current state as the initial condition. weDefault: 0.0 s"""
        _cantera.reactor_setInitialTime(self.__reactor_id, t0)

    def setInitialVolume(self, t0):
        """Set the initial reactor volume. Default: 1.0 m^3."""
        _cantera.reactor_setInitialVolume(self.__reactor_id, t0)

    def setEnergy(self, e):
        """Turn the energy equation on or off. If off, the reactor
        temperature is held constant."""
        ie = 1
        if e == 'off' or e == 0:
            ie = 0
        if self._verbose:
            if ie:
                print 'enabling energy equation for reactor',self._name
            else:
                print 'disabling energy equation for reactor',self._name                
        _cantera.reactor_setEnergy(self.__reactor_id, ie)

    def temperature(self):
        """Temperature [K]."""
        return _cantera.reactor_temperature(self.__reactor_id)

    def density(self):
        """Density [kg/m^3]."""
        return _cantera.reactor_density(self.__reactor_id)

    def volume(self):
        """Volume [m^3]."""
        return _cantera.reactor_volume(self.__reactor_id)            

    def time(self):
        """Time [s]. The reactor time is set by method advance."""
        return _cantera.reactor_time(self.__reactor_id)

    def mass(self):
        """The total mass of the reactor contents [kg]."""
        return _cantera.reactor_mass(self.__reactor_id)

    def enthalpy_mass(self):
        """The specific enthalpy [J/kg]."""
        return _cantera.reactor_enthalpy_mass(self.__reactor_id)

    def intEnergy_mass(self):
        """The specific interhal energy [J/kg]."""
        return _cantera.reactor_intEnergy_mass(self.__reactor_id)

    def pressure(self):
        """The pressure [Pa]."""
        return _cantera.reactor_pressure(self.__reactor_id)            

    def advance(self, time):
        """Advance the state of the reactor in time from the current
        time to time 'time'. Note: this method is deprecated. See class ReactorNet."""
        return _cantera.reactor_advance(self.__reactor_id, time)

    def step(self, time):
        """Advance the state of the reactor in time from the current
        time to time 'time'. Note: this method is deprecated. See class ReactorNet."""
        return _cantera.reactor_step(self.__reactor_id, time)    
    
    def massFraction(self, k):
        """Mass fraction of species k."""
        if type(k) == types.StringType:
            kk = self.contents.speciesIndex(k)
        else:
            kk = k
        return _cantera.reactor_massFraction(self.__reactor_id, kk)

    def massFractions(self):
        nsp = self.contents.nSpecies()
        y = zeros(nsp,'d')
        for k in range(nsp):
            y[k] = self.massFraction(k)
        return y

    def moleFractions(self):
        y = self.massFractions()
        self.contents.setMassFractions(y)
        return self.contents.moleFractions()

    def inlets(self):
        return self._inlets

    def outlets(self):
        return self._outlets

    def walls(self):
        return self._walls
    
    def _addInlet(self, inlet):
        """For internal use"""
        self._inlets.append(inlet)

    def _addOutlet(self, outlet):
        self._outlets.append(outlet)

    def _addWall(self, wall):
        self._walls.append(wall)

    
class Reactor(ReactorBase):
    """
    A reactor.
    """
    def __init__(self, contents = None, name = '<reactor>',
                 volume = 1.0, energy = 'on',
                 verbose = 0):
        """
        Create a Reactor instance, and if 'contents' is specified,
        insert it.
        """
        ReactorBase.__init__(self, contents = contents, name = name,
                             volume = volume, energy = energy,
                             verbose = verbose, type = 1)
            

class Reservoir(ReactorBase):
    """
    A reservoir is a reactor with a constant state. Class Reservoir
    derives from class ReactorBase, and overloads method advance to do
    nothing.
    """
    def __init__(self, contents = None, name = '<reservoir>', verbose = 0):
        ReactorBase.__init__(self, contents = contents,
                             name = name, verbose = verbose, type = 2)
            
    def advance(self, time):
        """Do nothing."""
        pass



#------------------ FlowDevice ---------------------------------

class FlowDevice:
    """
    Base class for devices that regulate the flow rate in a fluid line.
    """
    def __init__(self, type, name, verbose):
        """
        Create a new instance of type 'type'
        """
        self._name = name
        self._verbose = verbose
        self.__fdev_id = _cantera.flowdev_new(type)

    def __del__(self):
        """
        Delete the instance.
        """
        if self._verbose:
            print 'deleting '+self._name
        _cantera.flowdev_del(self.__fdev_id)
        
    def name(self):
        return self._name
    
    def ready(self):
        """
        Returns true if the device is ready to use.
        """
        return _cantera.flowdev_ready(self.__fdev_id)

    def massFlowRate(self):
        """
        Mass flow rate (kg/s).
        """
        return _cantera.flowdev_massFlowRate(self.__fdev_id)

    def setSetpoint(self, v):
        """
        Set the set point.
        """
        _cantera.flowdev_setSetpoint(self.__fdev_id, v)

    def setpoint(self):
        """
        The setpoint value.
        """
        return _cantera.flowdev_setpoint(self.__fdev_id)

    def install(self, upstream, downstream):
        """
        Install the device between the upstream and downstream
        reactors.
        """
        if self._verbose:
            print
            print self._name+': installing between '+upstream.name()+' and '+downstream.name()
        upstream._addOutlet(self)
        downstream._addInlet(self)
        _cantera.flowdev_install(self.__fdev_id, upstream.reactor_id(),
                                  downstream.reactor_id())
    def setParameters(self, c):
        params = array(c,'d')
        n = len(params)
        return _cantera.flowdev_setParameters(self.__fdev_id, n, params)    

_mfccount = 0

class MassFlowController(FlowDevice):
    def __init__(self, upstream=None, downstream=None, name='', verbose=0):
        global _mfccount
        if name == '':
            name = 'MFC_'+`_mfccount`
        _mfccount += 1
        FlowDevice.__init__(self,1,name,verbose)
        if upstream and downstream:
            self.install(upstream, downstream)

    def setMassFlowRate(self, mdot):
        if self._verbose:
            print self._name+': setting mdot to '+`mdot`+' kg/s'
        self.setSetpoint(mdot)


_valvecount = 0

class Valve(FlowDevice):
    def __init__(self, upstream=None, downstream=None, name='', verbose=0):
        global _valvecount
        if name == '':
            name = 'Valve_'+`_valvecount`
        _valvecount += 1
        FlowDevice.__init__(self,3,name,verbose)
        if upstream and downstream:
            self.install(upstream, downstream)                

    def setValveCoeff(self, v):
        vv = zeros(1,'d')
        vv[0] = v
        if self._verbose:
            print
            print self._name+': setting valve coefficient to '+`v`+' kg/Pa-s'
        self.setParameters(vv)


#------------- Wall ---------------------------

class Wall:
    """
    A Wall separates two reactors. Any number of walls may be created
    between any pair of reactors.
    """
    def __init__(self, left=None, right=None, area=1.0):
        typ = 0
        self.__wall_id = _cantera.wall_new(typ)
        if left and right:
            self.install(left, right)
        self.setArea(area)
        self.setExpansionRateCoeff(0.0)
        self.setExpansionRate()        
        self.setHeatFlux()

    def __del__(self):
        """
        Delete the Wall instance.
        """
        _cantera.wall_del(self.__wall_id)
        
    def ready(self):
        """
        Return 1 if the wall instance is ready for use, 0 otherwise.
        """
        return _cantera.wall_ready(self.__wall_id)

    def area(self):
        """
        The wall area (m^2).
        """
        return _cantera.wall_area(self.__wall_id)

    def setArea(self, a):
        """
        Set the area (m^2).
        """
        _cantera.wall_setArea(self.__wall_id, a)

    def setThermalResistance(self, rth):
        """Deprecated."""
        return _cantera.wall_setThermalResistance(self.__wall_id, rth)

    def setHeatTransferCoeff(self, u):
        """
        Set the overall heat transfer coefficient [W/m^2/K]
        """
        return _cantera.wall_setHeatTransferCoeff(self.__wall_id, u)

    def setHeatFlux(self, qfunc=None):
        """
        Specify the time-dependent heat flux function [W/m2].
        'qfunc' must be a functor.
        """
        n = 0
        if qfunc: n = qfunc.func_id()
        return _cantera.wall_setHeatFlux(self.__wall_id, n)

    def setExpansionRateCoeff(self, k):
        _cantera.wall_setExpansionRateCoeff(self.__wall_id, k)        
        
    def setExpansionRate(self, vfunc=None):
        """
        Specify the volumetric expansion rate function [m^3/s].
        """
        n = 0
        if vfunc: n = vfunc.func_id()
        _cantera.wall_setExpansionRate(self.__wall_id, n)
            
    def install(self, left, right):
        #self.left = left
        #self.right = right
        left._addWall(this)
        right._addWall(this)
        _cantera.wall_install(self.__wall_id, left.reactor_id(),
                               right.reactor_id())

    def setKinetics(self, left, right):
        """Specify surface reaction mechanisms for the left and right sides of the wall."""
        ileft = 0
        iright = 0
        if left:
            ileft = left.kin_index()
        if right:
            iright = right.kin_index()
        _cantera.wall_setkinetics(self.__wall_id, ileft, iright)
                                  
    def set(self, **p):
        for item in p.keys():
            if item == 'A' or item == 'area':
                self.setArea(p[item])
            elif item == 'R':
                self.setThermalResistance(p[item])
            elif item == 'U':
                self.setHeatTransferCoeff(p[item])                
            elif item == 'K':
                self.setExpansionRateCoeff(p[item])
            elif item == 'Q':
                self.setHeatFlux(p[item])
            elif item == 'Vdot':
                self.setExpansionRate(p[item])
            else:
                raise 'unknown parameter: ',item
                




        

    
    
    
        

    
    
    
        
