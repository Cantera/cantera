"""
Zero-dimensional reactors.
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
        self._type = type
        self._inlets = []
        self._outlets = []
        self._walls = []
        self._reservoirs = []
        self._name = name
        self._verbose = verbose
        self.insert(contents)
        self._setInitialVolume(volume)
        self._setEnergy(energy)
        if self._verbose:
            print 'Created '+self._name
            print '    Volume = ',volume,' m^3'
            if energy <> 'on':
                print '    Temperature will be held constant'
            print '    Initial State:'
            print contents


    def __del__(self):
        """Delete the reactor instance."""
        if self._verbose:
            print 'Deleting '+self._name
        _cantera.reactor_del(self.__reactor_id)

    def __str__(self):
        s = self._name
        s += ':\n     Volume = '+`self.volume()`
        if self._contents:
            s += "\n"+`self._contents`
        return s
    
    def __repr__(self):
        s = self._name
        s += ':\n     Volume = '+`self.volume()`        
        if self._contents:
            s += ": \n"+`self._contents`
        return s
        
    def name(self):
        """The name of the reactor specified when it was constructed."""
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
        self._contents = contents
        if contents:
            _cantera.reactor_setThermoMgr(self.__reactor_id, contents._phase_id)
            _cantera.reactor_setKineticsMgr(self.__reactor_id, contents.ckin)

        
    def setInitialTime(self, T0):
        """Deprecated.
        Set the initial time. Restarts integration from this time
        using the current state as the initial condition. Default: 0.0 s"""
        _cantera.reactor_setInitialTime(self.__reactor_id, T0)

    def _setInitialVolume(self, V0):
        """Set the initial reactor volume. """
        _cantera.reactor_setInitialVolume(self.__reactor_id, V0)

    def _setEnergy(self, eflag):
        """Turn the energy equation on or off. If the argument is the
        string 'off' or the number 0, the energy equation is disabled,
        and the reactor temperature is held constant at its initial
        value."""
        ie = 1
        if eflag == 'off' or eflag == 0:
            ie = 0
        if self._verbose:
            if ie:
                print 'enabling energy equation for reactor',self._name
            else:
                print 'disabling energy equation for reactor',self._name                
        _cantera.reactor_setEnergy(self.__reactor_id, ie)

    def temperature(self):
        """The temperature in the reactor [K]."""
        return _cantera.reactor_temperature(self.__reactor_id)

    def density(self):
        """The density of the fluid in the reactor [kg/m^3]."""
        return _cantera.reactor_density(self.__reactor_id)

    def volume(self):
        """The total reactor volume [m^3]. The volume may change with time
        if non-rigid walls are installed on the reactor."""
        return _cantera.reactor_volume(self.__reactor_id)            

    def time(self):
        """Deprecated. The current time [s]."""
        return _cantera.reactor_time(self.__reactor_id)

    def mass(self):
        """The total mass of fluid in the reactor [kg]."""
        return _cantera.reactor_mass(self.__reactor_id)

    def enthalpy_mass(self):
        """The specific enthalpy of the fluid in the reactor [J/kg]."""
        return _cantera.reactor_enthalpy_mass(self.__reactor_id)

    def intEnergy_mass(self):
        """The specific internal energy of the fluid in the reactor [J/kg]."""
        return _cantera.reactor_intEnergy_mass(self.__reactor_id)

    def pressure(self):
        """The pressure in the reactor [Pa]."""
        return _cantera.reactor_pressure(self.__reactor_id)            

    def advance(self, time):
        """Deprecated.
        Advance the state of the reactor in time from the current
        time to time 'time'. Note: this method is deprecated. See
        class ReactorNet."""
        return _cantera.reactor_advance(self.__reactor_id, time)

    def step(self, time):
        """Deprecated.
        Take one internal time step from the current time toward
        time 'time'. Note: this method is deprecated. See class
        ReactorNet."""
        return _cantera.reactor_step(self.__reactor_id, time)    
    
    def massFraction(self, s):
        """The mass fraction of species s, specified either by name or
        index number.
        >>> y1 = r.massFraction(7)
        >>> y2 = r.massFraction('CH3O')
        """
        if type(s) == types.StringType:
            kk = self._contents.speciesIndex(s)
        else:
            kk = s
        return _cantera.reactor_massFraction(self.__reactor_id, kk)

    def massFractions(self):
        """Return an array of the species mass fractions."""
        nsp = self._contents.nSpecies()
        y = zeros(nsp,'d')
        for k in range(nsp):
            y[k] = self.massFraction(k)
        return y

    def moleFractions(self):
        """Return an array of the species mole fractions."""
        y = self.massFractions()
        self._contents.setMassFractions(y)
        return self._contents.moleFractions()

    def moleFraction(self, s):
        """The mole fraction of species s, specified either by name or
        index number.
        >>> x1 = r.moleFraction(7)
        >>> x2 = r.moleFraction('CH3O')
        """        
        if type(s) == types.StringType:
            kk = self._contents.speciesIndex(s)
        else:
            kk = s
        x = self.moleFractions()
        return x[kk]
        
    def inlets(self):
        """Return the list of flow devices installed on inlets to this reactor.
        This method can be used to access information about the flows entering
        the reactor:
        >>> for n in r.inlets():
        ...    print n.name(), n.massFlowRate()
        See MassFlowController, Valve.
        """        
        return self._inlets

    def outlets(self):
        """Return the list of flow devices installed on outlets
        on this reactor.
        >>> for o in r.outlets():
        ...    print o.name(), o.massFlowRate()
        See MassFlowController, Valve.        
        """
        return self._outlets

    def walls(self):
        """Return the list of walls installed on this reactor.
        >>> for w in r.walls():
        ...    print w.name()
        See Wall.
        """
        return self._walls
    
    def _addInlet(self, inlet, other):
        """For internal use. Store a reference to 'inlet'
        so that it will not be deleted before this object."""        
        self._inlets.append(inlet)
        if self._type == 1 and other._type == 2:
            self._reservoirs.append(other)        

    def _addOutlet(self, outlet, other):
        """For internal use. Store a reference to 'outlet'
        so that it will not be deleted before this object."""        
        self._outlets.append(outlet)
        if self._type == 1 and other._type == 2:
            self._reservoirs.append(other)        

    def _addWall(self, wall, other):
        """For internal use. Store a reference to 'wall'
        so that it will not be deleted before this object."""
        self._walls.append(wall)
        if self._type == 1 and other._type == 2:
            self._reservoirs.append(other)

    def syncContents(self):
        """Set the state of the object representing the reactor contents
        to the current reactor state.
        >>> r = Reactor(gas)
        >>> (statements that change the state of object 'gas')
        >>> r.syncContents()
        After this statement, the state of object 'gas' is synchronized
        with the reactor state.
        See 'contents'.
        """
        self._contents.setState_TRY(self.temperature(),
                                   self.density(),
                                   self.massFractions())
        
    def contents(self):
        """Return an object representing the reactor contents, after first
        synchronizing its state with the current reactor state. This method
        is useful when some property of the fluid in the reactor is
        needed that is not provided by a method of class Reactor.
        >>> r = Reactor(gas)
        >>> (statements that change the state of object 'gas')
        >>> c = r.contents()
        >>> print c.gibbs_mole(), c.chemPotentials()
        
        Note that after calling method 'contents', object 'c'
        references the same underlying kernel object as object 'gas'
        does. Therefore, all properties of 'c' and 'gas' are
        identical. (Remember that Python objects are really C
        pointers; at the C level, both point to the same data
        structure.)
        It is also allowed to write
        >>> gas = r.contents()
        """
        self.syncContents()
        return self._contents
    

_reactorcount = 0
_reservoircount = 0

class Reactor(ReactorBase):
    """
    Zero-dimensional reactors. Instances of class Reactor represent
    zero-dimensional reactors. By default, they are closed (no inlets
    or outlets), have fixed volume, and have adiabatic, chemically-intert
    walls. These properties may all be changed by adding appropriate
    components.
    See classes 'Wall', 'MassFlowController', and 'Valve'.
    """
    def __init__(self, contents = None, name = '',
                 volume = 1.0, energy = 'on',
                 verbose = 0):
        """
        contents - Reactor contents. If not specified, the reactor is
        initially empty. In this case, call method insert to specify
        the contents.

        name - Used only to identify this reactor in output. If not
        specified, defaults to 'Reactor_n', where n is an integer
        assigned in the order Reactor objects are created.

        volume - Initial reactor volume. Defaults to 1 m^3.

        energy - Set to 'on' or 'off'. If set to 'off', the energy
        equation is not solved, and the temperature is held at its
        initial value. The default in 'on'.

        verbose - if set to a non-zero value, additional diagnostic
        information will be printed.

        Some examples showing how to create Reactor objects are shown below.
        >>> gas = GRI30()
        >>> r1 = Reactor(gas)
        This is equivalent to:
        >>> r1 = Reactor()
        >>> r1.insert(gas)
        Arguments may be specified using keywords in any order:
        >>> r2 = Reactor(contents = gas, energy = 'off',
        ...              name = 'isothermal_reactor')
        >>> r3 = Reactor(contents = gas, name = 'adiabatic_reactor')
        Here's an array of reactors:
        >>> reactor_array = [Reactor(), Reactor(gas), Reactor(Air())]
        """
        global _reactorcount
        if name == '':
            name = 'Reactor_'+`_reactorcount`
        _reactorcount += 1
        ReactorBase.__init__(self, contents = contents, name = name,
                             volume = volume, energy = energy,
                             verbose = verbose, type = 1)
            

class Reservoir(ReactorBase):
    """
    A reservoir is a reactor with a constant state. The temperature,
    pressure, and chemical composition in a reservoir never change from
    their initial values.
    """
    def __init__(self, contents = None, name = '', verbose = 0):
        """
        contents - Reservoir contents. If not specified, the reservoir is
        initially empty. In this case, call method insert to specify
        the contents.

        name - Used only to identify this reservoir in output. If not
        specified, defaults to 'Reservoir_n', where n is an integer
        assigned in the order Reservoir objects are created.

        verbose - if set to a non-zero value, additional diagnostic
        information will be printed.

        Some examples showing how to create Reservoir objects are shown below.
        >>> gas = GRI30()
        >>> res1 = Reservoir(gas)
        This is equivalent to:
        >>> res1 = Reactor()
        >>> res1.insert(gas)
        Arguments may be specified using keywords in any order:
        >>> res2 = Reservoir(contents = Air(), 
        ...                  name = 'environment')
        >>> res3 = Reservoir(contents = gas, name = 'upstream_state')
        """        
        global _reservoircount
        if name == '':
            name = 'Reservoir_'+`_reservoircount`
        _reservoircount += 1
        ReactorBase.__init__(self, contents = contents,
                             name = name, verbose = verbose, type = 2)
            
    def advance(self, time):
        """Deprecated. Do nothing."""
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
        """The name specified when initially constructed."""
        return self._name
    
    def ready(self):
        """
        Deprecated. Returns true if the device is ready to use.
        """
        return _cantera.flowdev_ready(self.__fdev_id)

    def massFlowRate(self, time = -999.0):
        """Mass flow rate (kg/s). """
        return _cantera.flowdev_massFlowRate(self.__fdev_id, time)

    def install(self, upstream, downstream):
        """
        Install the device between the upstream and downstream
        reactors or reservoirs.  
        >>> f.install(upstream = reactor1, downstream = reservoir2)
        """
        if self._verbose:
            print
            print self._name+': installing between '+upstream.name()+' and '+downstream.name()
        upstream._addOutlet(self, downstream)
        downstream._addInlet(self, upstream)
        _cantera.flowdev_install(self.__fdev_id, upstream.reactor_id(),
                                  downstream.reactor_id())
    def _setParameters(self, c):
        params = array(c,'d')
        n = len(params)
        return _cantera.flowdev_setParameters(self.__fdev_id, n, params)    

    def setFunction(self, f):
        _cantera.flowdev_setFunction(self.__fdev_id, f.func_id())

    def flowdev_id(self):
        return self.__fdev_id
        
_mfccount = 0

class MassFlowController(FlowDevice):
    
    """Mass flow controllers. A mass flow controller maintains a
    constant mass flow rate independent of upstream and downstream
    conditions. The equation used to compute the mass flow rate is
    \f[ \dot m = \dot m_0, \f] where \f$ \dot m_0 \f$ is a 
    non-negative value specified when the object is constructed or set
    by calling method setMassFlowRate.
    
    Unlike a real mass flow controller, a MassFlowController object
    will maintain the flow even if the downstream pressure is greater
    than the upstream pressure.  This allows simple implementation of
    loops, in which exhaust gas from a reactor is fed back into it
    through an inlet. But note that this capability should be used
    with caution, since no account is taken of the work required to do
    this.

    A mass flow controller is assumed to be adiabatic, non-reactive,
    and have negligible volume, so that it is internally always in
    steady-state even if the upstream and downstream reactors are
    not. The fluid enthalpy, chemical composition, and mass flow rate
    are constant across a mass flow controller, and the pressure
    difference equals the difference in pressure between the upstream
    and downstream reactors.
    """
    def __init__(self, upstream=None,
                 downstream=None,
                 name='',
                 verbose=0, mdot = 0.0):
        """
        upstream - upstream reactor or reservoir.

        downstream - downstream reactor or reservoir.

        name - name used to identify the mass flow controller in output.
        If no name is specified, it defaults to 'MFC_n', where n is an
        integer assigned in the order the MassFlowController object
        was created.

        mdot - Mass flow rate [kg/s]. This mass flow rate will be
        maintained, independent of unstream and downstream conditions,
        unless reset by calling method 'setMassFlowRate'.

        verbose - if set to a positive integer, additional diagnostic
        information will be printed.
        
        """
        global _mfccount
        if name == '':
            name = 'MFC_'+`_mfccount`
        _mfccount += 1
        FlowDevice.__init__(self,1,name,verbose)
        if upstream and downstream:
            self.install(upstream, downstream)
        if mdot:
            self.set(mdot = mdot)

    def _setMassFlowRate(self, mdot):
        """Set or reset the mass flow rate to 'mdot' [kg/s].
        """
        if self._verbose:
            print self._name+': setting mdot to '+`mdot`+' kg/s'
        if type(mdot) == types.InstanceType:
            self.setFunction(mdot)
        else:
            _cantera.flowdev_setMassFlowRate(self.flowdev_id(), mdot)


    def set(self, mdot = 0.0):
        """Set the mass flow rate [kg/s].

        >>> mfc.set(mdot = 0.2)
        """
        self._setMassFlowRate(mdot)


_valvecount = 0

class Valve(FlowDevice):
    """Valves. In Cantera, a Valve object is a flow devices with mass
    flow rate proportional to the pressure drop across it. The equation
    used to compute the mass flow rate is
    \f[ \dot m = K_v (P_1 - P_2) \f]
    if \f$ P_1 > P_2. \f$
    Otherwise,
    \f$ \dot m = 0 \f$. It is never possible for the flow to reverse
    and go from the downstream to the upstream reactor/reservoir through
    a line containing a Valve object. 

    'Valve' objects are often used between an upstream reactor and a
    downstream reactor or reservoir to maintain them both at nearly the
    same pressure. By setting the constant \f$ K_v \f$ to a
    sufficiently large value, very small pressure differences will
    result in flow between the reactors that counteracts the pressure
    difference.

    Since the mass flow rate is assumed to be linear in \f$ \Delta P \f$,
    these objects do not model real, physical valves, in which the flow rate 
    is proportional to \f$ \sqrt(\Delta P) \f$ for small pressure
    differences, and becomes independent of \f$ \Delta P \f$ when
    it becomes large (choked flow). Perhaps the name of this class should
    be changed to avoid confusion with real valves -- if you have suggestions,
    post a comment at  the Cantera User's Group site.

    A Valve is assumed to be adiabatic, non-reactive, and have
    negligible internal volume, so that it is internally always in
    steady-state even if the upstream and downstream reactors are
    not. The fluid enthalpy, chemical composition, and mass flow rate
    are constant across a Valve, and the pressure difference equals
    the difference in pressure between the upstream and downstream
    reactors.
    
    """
    def __init__(self, upstream=None, downstream=None,
                 name='', Kv = 0.0, mdot0 = 0.0, verbose=0):
        """
        upstream - upstream reactor or reservoir.
        
        downstream - downstream reactor or reservoir.
        
        name - name used to identify the valve in output.
        If no name is specified, it defaults to 'Valve_n', where n is an
        integer assigned in the order the Valve object
        was created.

        Kv - the constant in the mass flow rate equation.

        verbose - if set to a positive integer, additional diagnostic
        information will be printed.
        
        """        
        global _valvecount
        if name == '':
            name = 'Valve_'+`_valvecount`
        _valvecount += 1
        FlowDevice.__init__(self,3,name,verbose)
        if upstream and downstream:
            self.install(upstream, downstream)
        self.setValveCoeff(Kv, mdot0)


    def setValveCoeff(self, Kv = -1.0, mdot0 = 0.0):
        """Set or reset the valve coefficient \f$ K_v \f$."""
        vv = zeros(2,'d')
        vv[0] = Kv
        vv[1] = mdot0
        if self._verbose:
            print
            print self._name+': setting valve coefficient to '+`Kv`+' kg/Pa-s'
        self._setParameters(vv)

    def _setValveCharacteristic(self, f):
        """Set or reset the valve characteristics.
        """
        if type(f) == types.InstanceType:
            self.setFunction(f)
        else:
            raise CanteraError("Wrong type for valve characteristic function.")

    def set(self, Kv = -1.0, mdot = 0.0, F = None):
        if F:
            self.setFunction(F)
        if Kv > 0.0:
            self.setValveCoeff(Kv, mdot0 = mdot)
            


_pccount = 0

class PressureController(FlowDevice):

    def __init__(self, upstream=None, downstream=None,
                 name='', master = None, Kv = 0.0, verbose=0):
        """
        upstream - upstream reactor or reservoir.
        
        downstream - downstream reactor or reservoir.
        
        name - name used to identify the valve in output.
        If no name is specified, it defaults to 'Valve_n', where n is an
        integer assigned in the order the Valve object
        was created.

        Kv - the constant in the mass flow rate equation.

        verbose - if set to a positive integer, additional diagnostic
        information will be printed.
        
        """        
        global _pccount
        if name == '':
            name = 'PressureController_'+`_pccount`
        _pccount += 1
        FlowDevice.__init__(self,2,name,verbose)
        if upstream and downstream:
            self.install(upstream, downstream)
        self.setPressureCoeff(Kv)
        self.setMaster(master)


    def setPressureCoeff(self, Kv):
        """Set or reset the pressure coefficient \f$ K_v \f$."""
        vv = zeros(1,'d')
        vv[0] = Kv
        if self._verbose:
            print
            print self._name+': setting pressure coefficient to '+`Kv`+' kg/Pa-s'
        self._setParameters(vv)

    def setMaster(self, master):
        _cantera.flowdev_setMaster(self.flowdev_id(),
                                   master.flowdev_id())
            
    def set(self, Kv = -1.0, master = None):
        if master:
            self.setMaster(master)
        if Kv > 0.0:
            self.setPressureCoeff(Kv)
            

            
#------------- Wall ---------------------------

_wallcount = 0

class Wall:
    """
    Reactor walls.
    A Wall separates two reactors, or a reactor and a reservoir.
    """
    def __init__(self, left=None, right=None, name = '',
                 A = 1.0, K = 0.0, U = 0.0,
                 Q = None, velocity = None,
                 kinetics = [None, None]):
        typ = 0
        self.__wall_id = _cantera.wall_new(typ)

        global _wallcount
        if name == '':
            _nm = 'Wall_'+`_wallcount`
        else:
            _nm = name
        _wallcount += 1
        
        if left and right:
            self.install(left, right)
        elif left or right:
            raise CanteraError('both left and right reactors must be specified.')
        self.setArea(A)
        self.setExpansionRateCoeff(K)
        self.setVelocity(velocity)        
        self.setHeatTransferCoeff(U)
        self.setHeatFlux(Q)

        self.setKinetics(kinetics[0],kinetics[1])

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

    def setEmissivity(self, epsilon):
        """
        Set the emissivity.
        The radiative heat flux through the wall is computed from
        \f[ q_r = \epsion \sigma (T_\ell^4 - T_r^4) \f]
        """
        
    def setHeatFlux(self, qfunc=None):
        """
        Specify the time-dependent heat flux function [W/m2].
        'qfunc' must be a functor.
        """
        n = 0
        if qfunc: n = qfunc.func_id()
        return _cantera.wall_setHeatFlux(self.__wall_id, n)

    def setExpansionRateCoeff(self, k):
        """Set the coefficient K that determines the expansion rate
        resulting from a unit pressure drop."""
        _cantera.wall_setExpansionRateCoeff(self.__wall_id, k)        
        
    def setVelocity(self, vfunc=None):
        """
        Specify the velocity function [m/s].
        """
        n = 0
        if vfunc: n = vfunc.func_id()
        _cantera.wall_setVelocity(self.__wall_id, n)

    def vdot(self):
        """Rate of volume change [m^3]. A positive value corresponds
        to the left-hand reactor volume increasing, and the right-hand
        reactor volume decreasing."""
        return _cantera.wall_vdot(self.__wall_id)

    def heatFlowRate(self):
        """Rate of heat flow through the wall. A positive value
        corresponds to heat flowing from the left-hand reactor to the
        right-hand one."""
        return _cantera.wall_Q(self.__wall_id)
    
    def install(self, left, right):
        left._addWall(self, right)
        right._addWall(self, left)
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
                

class ReactorNet:
    
    """Networks of reactors. ReactorNet objects are used to
    simultaneously advance the state of a set of coupled reactors.

    Example:

    >>> r1 = Reactor(gas1)
    >>> r2 = Reactor(gas2)
    >>> <... install walls, inlets, outlets, etc...>

    >>> reactor_network = ReactorNet([r1, r2])
    >>> reactor_network.advance(time)
    
    """


    def __init__(self, reactorlist = None):
        """
        Create a new ReactorNet instance. If a list of reactors is supplied,
        these will be added to the network.
        """
        self._reactors = []
        self.__reactornet_id = _cantera.reactornet_new()
        if reactorlist:
            for r in reactorlist:
                self.add(r)


    def __del__(self):
        """Delete the reactor network instance. The reactors in the
        network are not deleted."""
        _cantera.reactornet_del(self.__reactornet_id)


    def reactornet_id(self):
        """ The integer index used to access the
        kernel reactornet object. For internal use.  """
        return self.__reactornet_id

    
    def add(self, reactor):
        """
        Add a reactor to the network.
        """
        self._reactors.append(reactor)
        _cantera.reactornet_addreactor(self.__reactornet_id,
                                       reactor.reactor_id())

        
    def setInitialTime(self, t0):
        """Set the initial time. Restarts integration from this time
        using the current state as the initial condition. Default: 0.0 s"""
        _cantera.reactornet_setInitialTime(self.__reactornet_id, t0)

    def advance(self, time):
        """Advance the state of the reactor network in time from the current
        time to time 'time'."""
        return _cantera.reactornet_advance(self.__reactornet_id, time)

    def step(self, time):
        """Take a single internal time step toward time 'time'.
        The time after taking the step is returned."""
        return _cantera.reactornet_step(self.__reactornet_id, time)    

    def reactors(self):
        return self._reactors
    
