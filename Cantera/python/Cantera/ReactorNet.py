"""
Reactor networks.

"""

import _cantera

class ReactorNet:
    
    """Networks of reactors. ReactorNet objects are used to
    simultaneously advance the state of a set of coupled reactors.

    Example:

    r1 = Reactor(gas1)
    r2 = Reactor(gas2)
    <... install walls, inlets, outlets, etc...>

    reactor_network = ReactorNet([r1, r2])
    reactor_network.advance(time)
    
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
        _cantera.reactornet_addreactor(self.__reactornet_id, reactor.reactor_id())

        
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
    
