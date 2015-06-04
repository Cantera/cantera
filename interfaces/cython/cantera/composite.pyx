class Solution(ThermoPhase, Kinetics, Transport):
    """
    A class for chemically-reacting solutions. Instances can be created to
    represent any type of solution -- a mixture of gases, a liquid solution, or
    a solid solution, for example.

    Class `Solution` derives from classes `ThermoPhase`, `Kinetics`, and
    `Transport`.  It defines no methods of its own, and is provided so that a
    single object can be used to compute thermodynamic, kinetic, and transport
    properties of a solution.

    To skip initialization of the Transport object, pass the keyword argument
    ``transport_model=None`` to the `Solution` constructor.

    The most common way to instantiate `Solution` objects is by using a phase
    definition, species and reactions defined in an input file::

        gas = ct.Solution('gri30.cti')

    If an input file defines multiple phases, the phase *name* (in CTI) or *id*
    (in XML) can be used to specify the desired phase::

        gas = ct.Solution('diamond.cti', 'gas')
        diamond = ct.Solution('diamond.cti', 'diamond')

    `Solution` objects can also be constructed using `Species` and `Reaction`
    objects which can themselves either be imported from input files or defined
    directly in Python::

        spec = ct.Species.listFromFile('gri30.cti')
        rxns = ct.Reaction.listFromFile('gri30.cti')
        gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                          species=spec, reactions=rxns)

    where the ``thermo`` and ``kinetics`` keyword arguments are strings
    specifying the thermodynamic and kinetics model, respectively, and
    ``species`` and ``reactions`` keyword arguments are lists of `Species` and
    `Reaction` objects, respectively.

    For non-trivial uses cases of this functionality, see the examples
    :ref:`py-example-extract_submechanism.py` and
    :ref:`py-example-mechanism_reduction.py`.

    In addition, `Solution` objects can be constructed by passing the text of
    the CTI or XML phase definition in directly, using the ``source`` keyword
    argument::

        cti_def = '''
            ideal_gas(name='gas', elements='O H Ar',
                      species='gri30: all',
                      reactions='gri30: all',
                      options=['skip_undeclared_elements', 'skip_undeclared_species', 'skip_undeclared_third_bodies'],
                      initial_state=state(temperature=300, pressure=101325))'''
        gas = ct.Solution(source=cti_def)
    """
    __slots__ = ()

class Interface(InterfacePhase, InterfaceKinetics):
    """
    Two-dimensional interfaces.

    Instances of class `Interface` represent reacting 2D interfaces between bulk
    3D phases. Class `Interface` defines no methods of its own. All of its
    methods derive from either `InterfacePhase` or `InterfaceKinetics`.

    To construct an `Interface` object, adjacent bulk phases which participate
    in reactions need to be created and then passed in as a list in the
    ``phases`` argument to the constructor::

        gas = ct.Solution('diamond.cti', 'gas')
        diamond = ct.Solution('diamond.cti', 'diamond')
        diamond_surf = ct.Interface('diamond.cti', 'diamond_100', [gas, diamond])
    """
    __slots__ = ('_phase_indices',)

class DustyGas(ThermoPhase, Kinetics, DustyGasTransport):
    """
    A composite class which models a gas in a stationary, solid, porous medium.

    The only transport properties computed are the multicomponent diffusion
    coefficients. The model does not compute viscosity or thermal conductivity.

    """
    __slots__ = ()
