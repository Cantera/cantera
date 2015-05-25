.. _sec-python-migration:

Migrating from the Old Python Module
************************************

With the introduction of the new Cython-based Python module in Cantera 2.1,
there are a number of changes to the interface which require modifications to
scripts in order for them to work with the new module. Broadly speaking, the
changes to the interface are intended to make the Cantera Python module easier
to use, and provide a more "Pythonic" interface by making use of common Python
language idioms, language features, and style guidelines.

This document describes the changes to the Python module which are likely to
require modifications to existing code.

Importing the Python Module
---------------------------

The name of the Python module is now ``cantera`` with a lowercase "c". This
change is made partly for compliance with `PEP8
<http://www.python.org/dev/peps/pep-0008/#package-and-module-names>`_.

Furthermore, the various submodules, e.g. ``Cantera.Reactor`` have been
eliminated. All classes and functions are available directly in the
``cantera`` module.

To avoid the namespace clutter introduced by using ``import *``, the following
syntax is preferred::

    >>> import cantera as ct

Naming Conventions
------------------

Generally, the names used in the Cantera Python module have been changed to
follow the recommendations of PEP8. This means that the names of methods and
properties are generally written as ``lowercase_with_underscores`` instead of
``capitalizingEachWord``. Also, some abbreviated names have been expanded. For
example, the following function calls::

    >>> gas.speciesName(0)
    >>> gas.nAtoms('H2', 'H')
    >>> gas.reactionEqn(3)

should be replaced with::

    >>> gas.species_name(0)
    >>> gas.n_atoms('H2', 'H')
    >>> gas.reaction_equation(3)

Importing Phases
----------------

The functions ``importPhase`` and ``IdealGasMix`` have been removed.
`Solution` objects, which represent the phase (regardless of the underlying
thermodynamic model) as well as providing access to kinetics and transport
properties, are created directly using the `Solution` class. For example::

    >>> gas = Solution('h2o2.xml')

Creates an object which represents an ``IdealGasPhase`` mixture with a
``GasKinetics`` reaction mechansm and a ``MixTransport`` transport model,
based on the parameters specified in the input file.

For importing multiple phases from a single file, the ``importPhases`` function
has been retained with the new name ``import_phases``::

    >>> gas, anode_bulk, oxide = ct.import_phases('sofc.cti',
                                                  ['gas', 'metal', 'oxide_bulk'])

Interfaces and edges are created using the `Interface` class, which represents
both 1D and 2D interfaces, rather than using the ``importEdge`` and
``importInterface`` functions::

    >>> anode_surf = ct.Interface('sofc.cti', 'metal_surface', [gas])
    >>> oxide_surf = ct.Interface('sofc.cti', 'oxide_surface', [gas, oxide])
    >>> tpb = ct.Interface('sofc.cti', 'tpb', [anode_bulk, anode_surf, oxide_surf])


Accessing Properties
--------------------

Most methods for accessing and setting the properties of objects have been
replaced with Python "properties" which do not need to be "called" in order to
accessed or changed. For example, the following::

    >>> u = gas.intEnergy_mass()
    >>> Wmx = gas.meanMolecularWeight()
    >>> kf = gas.fwdRateConstants()
    >>> gas.setName('foo')

should be replaced with::

    >>> u = gas.int_energy_mass
    >>> Wmx = gas.mean_molecular_weight
    >>> kf = gas.forward_rate_constants
    >>> gas.name = 'foo'

Some common properties have been renamed according to the variable that is
typically used to represent them::

    >>> gas.temperature()
    >>> gas.pressure()
    >>> gas.massFractions()

should be replaced with::

    >>> gas.T
    >>> gas.P
    >>> gas.Y

For pure fluid phases, the property ``X`` refers to the vapor mass fraction or
"quality" of the phase. The following::

    >>> w = Cantera.liquidvapor.Water()
    >>> w.set(T=400, Vapor=0.5)

should be replaced with::

    >>> w = ct.Water()
    >>> w.TX = 400, 0.5

Setting Thermodynamic State
---------------------------

The ``set`` method has been removed in favor of property pairs or triplets. The
following::

    >>> gas.setMoleFractions('CH4:1.0, O2:0.1')
    >>> gas.set(X='CH4:1.0, O2:0.1')
    >>> gas.set(U=-1.1e6, V=5.5)
    >>> gas.set(T=300, P=101325, Y='H2:1.0')

should be replaced with::

    >>> gas.X = 'CH4:1.0, O2:0.1'
    >>> gas.X = 'CH4:1.0, O2:0.1'
    >>> gas.UV = -1.1e6, 5.5
    >>> gas.TPY = 300, 101325, 'H2:1.0'

The ``saveState`` and ``restoreState`` methods have been removed. Their
functionality can be replicated as follows::

    >>> state = gas.TDY
    >>> # (operations that modify gas)
    >>> gas.TDY = state

Printing Phase Summaries
------------------------

`Solution` objects no longer print out a verbose summary as their string
representation. Instead, the summary report can be generated using the
`report()` method, which returns a string, or by calling the `Solution` object
to print the report to the screen. The following are equivalent::

    >>> print(gas.report())
    >>> gas()

Getting Properties for a Subset of Species
------------------------------------------

Some methods previously accepted an optional list of species as a filter, e.g.::

    >>> gas.massFractions(['OH','H'])

This is not compatible with the Python "property" syntax, so the following
alternative is used instead::

    >>> gas['OH','H2'].Y
    array([ 0.,  1.])

This works for any property which returns a value for each species, and works
with species names, indices, and index ranges::

    >>> gas[1,2,6].partial_molar_cp
    array([ 20786.15525072,  21900.30946418,  34929.99146762])

    >>> gas[3:6].species_names
    ['O2', 'OH', 'H2O']

Furthermore, the "sliced" object itself can be saved and used without needing
to specify the species list again::

    >>> reactants = gas['H2','O2']
    >>> reactants.X
    array([ 1.,  0.])

Transport Models
----------------

The old method for setting the transport model, `switchTransportModel` has been
replaced with the `transport_model` property. To use the multicomponent
transport model::

    >>> gas.transport_model = 'Multi'

Note that unlike the previous implementation, only one transport model can be
associated with a `Solution` object at a time, so there is a larger cost with
switching models. If you need to alternate between transport models, it is
generally better to use two different `Solution` objects.

Reactor Networks
----------------

As with the `Solution` class, properties are now used to get and set most
parameters of reactors, flow devices, walls, etc. The following old code::

    >>> Y = reactor.massFractions()
    >>> X = reactor.contents().moleFractions()
    >>> wall.setArea(2.0)

    >>> net.setTolerances(1e-8, 1e-14)

should be replaced with::

    >>> Y = reactor.Y
    >>> X = reactor.thermo.X
    >>> wall.area = 2.0

    >>> net.rtol = 1e-8
    >>> net.atol = 1e-14

Time-varying parameters have not been replaced with properties, since they
need to be evaluated at a particular time.

Elimination of the ``Func`` Module
----------------------------------

The ``Func`` module is no longer necessary, as the Cython module allows any
callable Python object (lambda, function, or class) to be used in places where
a function of a single variable are needed. For example, to set the velocity
of a wall as a function of time, the following are equivalent::

    >>> wall.set_velocity(lambda t: np.cos(3*t))

    >>> def myfunc(z):
    ...     return np.cos(3*z)
    >>> wall.set_velocity(myfunc)

One-Dimensional Reacting Flows
------------------------------

As elsewhere, the ``set`` method has been eliminated. The following old usage::

    >>> f.fuel_inlet.set(massflux=mdot_f,
    >>>                  mole_fractions=comp_f,
    >>>                  temperature=tin_f)

    >>> f.set(energy = 'off')

should be replaced with::

    >>> f.fuel_inlet.mdot = mdot_f
    >>> f.fuel_inlet.X = comp_f
    >>> f.fuel_inlet.T = tin_f

    >>> f.energy_enabled = False

However, the methods for setting tolerances and refinement criteria have been
retained in slightly modified forms. The following::

    >>> f.set(tol=tol_ss, tol_time=tol_ts)
    >>> f.setRefineCriteria(ratio=4, slope=0.2, curve=0.3, prune=0.04)

should be replaced with::

    >>> f.flame.set_steady_tolerances(default=tol_ss)
    >>> f.flame.set_transient_tolerances(default=tol_ts)
    >>> f.set_refine_criteria(ratio=4, slope=0.2, curve=0.3, prune=0.04)

To change the transport model and enable calculation of the Soret diffusion
term, the following::

    >>> gas.addTransportModel('Multi')
    >>> gas.switchTransportModel('Multi')
    >>> f.flame.setTransportModel(gas)
    >>> f.flame.enableSoret()

should be replaced with::

    >>> f.transport_model = 'Multi'
    >>> f.soret_enabled = True
