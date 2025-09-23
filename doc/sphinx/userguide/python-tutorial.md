---
file_format: mystnb
kernelspec:
  name: python3
---

```{py:currentmodule} cantera
```

# Python Tutorial

## Getting Started

```{seealso}
If you haven't already installed Cantera, we have instructions for
many platforms in our [Installation section](/install/index).
```

Start by opening an interactive Python session, for example by running
[IPython](http://ipython.org/). Import the Cantera Python module and NumPy by running:

```{code-cell} python
import cantera as ct
import numpy as np
```

When using Cantera, the first thing you usually need is an object representing: some
phase of matter. Here, we'll create a gas mixture

```{code-cell} python
gas1 = ct.Solution('gri30.yaml')
```

To view the state of the mixture, *call* the `gas1` object as if it were
a function:

```{code-cell}
gas1()
```

What you have just done is to create an object, `gas1` that implements GRI-Mech 3.0, the
53-species, 325-reaction natural gas combustion mechanism developed by
{cite:t}`smith1999`.

The `gas1` object has properties you would expect for a gas mixture - it has a
temperature, a pressure, species mole and mass fractions, etc. As we'll soon see, it has
many more properties.

The summary of the state of `gas1` printed above shows that new objects created from the
`gri30.yaml` input file start out with a temperature of 300 K, a pressure of 1 atm, and
have a composition that consists of only one species, in this case hydrogen. There is
nothing special about H2 - it just happens to be the first species listed in the input
file defining GRI-Mech 3.0. In general, whichever species is listed first will initially
have a mole fraction of 1.0, and all of the others will be zero.

## Setting the State

The state of the object can easily be changed. For example:

```{code-cell}
gas1.TP = 1200, 101325
```

sets the temperature to 1200 K and the pressure to 101325 Pa (Cantera always uses SI
units). After this statement, calling `gas1()` results in:

```{code-cell}
gas1()
```

Thermodynamics generally requires that *two* properties in addition to composition
information be specified to fix the intensive state of a substance (or mixture). The
state of the mixture can be set using several combinations of two properties. The
following are all equivalent:

```{code-cell} pythonconsole
gas1.TP = 1200, 101325          # temperature, pressure
gas1.TD = 1200, 0.020473        # temperature, density
gas1.HP = 1.3295e7, 101325      # specific enthalpy, pressure
gas1.UV = 8.3457e6, 1/0.020473  # specific internal energy, specific volume
gas1.SP = 85222, 101325         # specific entropy, pressure
gas1.SV = 85222, 1/0.020473     # specific entropy, specific volume
```

In each case, the values of the extensive properties must be entered *per unit mass*.

Properties may be read independently or together:

```{code-cell}
gas1.T
```

```{code-cell}
gas1.h
```

```{code-cell}
gas1.UV
```

The composition can be set in terms of either mole fractions (`X`) or mass fractions
(`Y`):

```{code-cell}
gas1.X = 'CH4:1, O2:2, N2:7.52'
```

Mass and mole fractions can also be set using `dict` objects, which is convenient in
cases where the composition is stored in a variable or being computed:

```{code-cell} python
phi = 0.8
gas1.X = {'CH4':1, 'O2':2/phi, 'N2':2*3.76/phi}
```

When the composition alone is changed, the temperature and density are held constant.
This means that the pressure and other intensive properties will change. The composition
can also be set in conjunction with the intensive properties of the mixture:

```{code-cell} python
gas1.TPX = 1200, 101325, 'CH4:1, O2:2, N2:7.52'
gas1()
```

The composition above was specified using a string. The format is a comma-separated list
of `<species name>:<relative mole numbers>` pairs. The mole numbers will be normalized
to produce the mole fractions, and therefore they are "relative" mole numbers. Mass
fractions can be set in this way too by changing `X` to `Y` in the above statements.

The composition can also be set using an array, which must have the same size as the
number of species. For example, to set all 53 mole fractions to the same value, do this:

```{code-cell} python
gas1.X = np.ones(53)  # NumPy array of 53 ones
```

Or, to set all the mass fractions to equal values:

```{code-cell}
gas1.Y = np.ones(53)
```

When setting the state, you can control what properties are held constant by passing the
special value `None` to the property setter. For example, to change the specific volume
to 2.1 m$^3$/kg while holding entropy constant:

```{code-cell}
gas1.SV = None, 2.1
```

Or to set the mass fractions while holding temperature and pressure constant:

```{code-cell}
gas1.TPX = None, None, 'CH4:1.0, O2:0.5'
```

## Working with a Subset of Species

Many properties of a {py:class}`Solution` provide values for each species
present in the phase. If you want to get values only for a subset of these species, you
can use Python's "slicing" syntax to select data for just the species of interest. To
get the mole fractions of just the major species in `gas1`, in the order specified, you
can write:

```{code-cell} python
Xmajor = gas1['CH4','O2','CO2','H2O','N2'].X
```

If you want to use the same set of species repeatedly, you can keep a reference to the
sliced phase object:

```{code-cell} python
major = gas1['CH4','O2','CO2','H2O','N2']
cp_major = major.partial_molar_cp
wdot_major = major.net_production_rates
```

The slice object and the original object share the same internal state, so modifications
to one will affect the other.

## Working With Mechanism Files

In previous example, we created an object that models an ideal gas mixture with the
species and reactions of GRI-Mech 3.0, using the `gri30.yaml` input file included with
Cantera. Several other reaction mechanism files are included with Cantera, including
ones that model high- temperature air, a hydrogen/oxygen reaction mechanism, and a few
surface reaction mechanisms. These files are usually located in the `data` subdirectory
of the Cantera installation directory, for example `C:\Program Files\Cantera\data` on
Windows or `/usr/local/cantera/data/` on Unix/Linux/Mac OS X machines, depending on how
you installed Cantera and the options you specified.

If for some reason Cantera has difficulty finding where these files are on your system,
set environment variable `CANTERA_DATA` to the directory or directories (separated using
`;` on Windows or `:` on other operating systems) where they are located. Alternatively,
you can call function {py:func}`add_data_directory` to add a directory to the
Cantera search path:

```{code-cell} python
ct.add_data_directory('~/cantera/my_data_files')
```

Cantera input files are plain text files, and can be created with any text editor. See
the page [](input-tutorial) for more information.

A Cantera input file may contain more than one phase specification, and may contain
specifications of interfaces (surfaces and edges). Here we import a surface phase and
the two adjacent bulk phases from the file `diamond.yaml`:

```{code-cell} python
diamond_surf = ct.Interface('diamond.yaml' , 'diamond_100')
diamond_surf.adjacent
```

## Converting CK-format files

See the page [](ck2yaml-tutorial) for information on how to convert from CK-format to
Cantera's YAML format.

## Getting Help

In addition to the Sphinx-generated [Python documentation](/python/index), documentation
of the Python classes and their methods can be accessed from within the Python
interpreter as well.

Suppose you have created a Cantera object and want to know what methods are available
for it, and get help on using the methods:

```{code-cell}
g = ct.Solution('gri30.yaml')
```

To get help on the Python class that this object is an instance of:

```{code-cell}
:tags: [hide-output]
help(g)
```

For a simple list of the properties and methods of this object:

```{code-cell}
:tags: [hide-output]
dir(g)
```

To get help on a specific method, such as the `species_index` method:

```{code-cell}
help(g.species_index)
```

For properties, getting the documentation is slightly trickier, as the usual method will
give you the help for the *result*. For example:

```{code-cell}
:tags: [hide-output]
help(g.T)
```

will provide help on Python's `float` class. To get the help for the temperature
property, ask for the attribute of the class object itself:

```{code-cell}
help(g.__class__.T)
```

If you are using the IPython shell, help can also be obtained using the `?` syntax:

```{code-cell}
In[1]: g.species_index?
```

## Chemical Equilibrium

To set a gas mixture to a state of chemical equilibrium, use the
{py:meth}`ThermoPhase.equilibrate` method:

```{code-cell} python
import cantera as ct
g = ct.Solution('gri30.yaml')
g.TPX = 300.0, ct.one_atm, 'CH4:0.95,O2:2,N2:7.52'
g.equilibrate('TP')
```

The above statement sets the state of object `g` to the state of chemical equilibrium
holding temperature and pressure fixed. Alternatively, the specific enthalpy and
pressure can be held fixed:

```{code-cell} python
g.TPX = 300.0, ct.one_atm, 'CH4:0.95,O2:2,N2:7.52'
g.equilibrate('HP')
```

Other options are:

-   `UV`: fixed specific internal energy and specific volume
-   `SV`: fixed specific entropy and specific volume
-   `SP`: fixed specific entropy and pressure

How can you tell if `equilibrate` has correctly found the chemical equilibrium state?
One way is verify that the net rates of progress of all reversible reactions are zero.
Here is the code to do this:

```{code-cell} python
:tags: [hide-output]
g.TPX = 300.0, ct.one_atm, 'CH4:0.95,O2:2,N2:7.52'
g.equilibrate('HP')

rf = g.forward_rates_of_progress
rr = g.reverse_rates_of_progress
for i in range(g.n_reactions):
    if g.reaction(i).reversible and rf[i] != 0.0:
        print(' %4i  %10.4g  ' % (i, (rf[i] - rr[i])/rf[i]))
```

If the magnitudes of the numbers in this list are all very small, then each reversible
reaction is very nearly equilibrated, which only occurs if the gas is in chemical
equilibrium.

You might be wondering how `equilibrate` works. (Then again, you might not). Method
`equilibrate` invokes Cantera's chemical equilibrium solver, which uses an element
potential method. The element potential method is one of a class of equivalent
*nonstoichiometric* methods that all have the characteristic that the problem reduces to
solving a set of $M$ nonlinear algebraic equations, where $M$ is the number of elements
(not species). The so-called *stoichiometric* methods, on the other hand, (including
Gibbs minimization), require solving $K$ nonlinear equations, where $K$ is the number of
species (usually $K >> M$). See {cite:t}`smith1982` for more information on the various
algorithms and their characteristics.

Cantera uses a damped Newton method to solve these equations, and does a few other
things to generate a good starting guess and to produce a reasonably robust algorithm.
If you want to know more about the details, look at the C++ class {ct}`ChemEquil`.

## Chemical Kinetics

{py:class}`Solution` objects are also {py:class}`Kinetics` objects, and provide all of
the methods necessary to compute the thermodynamic quantities associated with each
reaction, reaction rates, and species creation and destruction rates. They also provide
methods to inspect the quantities that define each reaction such as the rate constants
and the stoichiometric coefficients. The rate calculation functions are used extensively
within Cantera's [reactor network model](../python/zerodim) and [1D flame
model](../python/onedim).

Information about individual reactions that is independent of the thermodynamic state
can be obtained by accessing {py:class}`Reaction` objects with the
{py:meth}`Kinetics.reaction` method:

```{code-cell} python
g = ct.Solution('gri30.yaml')
r = g.reaction(2) # get a Reaction object
r
```

```{code-cell}
r.reactants
```

```{code-cell}
r.products
```

Information about specific reaction rate parameterizations should be queried using the
`input_data` property, which returns a YAML-compatible dictionary that represents input
data needed to create the corresponding rate object:

```{code-cell}
r.rate
```

```{code-cell}
r.rate.input_data
```

If we are interested in only certain types of reactions, we can use this information to
filter the full list of reactions to find the just the ones of interest. For example,
here we find the indices of just those reactions which convert `CO` into `CO2`:

```{code-cell} python
II = [i for i,r in enumerate(g.reactions())
      if 'CO' in r.reactants and 'CO2' in r.products]
for i in II:
    print(g.reaction(i).equation)
```

(Actually, we should also include reactions where the reaction is written such that
`CO2` is a reactant and `CO` is a product, but for this example, we'll just stick to
this smaller set of reactions.) Now, let's set the composition to an interesting
equilibrium state:

```{code-cell}
g.TPX = 300, 101325, {'CH4':0.6, 'O2':1.0, 'N2':3.76}
g.equilibrate('HP')
```

We can verify that this is an equilibrium state by seeing that the net reaction rates
are essentially zero:

```{code-cell}
g.net_rates_of_progress[II]
```

Now, let's see what happens if we decrease the temperature of the mixture:

```{code-cell} python
g.TP = g.T-100, None
g.net_rates_of_progress[II]
```

All of the reaction rates are positive, favoring the formation of `CO2` from `CO`, with
the third reaction, `CO + OH <=> CO2 + H` proceeding the fastest. If we look at the
enthalpy change associated with each of these reactions:

```{code-cell}
g.delta_enthalpy[II]
```

we see that the change is negative in each case, indicating a net release of thermal
energy. The total heat release rate can be computed either from the reaction rates:

```{code-cell}
np.dot(g.net_rates_of_progress, g.delta_enthalpy)
```

or from the species production rates:

```{code-cell}
np.dot(g.net_production_rates, g.partial_molar_enthalpies)
```

The contribution from just the selected reactions is:

```{code-cell}
np.dot(g.net_rates_of_progress[II], g.delta_enthalpy[II])
```

Or about 16% of the total heat release rate.

## Next Steps

Congratulations! You have finished the Cantera Python tutorial. You should now be ready
to begin using Cantera on your own. Please see additional sections of the [Cantera User
Guide](/userguide/index) for assistance with intermediate and advanced Cantera
functionality. You can also take a look at our [Example Gallery](/examples/python/index)
for ideas on how you can use Cantera. Here are a few examples that may be of particular
interest:

* [](/userguide/flame-temperature)
* [](/userguide/heating-value)
* [](/examples/python/thermo/equations_of_state)
* [](/examples/python/kinetics/lithium_ion_battery)
* [](/examples/python/reactors/continuous_reactor)
* [](/examples/python/reactors/non_ideal_shock_tube)
* [](/examples/python/onedim/flamespeed_sensitivity)
* [](/examples/python/onedim/premixed_counterflow_twin_flame)
