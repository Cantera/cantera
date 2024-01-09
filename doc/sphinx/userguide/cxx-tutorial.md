# Using Cantera from C++

One of the main reasons to use the Cantera C++ interface directly is to be able to
couple it with a higher-level simulation code, where Cantera is used to calculate
thermodynamic properties, reaction rates, and transport properties for the host code.
This tutorial demonstrates how to initialize Cantera objects from
[YAML input files](/userguide/input-tutorial.html) and use those objects to calculate
properties for mixtures with different thermodynamic states.

## A very simple C++ program

A short C++ program that uses Cantera is shown below. This program reads in a
specification of a gas mixture from an input file, and then builds a new object
representing the mixture. It then sets the thermodynamic state and composition of the
gas mixture, and prints out a summary of its properties.

```{literalinclude} demo1a.cpp
:language: c++
```

:::{tip}
Before you can run this program, it first needs to be compiled. On a Linux system using
the GCC compiler, a typical command line for compiling this program might look like
this:

```bash
g++ combustor.cpp -o combustor -O3 $(pkg-config --cflags --libs cantera)
```

This example relies on the `pkg-config` tool to determine the appropriate compiler
flags, such as those specifying the Cantera header and library files. For more advanced
and flexible methods of compiling programs that use the Cantera C++ library, see
[](compiling-cxx).
:::

Running this program produces the output below:

```none
 ohmech:

      temperature   500 K
         pressure   2.0265e+05 Pa
          density   0.36118 kg/m^3
 mean mol. weight   7.4093 kg/kmol
  phase of matter   gas

                         1 kg             1 kmol
                    ---------------   ---------------
         enthalpy       -2.4772e+06       -1.8354e+07  J
  internal energy       -3.0382e+06       -2.2511e+07  J
          entropy             20699        1.5337e+05  J/K
   Gibbs function       -1.2827e+07       -9.5038e+07  J
heat capacity c_p            3919.1             29038  J/K
heat capacity c_v              2797             20724  J/K

                     mass frac. Y      mole frac. X     chem. pot. / RT
                    ---------------   ---------------   ---------------
               H2           0.21767               0.8           -15.644
                H                 0                 0
                O                 0                 0
               O2                 0                 0
               OH                 0                 0
              H2O           0.24314               0.1           -82.953
              HO2                 0                 0
             H2O2                 0                 0
               AR           0.53919               0.1           -20.503
               N2                 0                 0
```

As C++ programs go, this one is quite short. It is the Cantera equivalent of the "Hello,
World" program most programming textbooks begin with. But it illustrates some important
points in writing Cantera C++ programs.

### Making Cantera classes and functions available

Cantera provides several header files designed for use in C++ application programs.
These are designed to include those portions of Cantera needed for particular types of
calculations. The headers and their functions are:

`core.h`
: Base classes and functions for creating {ct}`Solution` and {ct}`Interface` objects
  from input files, as well as associated classes accessed through these objects such as
  {ct}`ThermoPhase`, {ct}`Kinetics`, and {ct}`Transport`. *(New in Cantera 3.0)*.

`zerodim.h`
: Zero-dimensional reactor networks.

`onedim.h`
: One-dimensional reacting flows.

`reactionpaths.h`
: Reaction path diagrams.

`thermo.h`
: Base thermodynamic classes and functions for creating {ct}`ThermoPhase` objects from
  input files *(Superseded by `core.h` in Cantera 3.0)*.

`kinetics.h`
: Base kinetics classes and functions for creating {ct}`Kinetics` objects from input
  files *(Superseded by `core.h` in Cantera 3.0)*.

`transport.h`
: Base transport property classes and functions for creating {ct}`Transport` objects
  from input files *(Superseded by `core.h` in Cantera 3.0)*.

### Creating a `Solution` object from an input file

Here, {ct}`newSolution` imports all information held by a YAML input file into a Cantera
{ct}`Solution` object, which is accessed by the pointer `sol`. The thermodynamic
information is accessible via `sol->thermo()`, which itself returns a pointer to a
{ct}`ThermoPhase` object.

Class {ct}`ThermoPhase` is the base class for Cantera classes that represent phases of
matter. It defines the public interface for all classes that represent phases. For
example, it specifies that they all have a method `temperature()` that returns the current
temperature, a method `setTemperature(double T)` that sets the temperature, a method
`getChemPotentials(double* mu)` that writes the species chemical potentials into array
`mu`, and so on.

{ct}`ThermoPhase` objects can represent the intensive state of any single-phase solution
of multiple species. The phase may be a bulk, three-dimensional phase (a gas, a liquid,
or a solid), or it may be a two-dimensional surface phase, or even a one-dimensional
"edge" phase. Many different [phase thermodynamic models](/reference/thermo/phase-thermo)
are implemented as classes derived from {ct}`ThermoPhase`.

### The `report` function

The {ct}`ThermoPhase::report` method generates a nicely-formatted report of the
properties of a phase, including its composition in both mole (X) and mass (Y) units.
For each species present, the non-dimensional chemical potential is also printed. This
is handy particularly when doing equilibrium calculations. This function is very useful
to see at a glance the state of some phase.

### Handling errors

The entire body of the program is put inside a function that is invoked within a `try`
block in the main program. In this way, exceptions thrown in the function or in any
procedure it calls may be caught. In this program, a `catch` block is defined for
exceptions of type {ct}`CanteraError`. Cantera throws exceptions of this type, so it is
always a good idea to catch them.

## Accessing thermodynamic properties

In the program below, a gas mixture object is created, and a few thermodynamic
properties are computed and printed out:

```{literalinclude} thermodemo.cpp
:language: c++
```

Methods that compute scalar properties take no input parameters. The properties are
computed for the state that has been previously set and stored internally within the
object. Methods that return *molar* properties have names that end in `_mole`, while
methods that return properties *per unit mass* have names that end in `_mass`.

Methods that compute an array of values have names beginning with `get` and write the
computed values into an array provided as an input argument. For example, the method
`getChemPotentials(double* mu)` writes the species chemical potentials into the output
array `mu`.

All thermodynamic property methods are declared in class {ct}`ThermoPhase`, which is the
base class from which all classes that represent any type of phase of matter derive. See
the documentation for class {ct}`ThermoPhase` for the full list of available
thermodynamic properties.

## Computing chemical equilibrium

In the program below, the {ct}`ThermoPhase::equilibrate` method is called to set a gas
to a state of chemical equilibrium, holding the temperature and pressure fixed.

```{literalinclude} demoequil.cpp
:language: c++
```

The program output is:

```none
 ohmech:

      temperature   1500 K
         pressure   2.0265e+05 Pa
          density   0.31683 kg/m^3
 mean mol. weight   19.499 kg/kmol
  phase of matter   gas

                         1 kg             1 kmol
                    ---------------   ---------------
         enthalpy       -4.1789e+06       -8.1485e+07  J
  internal energy       -4.8186e+06       -9.3957e+07  J
          entropy             11283        2.2001e+05  J/K
   Gibbs function       -2.1104e+07        -4.115e+08  J
heat capacity c_p              1893             36912  J/K
heat capacity c_v            1466.6             28597  J/K

                     mass frac. Y      mole frac. X     chem. pot. / RT
                    ---------------   ---------------   ---------------
               H2          0.025847              0.25           -19.295
                H        3.2181e-07        6.2252e-06           -9.6477
                O        6.2927e-12        7.6693e-12           -26.377
               O2        1.1747e-11        7.1586e-12           -52.753
               OH        3.0994e-07        3.5535e-07           -36.024
              H2O           0.46195               0.5           -45.672
              HO2        1.2362e-14        7.3034e-15           -62.401
             H2O2         6.904e-13        3.9578e-13           -72.049
               AR           0.51221              0.25           -21.339
               N2                 0                 0
```

How can we tell that this is really a state of chemical equilibrium? Well, by applying
the equation of reaction equilibrium to formation reactions from the elements, it is
straightforward to show that

$$  \mu_k = \sum_m \lambda_m a_{km}  $$

where $\mu_k$ is the chemical potential of species $k$, $a_{km}$ is the number of atoms
of element $m$ in species $k$, and $\lambda_m$ is the chemical potential of the
elemental species per atom (the so-called "element potential"). In other words, the
chemical potential of each species in an equilibrium state is a linear sum of
contributions from each atom. We can see that this is true in the output above---the
chemical potential of H2 is exactly twice that of H, the chemical potential for OH is
the sum of the values for H and O, the value for H2O2 is twice as large as the value for
OH, and so on.

## Computing reaction rates and transport properties

The following program demonstrates the typical method for accessing the following object
types from a {ct}`Solution` object:

- {ct}`ThermoPhase`: Represents the thermodynamic properties of mixtures containing one
  or more species. Accessed using {ct}`Solution::thermo` method.
- {ct}`Kinetics`: Represents a kinetic mechanism involving one or more phases. Accessed
  using the {ct}`Solution::kinetics`.
- {ct}`Transport`: Computes transport properties for a phase. Accessed using the
  {ct}`Solution::transport` method.

```{literalinclude} kinetics_transport.cpp
:language: c++
```

This program produces the output below:

```none
Net reaction rates for reactions involving CO2
 11  CO + O (+M) <=> CO2 (+M)         3.54150687e-08
 13  HCO + O <=> CO2 + H              1.95679990e-11
 29  CH2CO + O <=> CH2 + CO2          3.45366954e-17
 30  CO + O2 <=> CO2 + O              2.70102741e-13
 98  CO + OH <=> CO2 + H              6.46935827e-03
119  CO + HO2 <=> CO2 + OH            1.86807592e-10
131  CH + CO2 <=> CO + HCO            9.41365868e-14
151  CH2(S) + CO2 <=> CH2 + CO2       3.11161343e-12
152  CH2(S) + CO2 <=> CH2O + CO       2.85339294e-11
225  NCO + O2 <=> CO2 + NO            3.74127381e-19
228  NCO + NO <=> CO2 + N2            6.25672710e-14
261  HNCO + O <=> CO2 + NH            6.84524918e-13
267  HNCO + OH <=> CO2 + NH2          7.78871222e-10
279  CO2 + NH <=> CO + HNO           -3.30333709e-09
281  NCO + NO2 <=> CO2 + N2O          2.14286657e-20
282  CO2 + N <=> CO + NO              6.42658345e-10
289  CH2 + O2 => CO2 + 2 H            1.51032305e-18
304  CH2CHO + O => CH2 + CO2 + H      1.00331721e-19

T        viscosity     thermal conductivity
------   -----------   --------------------
300.0    1.6701e-05    4.2143e-02
400.0    2.0896e-05    5.2797e-02
500.0    2.4704e-05    6.2827e-02
600.0    2.8230e-05    7.2625e-02
700.0    3.1536e-05    8.2311e-02
```

## Examples of additional C++ functionality

- [`combustor.cpp`](/examples/cxx/combustor): An example that shows how to set up a
  reactor network model, using classes {ct}`Reactor`, {ct}`ReactorNet`, {ct}`Reservoir`,
  {ct}`MassFlowController`, and {ct}`Valve`.
- [`flamespeed.cpp`](/examples/cxx/flamespeed): An example simulating a freely
  propagating laminar flame using the {ct}`Sim1D` solver.
