Creating ThermoPhase, Kinetics, and Transport objects
=====================================================

The following program demonstrates the general method for creating the following
object types:

- `ThermoPhase` - represents the thermodynamic properties of mixtures containing
  one or more species)
- `Kinetics` - represents a kinetic mechanism involving one or more phases)
- `Transport` - computes transport properties for a `ThermoPhase`

This program uses "factory" functions to create derived objects objects of the
appropriate type which are specified in the input file `gri30.cti`.

.. literalinclude:: demo1b.cpp
   :language: c++

This program produces the output below::

    Net reaction rates for reactions involving CO2
     11  CO + O (+M) <=> CO2 (+M)         3.54150724e-08
     13  HCO + O <=> CO2 + H              1.95680014e-11
     29  CH2CO + O <=> CH2 + CO2          3.45366988e-17
     30  CO + O2 <=> CO2 + O              2.70102522e-13
     41  CO2 + 2 H <=> CO2 + H2           3.45305359e-08
     98  CO + OH <=> CO2 + H              6.46935907e-03
    119  CO + HO2 <=> CO2 + OH            1.86807529e-10
    131  CH + CO2 <=> CO + HCO            9.41365695e-14
    151  CH2(S) + CO2 <=> CH2 + CO2       3.11161382e-12
    152  CH2(S) + CO2 <=> CH2O + CO       2.85339329e-11
    225  NCO + O2 <=> CO2 + NO            3.74127282e-19
    228  NCO + NO <=> CO2 + N2            6.25672779e-14
    261  HNCO + O <=> CO2 + NH            6.84524890e-13
    267  HNCO + OH <=> CO2 + NH2          7.78871264e-10
    279  CO2 + NH <=> CO + HNO           -3.30333658e-09
    281  NCO + NO2 <=> CO2 + N2O          2.14286686e-20
    282  CO2 + N <=> CO + NO              6.42658283e-10
    289  CH2 + O2 => CO2 + 2 H            1.51032319e-18
    304  CH2CHO + O => CH2 + CO2 + H      1.00331734e-19

    T        viscosity     thermal conductivity
    ------   -----------   --------------------
    300.0    1.6658e-05    4.2089e-02
    400.0    2.0861e-05    5.2537e-02
    500.0    2.4681e-05    6.2451e-02
    600.0    2.8218e-05    7.2157e-02
    700.0    3.1534e-05    8.1754e-02
