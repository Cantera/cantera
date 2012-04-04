********************************
Thermodynamic Properties Program
********************************

In the program below, a gas mixture object is created, and a few thermodynamic
properties are computed and printed out:

.. literalinclude:: thermodemo.cpp
   :language: c++

Note that the methods that compute the properties take no input parameters. The
properties are computed for the state that has been previously set and stored
internally within the object.
 
Naming Conventions
------------------

- methods that return *molar* properties have names that end in ``_mole``. 
- methods that return properties *per unit mass* have names that end in
  ``_mass``.
- methods that write an array of values into a supplied output array have names
  that begin with ``get``. For example, the method
  :ct:`ThermoPhase::getChemPotentials(double* mu)` writes the species chemical
  potentials into the output array ``mu``.

The thermodynamic property methods are declared in class :ct:`ThermoPhase`,
which is the base class from which all classes that represent any type of phase
of matter derive.

See :ct:`ThermoPhase` for the full list of available thermodynamic properties.
