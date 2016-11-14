**********************************
Computing Thermodynamic Properties
**********************************

Class ThermoPhase
=================

Cantera can be used to compute thermodynamic properties of pure substances,
solutions, and mixtures of various types, including ones containing multiple
phases. The first step is to create an object that represents each phase. A
simple, complete program that creates an object representing a gas mixture and
prints its temperature is shown below:

.. code-block:: c++

    #include "cantera/thermo.h"
    #include <iostream>

    int main(int argc, char** argv)
    {
        std::unique_ptr<Cantera::ThermoPhase> gas(
            Cantera::newPhase("h2o2.cti", "ohmech"));
        std::cout << gas->temperature() << std::endl;
        return 0;
    }

Class :ct:`ThermoPhase` is the base class for Cantera classes that represent
phases of matter. It defines the public interface for all classes that represent
phases. For example, it specifies that they all have a method :ct:`temperature
<Phase::temperature>` that returns the current temperature, a method
:ct:`setTemperature(double T) <Phase::setTemperature>` that sets the
temperature, a method :ct:`getChemPotentials(double* mu)
<ThermoPhase::getChemPotentials>` that writes the species chemical potentials
into array ``mu``, and so on.

Class ThermoPhase can be used to represent the intensive state of any
single-phase solution of multiple species. The phase may be a bulk,
three-dimensional phase (a gas, a liquid, or a solid), or it may be a
two-dimensional surface phase, or even a one-dimensional "edge" phase. The
specific attributes of each type of phase are specified by deriving a class from
:ct:`ThermoPhase` and providing implementations for its virtual methods.

Cantera has a wide variety of models for bulk phase currently. Special attention
(in terms of the speed of execution) has been paid to an ideal gas phase
implementation, where the species thermodynamic polynomial representations
adhere to either the NASA polynomial form or to the Shomate polynomial
form. This is widely used in combustion applications, the original application
that Cantera was designed for. Recently, a lot of effort has been placed into
constructing non-ideal liquid phase thermodynamics models that are used in
electrochemistry and battery applications. These models include a Pitzer
implementation for brines solutions and a Margules excess Gibbs free energy
implementation for molten salts.

The Intensive Thermodynamic State
---------------------------------

Class :ct:`ThermoPhase` and classes derived from it work only with the intensive
thermodynamic state. That is, all extensive properties (enthalpy, entropy,
internal energy, volume, etc.) are computed for a unit quantity (on a mass or
mole basis). For example, there is a method :ct:`enthalpy_mole()` that returns
the molar enthalpy (J/kmol), and a method :ct:`enthalpy_mass()` that returns the
specific enthalpy (J/kg), but no method *enthalpy()* that would return the total
enthalpy (J). This is because class ThermoPhase does not store the total amount
(mass or mole) of the phase.

The intensive state of a single-component phase in equilibrium is fully
specified by the values of any :math:`r+1` independent thermodynamic properties,
where :math:`r` is the number of reversible work modes. If the only reversible
work mode is compression (a "simple compressible substance"), then two
properties suffice to specify the intensive state. Class ThermoPhase stores
internally the values of the *temperature*, the *mass density*, and the *mass
fractions* of all species. These values are sufficient to fix the intensive
thermodynamic state of the phase, and to compute any other intensive properties.
This choice is arbitrary, and for most purposes you can't tell which properties
are stored and which are computed.

Derived Classes
---------------

Many of the methods of ThermoPhase are declared virtual, and are meant to be
overloaded in classes derived from ThermoPhase. For example, class
:ct:`IdealGasPhase` derives from :ct:`ThermoPhase`, and represents ideal gas
mixtures.

Although class ThermoPhase defines the interface for all classes representing
phases, it only provides implementations for a few of the methods. This is
because ThermoPhase does not actually know the equation of state of any
phase---this information is provided by classes that derive from ThermoPhase.
The methods implemented by ThermoPhase are ones that apply to all phases,
independent of the equation of state. For example, it implements methods
``temperature()`` and ``setTemperature()``, since the temperature value is
stored internally.

* `Classes which inherit from ThermoPhase <../../../doxygen/html/group__thermoprops.html>`_
* `Classes which handle standard states for species <../../../doxygen/html/group__spthermo.html>`_


Example Program
===============

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
