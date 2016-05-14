
****************
C++ Header Files
****************

Cantera provides some header files designed for use in C++ application
programs. These are designed to include those portions of Cantera needed for
particular types of calculations.

These headers are designed for use in C++ application programs, and are not
included by the Cantera core. The headers and their functions are:

``IdealGasMix.h``
    Provides class :ct:`IdealGasMix`.

``Interface.h``
    Provides class :ct:`Interface`.

``integrators.h``
    ODE Integrators.

``kinetics.h``
    Base kinetics classes and functions for creating :ct:`Kinetics` objects from
    input files.

``onedim.h``
    One-dimensional reacting flows.

``reactionpaths.h``
    Reaction path diagrams.

``thermo.h``
    Base thermodynamic classes and functions for creating :ct:`ThermoPhase`
    objects from input files.

``transport.h``
    Base transport property classes and functions for creating :ct:`Transport`
    objects from input files.

``zerodim.h``
    Zero-dimensional reactor networks.
