
****************
C++ Header Files
****************

Cantera provides some header files designed for use in C++ application
programs. These are designed to include those portions of Cantera needed for
particular types of calculations. For example, the header file ``equilibrium.h``
includes header files needed to do equilibrium calculations (specifically, files
``ChemEquil.h`` and ``MultiPhaseEquil.h``).

These headers are designed for use in C++ application programs, and are not
included by the Cantera core. The headers and their functions are:

``equilibrium.h``
    Chemical equilibrium.

``IdealGasMix.h``
    Provides class :ct:`IdealGasMix`.

``Interface.h``
    Provides class :ct:`Interface`.

``integrators.h``
    ODE Integrators.

``kinetics.h``
    Chemical kinetics.

``numerics.h``
    Classes for matrices.

``onedim.h``
    One-dimensional reacting flows.

``reactionpaths.h``
    Reaction path diagrams.

``transport.h``
    Transport properties.

``zerodim.h``
    Zero-dimensional reactor networks.
