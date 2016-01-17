class Solution(ThermoPhase, Kinetics, Transport):
    """
    A class for chemically-reacting solutions.

    Instances can be created to represent any type of solution -- a
    mixture of gases, a liquid solution, or a solid solution, for
    example.

    Class Solution derives from classes `ThermoPhase`, `Kinetics`, and
    `Transport`.  It defines very few methods of its own, and is provided so
    that a single object can be used to compute thermodynamic, kinetic, and
    transport properties of a solution.
    """

class Interface(InterfacePhase, InterfaceKinetics):
    """
    Two-dimensional interfaces.

    Instances of class Interface represent reacting 2D interfaces
    between bulk 3D phases. Class Interface defines no methods of its
    own. All of its methods derive from either `InterfacePhase` or
    `InterfaceKinetics`.
    """

class DustyGas(ThermoPhase, Kinetics, DustyGasTransport):
    """
    A composite class which models a gas in a stationary, solid, porous medium.

    The only transport properties computed are the multicomponent diffusion
    coefficients. The model does not compute viscosity or thermal conductivity.

    """
