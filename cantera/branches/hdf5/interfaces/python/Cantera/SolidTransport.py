"""
Transport properties for solids.

This class implements a simple model for the diffusion coefficients and
the thermal conductivity of a solid. The diffusion coefficients have
modified Arrhenius form, and the thermal conductivity is constant.
All parameters are user-specified, not computed from a physical model.

Examples:

    >>> tr = SolidTransport(solid_phase)
    >>> tr.setThermalConductivity(0.5)  # W/m/K
    >>> tr.setDiffCoeff(species = "OxygenIon", A = 2.0, n = 0.0, E = 700.0)

Note that the diffusion coefficient is computed from D = A * T^n *
exp(-E/t) in m^2/s. Diffusion coefficients for unspecified species are
set to zero.

"""

from Cantera.Transport import Transport

class SolidTransport(Transport):
    def __init__(self, phase = None):
        Transport.__init__(self, model = "Solid", phase = phase)

    def setThermalConductivity(self, lamb):
        self.setParameters(1, 0, [lamb, 0.0])

    def setDiffCoeff(self, species = "", A = 0.0, n = 0.0, E = 0.0):
        k = self._phase.speciesIndex(species)
        self.setParameters(0, k, [A, n, E])
