"""
Dusty Gas model for transport in porous media.

"""

from Cantera.Transport import Transport

class DustyGasTransport(Transport):
    """The Dusty Gas transport model.  This class implements a
    transport manager for the Dusty Gas model for the effective
    transport properties of a gas in a stationary, solid, porous
    medium. The only properties computed are the multicomponent
    diffusion coefficients. The model does not compute viscosity or
    thermal conductivity.

    This class is a Python shadow class for Cantera C++ class
    DustyGasTransport.
    """


    def __init__(self, phase = None):
        """
        phase - The object representing the gas phase within the
        pores.
        """
        Transport.__init__(self, model = "DustyGas", phase = phase)

    def setPorosity(self, porosity):
        """Set the porosity. Internal. See: set"""
        self.setParameters(0, 0, [porosity, 0.0])

    def setTortuosity(self, tortuosity):
        """Set the tortuosity. Internal. See: set"""
        self.setParameters(1, 0, [tortuosity, 0.0])

    def setMeanPoreRadius(self, pore_radius):
        """Set the mean pore radius [m]. Internal. See: set"""
        self.setParameters(2, 0, [pore_radius, 0.0])

    def setMeanParticleDiameter(self, diameter):
        """Set the mean particle diameter [m]. Internal. See: set"""
        self.setParameters(3, 0, [diameter, 0.0])

    def setPermeability(self, permeability):
        """Set the permeability. If not called, the value for close-packed
        spheres is used. Internal."""
        self.setParameters(4, 0, [permeability, 0.0])


    def set(self, **p):
        """Set model parameters. This is a convenience method that simply
        calls other methods depending on the keyword.

        porosity - Porosity. Volume fraction of pores.

        tortuosity - Tortuosity. A measure of the extent to which the
        pores are straight cylinders (tortuosity = 1), or are more
        tortuous.

        pore_radius - The pore radius [m].

        All keywords are optional.
        """
        for o in p.keys():
            if o == "porosity":
                self.setPorosity(p[o])
            elif o == "tortuosity":
                self.setTortuosity(p[o])
            elif o == "pore_radius":
                self.setMeanPoreRadius(p[o])
            elif o == "diameter":
                self.setMeanParticleDiameter(p[o])
            elif o == "permeability":
                self.setPermeability(p[o])
            else:
                raise 'unknown parameter'
