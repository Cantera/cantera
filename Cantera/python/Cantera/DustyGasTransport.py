"""
Dusty Gas model for transport in porous media.

"""

from Cantera.Transport import Transport

class DustyGasTransport(Transport):
    def __init__(self, phase = None):
        Transport.__init__(self, model = "DustyGas", phase = phase)

    def setPorosity(self, porosity):
        """Set the porosity."""
        self.setParameters(0, 0, [porosity, 0.0])

    def setTortuosity(self, tortuosity):
        """Set the tortuosity."""
        self.setParameters(1, 0, [tortuosity, 0.0])

    def setMeanPoreRadius(self, pore_radius):
        """Set the mean pore radius."""
        self.setParameters(2, 0, [pore_radius, 0.0])

    def setMeanParticleDiameter(self, diameter):
        """Set the mean particle diameter."""
        self.setParameters(3, 0, [diameter, 0.0])                        

    def setPermeability(self, permeability):
        """Set the permeability. If not called, the value for close-packed
        spheres is used."""
        self.setParameters(4, 0, [permeability, 0.0])

##     def molarFluxes(self,
##                     conc = None,
##                     gradConc = None,
##                     gradPressure = 0.0):
##         self.setConcentrations(concentrations)
                    
    def set(self, **p):
        """Set model parameters. This is a convenience method that simply
        calls other methods depending on the keyword.

        Keywords:
        
        - porosity
        - tortuosity
        - pore_radius
        - diameter
        - permeability
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

    

    
