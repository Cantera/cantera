"""
Dusty Gas model for transport in porous media.

"""

from Cantera.Transport import Transport

class DustyGasTransport(Transport):
    def __init__(self, phase = None):
        Transport.__init__(self, model = "DustyGas", phase = phase)

    def setPorosity(self, porosity):
        self.setParameters(0, 0, [porosity, 0.0])

    def setTortuosity(self, tortuosity):
        self.setParameters(1, 0, [tortuosity, 0.0])

    def setMeanPoreRadius(self, pore_radius):
        self.setParameters(2, 0, [pore_radius, 0.0])

    def setMeanParticleDiameter(self, diameter):
        self.setParameters(3, 0, [diameter, 0.0])                        

    def setPermeability(self, permeability):
        self.setParameters(4, 0, [permeability, 0.0])
        
    def set(self, **p):
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

    
    
