"""
Plotting of Flow1D solutions.
"""

from Numeric import shape, zeros, ones, array

class FlowPlotter:
    """FlowPlotter objects handle generating plots of Flow1D solutions.

    This class is primarily designed for internal use by class Flow1D.
    """
    
    def __init__(self, flow):
        self.flow = flow
        self.np = 0
        self.nv = 0
        self.names = ['z [m]', 'u [m/s]', 'V [1/s]',
                      'T [K]', 'lambda', 'P [Pa]'] + list(flow.gas.speciesNames())

    def value(self, j, n):
        if n == 0:
            return self.flow.z[j]
        elif n < 5:
            return self.flow.x[j,n-1]
        elif n == 5:
            return self.flow.gas.pressure()
        else:
            return self.mf[j,n-6]
        
        
    def plot(self, fmt = 'TECPLOT',
             fname = 'plot.dat', title = 'plot',
             zone = 'zone1',
             moles = 1,
             append = 0):

        self.np, self.nv = shape(self.flow.x)
        self.mf = zeros((self.np, self.nv - 4), 'd')
        if moles:
            for j in range(self.np):
                self.flow.setGas(j)
                self.mf[j,:] = self.flow.gas.moleFractions()
        else:
            for j in range(self.np):
                self.flow.setGas(j)
                self.mf[j,:] = self.flow.gas.massFractions()                
        
        if fmt == 'TECPLOT':
            from tecplot import write_TECPLOT_zone
            data = zeros((self.np,self.flow.gas.nSpecies()+6),'d')
            data[:,0] = self.flow.z
            data[:,1:5] = self.flow.x[:,0:4]
            data[:,5] = ones(self.np)*self.flow.gas.pressure()
            data[:,6:] = self.mf
            
            write_TECPLOT_zone(fname, title, zone, self.names,
                               self.np, self.nv+2, append, data)
            
        elif fmt == 'EXCEL':
            from excel import write_CSV_data
            write_CSV_data(fname, self.names,
                               self.np, self.nv+2, append, self)
        else:
            raise 'unknown format'+fmt
            
