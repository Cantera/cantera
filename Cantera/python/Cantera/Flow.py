"""
One-dimensional reacting flows.
"""

from Cantera import getCanteraError
from FlowPlotter import FlowPlotter

from exceptions import *
import refine
import sys, types, copy, tempfile
import _cantera
import interp
import math

from Numeric import array, zeros, ones, transpose, size, sort, shape, asarray

_flows = {'Stagnation':0, 'Stag':0,
          'OneDimensional':1, '1D':1, 'OneD':1, 'OneDim':1,
          'Free':2}

_geom = {'Axisymmetric':0, 'Axi':0, 'Planar':1}



class Flow1D:
    """ One-dimensional reacting flows.

    Class Flow1D simulates a one-dimensional flow domain. To use
    Flow1D objects, they must be installed in a container, which is an
    object of class OneDim. Each Flow1D domain must be terminated by
    boundary domains.
    
    Class Flow1D can model several types of steady 'one dimensional'
    reacting flows. The flows are one-dimensional in the sense that
    the governing equations for the steady-state solution can be cast
    in the form of a set of ordinary differential equations in one
    axial coordinate (z). For the case of stagnation flows, this
    results from a similarity transformation that reduces the
    physically two-dimensional problem to one that is mathematically
    one-dimensional.

    The types of flows that may be simulated are:

    - One-dimensional reacting flows, such as burner-stabilized
    premixed flames;
    - Axisymmetric stagnation-point flows
    - Planar stagnation-point flows

    """


    # Allowed option keywords are defined here so that only these
    # keywords may be added to the _opt dictionary.
    
    _timeint_options = ['ftime', 'min_timestep', 'max_timestep',
                        'nsteps', 'timestep']
    _newton_options = ['max_jac_age', 'rtol', 'atol']
    _output_options = ['loglevel', 'plotfile']
    
    _options = _newton_options + _timeint_options + _output_options


    
    def __init__(self,
                 flow_type = 'Stagnation',
                 flow_geom = 'Axisymmetric',
                 gas = None,
                 grid = None,
                 pressure = 1.01325e5):
        """Flow1D Constructor.
        """

        self.__flow_id = -1
        self.domainType = 0
        self.loglevel = 1
        self.initial = {}
        if not grid:
            raise CanteraError('Grid not specified!')

        if not gas:
            raise CanteraError('Gas mixture object not specified!')
        self.gas = gas
        self.nsp = self.gas.nSpecies()
        self.nv = self.nsp + 4
        
        if _flows.has_key(flow_type):
            self.type = _flows[flow_type]
        else:
            raise CanteraError('unsupported flow type: '+flow_type)
        
        if self.type == 0:
            if flow_geom == 'Planar': self.type = 3
                

        # Create the kernel object. This is an instance of a subclass of
        # C++ Cantera class 'StFlow'.
        self.__flow_id = _cantera.Flow(self.type, self.gas.phase_id(), len(grid))
        

        # Set the grid. Must be done _after_ creating the kernel
        # object, since it sets the grid there too.
        self.setGrid(asarray(grid))


        # Set the pressure. Note that the pressure is constant
        # throughout the flowfield, due to the assumption of low Mach
        # number.
        self.setPressure(pressure)

                
        # set the thermo, kinetics, and transport managers to those of
        # object self.gas.
        self.setThermo(self.gas)
        self.setKinetics(self.gas)
        self.setTransport(self.gas)

        
        # create a NumPy array to hold the solution.  Since NumPy
        # arrays are stored by row, while Cantera expects arrays
        # stored by column, this array is defined with the grid point
        # number as the first index. In this way, all solution
        # variables for one grid point will be stored in contiguous
        # locations.        
        self.x = zeros((self.npts, self.nsp + 4), 'd')
        self.xnew = zeros((self.npts, self.nsp + 4), 'd')        
        

        # Default error tolerances
        self.setTolerances(
            u = (1.e-8, 1.e-15),
            V = (1.e-8, 1.e-15),
            T = (1.e-8, 1.e-15),
            L = (1.e-8, 1.e-15),
            Y = (1.e-8, 1.e-15))

        self.energy = 0
        self.names = ['u','V','T','L']+ list(self.gas.speciesNames())

        # finish setting default parameter values
        self.restoreDefaults()

        self.time = 0.0

        # For refiner, only refine the grid based on species mass
        # fraction and velocity profiles until the energy equation is
        # enabled.
        self.refine_components = range(4,4+self.nsp)
        
        self.refiner = refine.Refiner(components = self.refine_components,
                                      names = self.names)


        # Create an object to handle plotting results.
        self.plotter = FlowPlotter(self)

        #================>  end of method '__init__'  <===================


        
    def __del__(self):
        """Delete the Flow1D instance."""
        if self.__flow_id >= 0:
            _cantera.flow_delete(self.__flow_id)

    def shape(self):
        """Return (rows, columns) of solution matrix."""
        return (len(self.z), self.gas.nSpecies() + 4)

    def nPoints(self):
        """Number of grid points."""
        return self.npts
    
    def flow_id(self):
        """ID used to accesss the kernel object."""
        return self.__flow_id

    def option(self, key):
        """Return the value of an option."""
        return self._opt[key]
                         
    def setGrid(self, z):
        """Set the grid to the values in sequence z.

        The values will be sorted, and so the input sequence
        does not need to be monotonic.
        """ 
        self.z = array(z)
        sort(self.z)
        self.npts = len(z)
        return _cantera.flow_setupgrid(self.flow_id(), z)        


    def setThermo(self, th):
        """Set the thermodynamic property manager."""
        id = _cantera.flow_setthermo(self.flow_id(), th.phase_id())
        return id
    
    def setKinetics(self, kin):
        """Set the kinetics manager."""
        id = _cantera.flow_setkinetics(self.__flow_id, kin.ckin)
        return id

    def setTransport(self, tr, soret=0):
        """Set the transport manager."""
        id = _cantera.flow_settransport(self.__flow_id,
                                      tr.transport_id(), soret)
        return id

    def setPressure(self, p):
        """Set the pressure [Pa].

        Since the flow has very nearly the same pressure everywhere,
        this pressure value is used in all computations involving the
        equation of state.
        """
        self.p = p
        _cantera.flow_setpressure(self.__flow_id, p)


    def holdTemperature(self, points, t0):
        """Hold the temperature at grid points 'points' to t0,
        and disable the energy equation.
        """
        self.x[points,2] = t0
        _cantera.flow_settemperature(self.__flow_id, points, t0)
        self.setEnergyEqn('off')


    def holdMassFraction(self, j, k, y0):
        self.x[j,4+k] = y0
        _cantera.flow_setmassfraction(self.__flow_id, j, k, y0)
                

    def setInitialProfiles(self, datatable=None):

        """Set initial velocity, temperature, and/or species profiles.

        datatable -- Dictionary mapping variable names to sequences of
        (position, value) pairs. The position is specified in
        relative terms, as a number in the range [0,1], where the
        value zero corresponds to the smallest grid value, and the
        value one to the largest.

        The keys of datatable must be 'u', 'V', 'T', or a species
        name. Velocity and temperature values are entered in SI units,
        and species values are entered in arbitrary molar units, and
        will be normalized to produce mole fractions.

        The profile will be linearly interpolated onto the grid from
        the data provided. Each variable may be specified at different
        locations.

        Example:

        data = {}
        data['T'] = [(0, 500), (0.3, 2000), (0.8, 2200), (1, 1500)]
        data['u'] = [(0, 0.0), (1, 2)]
        data['H2'] = [(0, 0.2), (1, 0.3)]
        data['O2'] = [(0, 0.8), (1, 0.7)]
        flow.setInitialProfiles(data)

        """

        if datatable:
            self.datatable = datatable
        else:
            datatable = self.datatable
            
        vars = datatable.keys()
        x = zeros((self.npts, self.nsp),'d')
        equil = 0
        
        for var in vars:
            data = datatable[var]
            if not var == 'equil':
                data.sort()
                zz = []
                v = []
                for item in data:
                    zz.append(self.z[0] + item[0]*(self.z[-1] - self.z[0]))
                    v.append(item[1])
                    self.initial[var] = (zz, v)
                    
                if (var == 'u'):
                    for j in range(self.npts):
                        self.x[j,0] = interp.interp(self.z[j],zz,v)

                elif (var == 'V'):
                    for j in range(self.npts):
                        self.x[j,1] = interp.interp(self.z[j],zz,v)                    
                elif (var == 'T'):
                    for j in range(self.npts):
                        self.holdTemperature(j,interp.interp(self.z[j],zz,v))
                else:
                    k = self.gas.speciesIndex(var)
                    if k < 0:
                        raise CanteraError('Unknown species name: '+var)
                    for j in range(self.npts):
                        x[j,k] = interp.interp(self.z[j],zz,v)
            else:
                equil = 1
                
        if equil == 1:
            xin = self.left.X
            x[0,:] = xin
            for j in range(1,self.npts):
                self.gas.setState_TPX(self.x[j,2],self.p,xin)
                try:
                    self.gas.equilibrate('TP')
                except:
                    pass
                x[j,:] = self.gas.moleFractions()
                
        # convert input mole fractions to mass fractions,
        # and set the mass fraction profiles
        
        for j in range(self.npts):
            self.gas.setMoleFractions(x[j,:])
            y = self.gas.massFractions()
            for k in range(self.nsp):
                self.x[j, k+4] = y[k]
                self.holdMassFraction(j,k,y[k])
        self.enableSpecies()


    def regrid(self, grid):
        oldx = self.x
        oldgrid = self.z
        np, nv = shape(oldx)
        v = range(nv)
        self.setGrid(grid)
        for j in range(self.npts):
            for n in v:
                self.x[j,n] = interp.interp(self.z[j], oldgrid, oldx)           

    def __repr__(self):
        return self.show()
    
    def show(self, x = None):
        fname = tempfile.mktemp('.dat')
        x = self.x
        _cantera.flow_showsolution(self.__flow_id, fname, x)
        f = open(fname,'r')
        y = f.readlines()
        print
        for line in y:
            print line,
        print
        f.close()

    def showResid(self):
        """Print the current residual values"""
        print transpose(self.resid())
    
    def T(self,j):
        """Temperature at grid point j [K]."""
        return self.x[j,2] 

    def u(self,j):
        """Axial velocity at grid point j [m/s]."""
        return self.x[j, 0] 

    
    def V(self,j):
        """Radial velocity divided by radius at grid point j [1/s]."""
        return self.x[j, 1]
    
    def lamb(self,j):
        """(1/r)(dP/dr) at grid point j.

        If the solution has converged, this will be the same at all
        grid points.
        """
        return self.x[j, 3]

    def massFraction(self,sp,j):
        """Mass fraction of species 'sp', which may be referenced by
        name or by index number."""
        k = self.gas.speciesIndex(sp)
        return self.x[j, k + 4]        

    def density(self, j):
        """Density [kg/m^3]."""
        self.setGas(j)
        return self.gas.density()

    def molWt(self, j):
        """Mean molecular weight [kg/kmol]."""
        self.setGas(j)
        return self.gas.meanMolecularWeight()    
    
    def setGas(self, j):
        """Set the state of the internal gas mixture object to be
        consistent with the solution at grid point j."""
        self.gas.setTemperature(self.T(j))
        y = self.x[j,4:]  
        self.gas.setMassFractions(y)
        self.gas.setPressure(self.p)


    def setTolerances(self, u = None, V = None, T = None,
                      L = None, Y = None):
        """Set error tolerances.

        The inputs are tuples of (relative, absolute) tolerances for
        each of u, V, T, L, and Y.
        """
        default = (1.e-7, 1.e-15)
        if u == None: u = default
        if V == None: V = default
        if T == None: T = default        
        if L == None: L = default
        if Y == None: Y = default
        
        rtol = zeros(self.nsp+4,'d')
        atol = zeros(self.nsp+4,'d')
        
        rtol[0] = u[0]
        rtol[1] = V[0]
        rtol[2] = T[0]
        rtol[3] = L[0]
        for k in range(self.nsp):
            rtol[4+k] = Y[0]
            
        atol[0] = u[1]
        atol[1] = V[1]
        atol[2] = T[1]
        atol[3] = L[1]
        for k in range(self.nsp):
            atol[4+k] = Y[1]            

        _cantera.flow_settolerances(self.__flow_id, len(rtol), rtol,
                                  len(atol), atol)

    def disableSpecies(self):
        off = zeros(self.nsp,'d')
        _cantera.flow_solvespecies(self.__flow_id, self.nsp, off)

    def enableSpecies(self):
        o = ones(self.nsp,'d')
        _cantera.flow_solvespecies(self.__flow_id, self.nsp, o)        

        
    def setEnergyEqn(self, o, loglevel = 0, pt = -1):
        """Enable or disable the energy equation."""
        if o == 'on':
            self.energy = 1
            _cantera.flow_energy(self.__flow_id, pt, 1)
            if loglevel > 0: print '\n%%%%%%%%%%%%  Enabling energy equation  %%%%%%%%%%%%\n'
        elif o == 'off':
            self.energy = 0
            _cantera.flow_energy(self.__flow_id, pt, 0)
            if loglevel > 0: print '\n%%%%%%%%%%%%  Disabling energy equation  %%%%%%%%%%%%\n'            

    def setEnergyFactor(self, e):
        _cantera.flow_setenergyfactor(self.__flow_id, e)
            
##     def resid(self, point = -1):
##         """Return the residual vector.

##         If 'point' is specified, the residual is only evaluated at
##         this point and those adjacent to it. Otherwise, it is
##         evaluated at all grid points.
##         """
        
##         r = zeros(shape(self.x),'d')
##         if point >= 0 and point < self.npts:
##             j = point
##         else:
##             j = -1
##         _cantera.flow_eval(self.__flow_id, j, self.x, r)
##         return r


    def ssnorm(self):
        """Absolute maximum value of the steady-state residual for any
        component at any point."""
        a = self.container.ssnorm(self.x,self.xnew)
        return a

    def integrateChem(self, dt, loglevel=1):
        """Condition the species profiles by integrating the
        constant-pressure kinetics rate equations
        \[
           \dot Y_k = \dot\omega_k M_k / \rho.
        \]
        at each grid point for time 'dt'. 
        """
        if loglevel > 0:
            print '\nIntegrating the chemical ource terms for %10.4g s...' % dt,
        _cantera.flow_integratechem(self.__flow_id, self.x, dt)
        print _cantera.readlog()


    def restoreDefaults(self):
        """Restore default options."""
        self._opt = {}
        self.setOptions(

            max_jac_age = 5,
            timestep = 1.e-6,
            min_timestep = 1.e-12,
            max_timestep = 0.1,
            nsteps = 20,
            ftime = 3.0,

            plotfile = ""
            )
        
    def setOptions(self, **options):
        """Set options."""
        for kw in options.keys():
            if kw in Flow1D._options:
                self._opt[kw] = options[kw]
            else:
                raise OptionError(kw)
        
        
    def refine(self, loglevel = 2):
        """Refine the grid.

        """
        r = self.refiner
        r.components = range(4,self.nsp+4)
        if self.energy:
            r.components.append(2)        
            if self.type == 0:
                r.components += [0,1]

        #dsave = r.delta
        #while 1 > 0:
        #    try:
        znew, zadded, xn, ok = r.refine(grid = self.z, solution = self.x)
        nin = len(znew) - len(self.z)
        #        break
        #    except:
        #        r.delta = (0.9*r.delta[0], 0.9*r.delta[1])
        #r.delta = dsave
        
        if not ok:
            self.setGrid(znew)
            self.x = array(xn,'d')
            self.xnew = zeros(shape(xn),'d') 

            # update the fixed temperature values if the energy
            # equation is not being solved

            if self.energy == 0:
                for j in range(self.npts):
                    zz, tt = self.initial['T']
                    t = interp.interp(self.z[j],zz,tt)
                    self.holdTemperature(j,t)
            else:
                for j in range(self.npts):
                    self.holdTemperature(j,self.x[j,2])
                self.setEnergyEqn('on')  

        if loglevel > 0:
            print 'Refine: ',
            print 'added',nin,'points.'
            print 'Grid size = ',len(self.z)

        return nin

        
    def prune(self, loglevel = 2):
        """Prune the grid.

        """
        r = self.refiner
        r.components = range(4,self.nsp+4)
        if self.energy:
            r.components.append(2)        

        znew, xn = r.prune(grid = self.z, solution = self.x)
        nout = len(self.z) - len(znew)
        
        if nout > 0:
            self.setGrid(znew)
            self.x = array(xn,'d')
            self.xnew = zeros(shape(xn),'d') 

            # update the fixed temperature values if the energy
            # equation is not being solved

            if self.energy == 0:
                for j in range(self.npts):
                    zz, tt = self.initial['T']
                    t = interp.interp(self.z[j],zz,tt)
                    self.holdTemperature(j,t)
            else:
                for j in range(self.npts):
                    self.holdTemperature(j,self.x[j,2])
                self.setEnergyEqn('on')  

        if loglevel > 0:
            print 'Prune: ',
            print 'removed',nout,'points.'
            print 'Grid size = ',len(self.z)

        return nout

        
    def save(self, filename, id, desc="", append=0):
        """Save a solution to a file.

        filename -- file name. If the file, does not exist, it will
                    be created. The save files are xml files, and the
                    filename should have the extension '.xml'. If it
                    does not, this extension will be appended to the
                    name.
                    
        id       -- the ID tag of the solution. Multiple solutions may
                    be saved to the same file. Specifying a unique ID
                    tag allows this solution to selected later by
                    method 'restore'.
                    
        append   -- If append > 0, the solution will be appended to the
                    file. Otherwise, the file will be overwritten if it
                    exists.
        """
        appnd = append
        fn = filename
        extn = filename[-4:]
        if extn <> '.xml' and extn <> '.XML':
            fn = filename + '.xml'
            
        #if self.loglevel > 0:
        #print 'Solution saved to file',filename,'as solution',`id`
            
        #_cantera.flow_save(self.__flow_id, fn, id, appnd, self.x)
        _cantera.flow_save(self.__flow_id, fn, id, desc, self.x)
        print _cantera.readlog()

        np, nv = shape(self.x)
        f = open('ctsoln.dat','w')
        for j in range(np):
            f.write('%14.6e %14.6e ' % (100.0*self.z[j], self.x[j,2]))
            for k in range(4,nv):
                f.write('%14.6e ' % (self.x[j,k],))
            f.write('\n')
        f.close()


    def restore(self, filename, id):
        """Restore a previously-saved solution.

        filename -- name of a file containing a solution previously
                    saved by a call to 'save'
        id       -- the ID tag of the solution
        """

        #try:
        (z, s) = _cantera.flow_restore(self.__flow_id, 0, filename, id)
        #except:
        #   raise CanteraError()
        
        self.setGrid(z)
        self.x = array(s,'d')
        self.xnew = array(self.x,'d')
        
        for j in range(self.npts):
            self.holdTemperature(j,self.x[j,2])
        self.initial['T'] = (self.z, self.x[:,2])
                
        self.setEnergyEqn('off')
        if self.loglevel > 0:
            print 'Solution ',`id`,'read from file',filename
            print _cantera.readlog()


    def outputTEC(self, plotfile="", title="", zone="c0", append=0):
        """Write the current solution to a file in TECPLOT format.

        plotfile   --  file name (required)
        title      --  plot title
        zone       --  zone name
        append     --  if append > 0, the output is appended to the file
        
        """
        self.plotter.plot(fname = plotfile, title = title,
                          zone = zone, append=append)


    def outputCSV(self, plotfile="", append=0):
        """Write the current solution to a file in CSV format.

        plotfile   --  file name (required)
        append     --  if append > 0, the output is appended to the file
        
        """
        self.plotter.plot(fname = plotfile, fmt = 'EXCEL', 
                          append=append)        



    def plot(self, i):
        """Plot solution component i. Requires the scipy package."""
        from scipy import gplt
        return gplt.plot(self.z, self.x[:,i])


    
    def setBoundaries(self, left = None, right = None):
        """Install the boundary objects.

        The type of boundary object determines the boundary conditions.
        """
        nleft = 0
        nright = 0
        if left:
            self.left = left
            nleft = left.bdry_id()
        if right:
            self.right = right
            nright = right.bdry_id()
        _cantera.flow_setboundaries(self.__flow_id, nleft, nright)

"""
 $Author$
 $Revision$
 $Date$

 Copyright 2001  California Institute of Technology
"""

