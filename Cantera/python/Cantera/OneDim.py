
from exceptions import *
import refine
import sys, types, copy, tempfile
import interp
import math

import _cantera

from Numeric import *

def print_heading(msg):
    print '\n\n\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n'
    print msg
    print '\n'

    

class OneDim:
    """
    One-dimensional, multi-domain problems.

    Class OneDim allows solving multi-domain, one-dimensional,
    steady-state problems implicitly. Each domain has a set of one or
    more grid points, and each may have a different number of solution
    components.

    At ...  N(I) algebraic residual
    equations are defined for the N(I) solution components. Each
    domain may have a different number of components.

    The domains are linked in a linear chain, and are of two
    types. Standard domains know nothing of their neighbors, and
    evaluate their residual functions using only information in their
    own domain. 'Connector' domains can also modify the residual
    equations of their immediate neighbors, but only at the nearest
    grid point. Every standard domain must be attached to a connector
    at both ends. Connectors also serve to terminate the ends of the
    chain.
    
    """

    _timeint_options = ['ftime', 'min_timestep', 'max_timestep',
                        'nsteps', 'timestep', 'ts_jac_age']
    _newton_options = ['max_jac_age', 'rtol', 'atol']
    _output_options = ['loglevel', 'plotfile']
    
    _options = _newton_options + _timeint_options + _output_options

    
    def __init__(self, domains):
        """Create a new one-didmensional model from a list of domains. """

        # instance variables
        self._size = []
        self._start = []
        self._end = []
        self._domain = []     # all domains
        self._flow = []       # extended domains
        self._shape = []
        self._loc = 0
        self._opt = {}
        self.time = 0.0
        self.x = array([0.0,],'d')
        self._surf = []
        self.npts = []

        # local variables
        dtype = []          # list of integer domain types
        dlist = []          # list of integer domain ids

        # add each domain
        for d in domains:
            
            if d.domainType == 0:
                self.addFlow(d)
                dlist.append(d.flow_id())
                self.npts.append(d.nPoints())

            elif d.domainType == 1:
                self.addSurface(d)
                dlist.append(d.surf_id())
                self.npts.append(1)

            elif d.domainType == 2:
                self.addBoundary(d)
                dlist.append(d.bndry_id())
                self.npts.append(1)

            else:
                raise 'unknown domain type'
            dtype.append(d.domainType)

        self.__onedim_id = _cantera.onedim_new(len(dlist),
                                               array(dlist,'i'),
                                               array(dtype,'i'))
        self.collect()        
        self.restoreDefaults();
        self.ienergy = 0
        self.ts_jac_age = 50

        
    def __del__(self):
        """Delete the kernel object.

        This does not delete the individual domains."""
        _cantera.onedim_del(self.__onedim_id)


    def addFlow(self, flow):

        # add the domain to the list of all domains and to the list of
        # extended domains
        self._domain.append(flow)
        self._flow.append(flow)

        # set the index of this domain
        flow.index = len(self._domain) - 1

        np, nv = flow.shape()
        self._shape.append((np,nv))
        self._size.append(np*nv)
        self._start.append(self._loc)
        self._loc += np*nv
        self._end.append(self._loc)        


    def index(self, n, j, i):
        np, nv = self._shape[i]
        return self._start[n] + nv*j + i


##     def resetEnergy(self):
##         i = 0
##         ilast = self.ienergy
##         for f in self._flow:
##             i += f.resetEnergy()
##         self.ienergy = 0
##         return 0# ilast

    
    def finish(self):
        """
        Update the solution in each domain based on the global solution.

        This method is called by function 'solve' when a converged
        solution has been found, just prior to grid refinement.
        """
        for i in range(len(self._domain)):
            self._domain[i].x = self.solution(i)

            
    def solution(self, i):
        """ Return the solution array for domain i.

        The returned array has the shape (points,
        components) appropriate for domain i.  """
        
        x = self.x[self._start[i]:self._end[i]]
        dx = reshape(x, self._domain[i].shape())
        return dx


    def resid(self, i):
        """
        The residual matrix for domain i.

        The returned array has the shape
        (points, components) appropriate for domain i.
        """
        self.ssnorm()
        x = self.xnew[self._start[i]:self._end[i]]
        dx = reshape(x, self._domain[i].shape())
        return dx

        
    def addSurface(self, surf):
        """Add a surface domain."""
        self._surf.append(surf)
        self._domain.append(surf)
        surf.index = len(self._domain) - 1
        nv = surf.kin.nSpecies()
        np = 1
        self._shape.append((np,nv))
        self._size.append(np*nv)
        self._start.append(self._loc)
        self._loc += np*nv
        self._end.append(self._loc)


    def addBoundary(self, b):
        """Add a boundary domain."""
        #self._surf.append(surf)
        self._domain.append(b)
        b.index = len(self._domain) - 1
        nv = 2 # surf.kin.nSpecies()
        np = 1
        self._shape.append((np,nv))
        self._size.append(np*nv)
        self._start.append(self._loc)
        self._loc += np*nv
        self._end.append(self._loc)        

        
    def setNewtonOptions(self, max_jac_age = 5):
        _cantera.onedim_setnewtonoptions(self.__onedim_id, max_jac_age)


    def newton_solve(self, loglevel = 0):
        """Damped Newton iteration.

        This method invokes C++ method 'solve' of kernel class
        'OneDim' on the current solution. The solution is only
        modified if the damped Newton process leads to a fully
        converged solution.  Otherwise, an exception is raised.
        """

        iok = _cantera.onedim_solve(self.__onedim_id, self.x,
                                   self.xnew, loglevel)
        #if loglevel > 0: print _cantera.readlog()
        if iok >= 0:
            _cantera.copy(size(self.x),self.xnew,self.x)
        elif iok > -10:
            raise CanteraError()
        else:
            raise 'iok = '+`iok`
        return iok


        
    def collect(self):
        """Collect the state information from each domain to
        construct the global solution vector."""
        n = 0
        strt = []        # list of start locations for each domain
        self.npts = []
        nd = len(self._domain)
        for d in self._domain:
            strt.append(n)
            n += size(d.x)
            self.npts.append(d.shape()[0])
        strt.append(n)

        self.x = zeros(n,'d')
        self.xnew = zeros(n,'d')
        
        # set the portion of the global solution vector corresponding
        # to each domain to the flattened solution matrix for that
        # domain
        for i in range(nd):
            self.x[strt[i]:strt[i+1]] = reshape(self._domain[i].x,(-1,))

        
    def ssnorm(self):
        """Max norm of the steady-state residual."""
        n = _cantera.onedim_ssnorm(self.__onedim_id, self.x, self.xnew)
        return n

    
    def setSteadyMode(self):
        """Prepare to solve the steady-state problem."""
        return _cantera.onedim_setsteadymode(self.__onedim_id)

    
    def setTransientMode(self, dt):
        """Prepare for time-stepping with timestep dt.

        Must be called before each step."""
        return _cantera.onedim_settransientmode(self.__onedim_id, dt, self.x)


    def option(self, key):
        """Return the value of an option."""
        return self._opt[key]


    def restoreDefaults(self):
        """Restore default options."""
        self._opt = {}
        self.setOptions(
            max_jac_age = 20,
            ts_jac_age = 30,
            timestep = 1.e-6,
            min_timestep = 1.e-12,
            max_timestep = 0.1,
            nsteps = [1,2,4,8,20],
            ftime = 3.0,
            plotfile = ""
            )


    def setOptions(self, **options):
        """
        Set options.

        Time stepping:
           nsteps       --  number of steps.
           min_timestep --  minimum timestep
           max_timestep --  maximum timestep
           ftime        --  factor by which to increase the timestep for next
                            set of 'nsteps' timesteps

        Newton solver:
           max_jac_age  --  maximum number of times Jacobian will be used
                            before re-evaluating
           rtol         --  relative error tolerance
           atol         --  absolute error tolerance

        Output:
           loglevel     --  controls amount of diagnostic output
           plotfile     --  file to write plot data for intermediate solutions
           
        """
        for kw in options.keys():
            if kw in OneDim._options:
                self._opt[kw] = options[kw]
            else:
                raise OptionError(kw)
            if kw in OneDim._newton_options:
                self.setNewtonOptions(max_jac_age =
                                      self._opt['max_jac_age'])


    def refine(self, loglevel = 2):
        """Refine the grid of every flow domain."""
        new_points = 0
        for f in self._flow:
            new_points += f.refine(loglevel)
        if new_points > 0:
            self.collect()
            _cantera.onedim_resize(self.__onedim_id)
            self._shape = []
            self._size = []
            self._start = []
            self._end = []
            self._loc = 0
            for d in self._domain:
                np, nv = d.shape()
                self._shape.append((np, nv))
                self._size.append(np*nv)
                self._start.append(self._loc)
                self._loc += np*nv
                self._end.append(self._loc)
        return new_points

    def prune(self, loglevel = 2):
        """Prune the grid of every flow domain."""
        rem_points = 0
        for f in self._flow:
            rem_points += f.prune(loglevel)
        if rem_points > 0:
            self.collect()
            _cantera.onedim_resize(self.__onedim_id)
            self._shape = []
            self._size = []
            self._start = []
            self._end = []
            self._loc = 0
            for d in self._domain:
                np, nv = d.shape()
                self._shape.append((np, nv))
                self._size.append(np*nv)
                self._start.append(self._loc)
                self._loc += np*nv
                self._end.append(self._loc)
        return rem_points

    def setEnergyFactor(self, e):
        for f in self._flow:
            f.setEnergyFactor(e)
            
    def restore(self, n, file, soln, loglevel = 2):
        """Read the solution for domain n from a file."""
        self._domain[n].restore(file, soln)
        self.collect()
        _cantera.onedim_resize(self.__onedim_id)
        self._shape = []
        self._size = []
        self._start = []
        self._end = []
        self._loc = 0
        for d in self._domain:
            np, nv = d.shape()
            self._shape.append((np, nv))
            self._size.append(np*nv)
            self._start.append(self._loc)
            self._loc += np*nv
            self._end.append(self._loc)    


    def c_timeStep(self, nsteps, dt, loglevel = 0):
        dtnew = _cantera.onedim_timestep(self.__onedim_id, nsteps, dt,
                               self.x, self.xnew, loglevel)
        #print _cantera.readlog()
        return dtnew


    def py_timeStep(self, nsteps, dt, loglevel = 0):
        """Take time steps using Backward Euler.

        nsteps   -- number of steps
        dt       -- initial step size
        loglevel -- controls amount of printed diagnostics
        """

        self.setNewtonOptions(max_jac_age = self._opt['ts_jac_age'])
        print 'max jac age = ',self._opt['ts_jac_age']
        
        if loglevel > 0:
            print_heading('Begin time integration.\n\n')
            print(' step    size (s)    log10(ss) ')
            print('===============================')

        n = 0
        maxdt = self._opt['max_timestep']
        while n < nsteps:
            if loglevel > 0:
                ss = self.ssnorm()
                str = ' %4d  %10.4g  %10.4g' % (n,dt,math.log10(ss))
                print str,
            try:
                self.setTransientMode(dt)
                m = self.newton_solve(loglevel-1)
                self.time += dt
                n += 1
                if m == 100:
                    dt *= 1.5
                if dt > maxdt: dt = maxdt
                if loglevel > 0: print
                
            except CanteraError:
                #print self.resid(1)[:,0]
                if loglevel > 0: print '...failure.'
                dt *= 0.5
                if dt < 1.e-16:
                    self._domain[0].show()
                    raise CanteraError('Time integration failed.')

        self.setSteadyMode()
        self.setNewtonOptions(max_jac_age =
                              self._opt['max_jac_age'])        
        return dt

    def show(self):
        for d in self._domain:
            d.show()


    def showStatistics(self):
        _cantera.onedim_writestats(self.__onedim_id)
        #print _cantera.readlog()


    def save(self, filename, id, desc=""):
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
                    
        """
        fn = filename
        extn = filename[-4:]
        if extn <> '.xml' and extn <> '.XML':
            fn = filename + '.xml'
            
        _cantera.onedim_save(self.__onedim_id, fn, id, desc, self.x)
        #print _cantera.readlog()


