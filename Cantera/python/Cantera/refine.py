"""Grid refinement.

Suppose you have a monotonic NumPy array of grid points 'z', and a
solution array soln[j,n] that contains 3 three solution components
denoted 'a', 'b', and 'c', evaluated at the grid points. To refine the
grid based on components 'a' and 'b' but not 'c', do the following.

>>> from refine import Refiner
>>> r = Refiner([(0, 'a'), (1, 'b')])
>>> new_grid, new_soln = r.refine(grid, soln)

"""


import Numeric
import math
from Cantera import CanteraError
from Cantera import interp

def eps():
    """Return the square root of machine precision."""
    e = 1.0
    while 1.0 + e <> 1.0: e = 0.5*e
    return math.sqrt(e)


def delta(f):
    """Given an array f, return an array of the difference in
    adjacent values."""
    n = len(f)
    d = Numeric.zeros(n-1,'d')
    for j in range(n-1):
        d[j] = f[j+1] - f[j]
    return d


def slope(z, f):
    """Given arrays z and f, return an array of the slopes df/dz in
    each interval."""
    n = len(z)
    s = Numeric.zeros(n-1,'d')
    for j in range(n-1):
        s[j] = ((f[j+1] - f[j])/(z[j+1] - z[j]))
    return Numeric.array(s,'d')


class RefineError(CanteraError):
    def __init__(self, msg):
        self.msg = 'Grid refinement error!\n'+msg


class Refiner:
    """Grid refiner.

    Attributes:

    components -- sequence of (number, name) pairs specifying the
    components of the solution to use for grid refinement. The number
    is used to access the component in the solution array, and the
    name is used only for diagnostic messages.

    max_delta -- Maximum tolerated difference in solution values
    between neighboring grid points, expressed as a fraction between 0
    and 1 of the total range of the component over all grid
    points. Default: 0.8 (minimal refinement).

    max_delta_slope -- Maximum tolerated difference in solution slopes
    between neighboring grid intervals, expressed as a fraction between 0
    and 1 of the total range of the component over all grid
    points. Default: 0.8 (minimal refinement).    
    
    """

    def __init__(self, components = [], delta = (2.0, 0.1, 0.2), names = []):
        self.components = components
        self.delta = delta
        self.names = names
        self.loglevel = 2
        self.eps = eps()
        self.min_range = 0.01
        self.direction = 1
        self.fctr = 1.0
        self.ok = 0


    def prune(self, grid = None, solution = None, threshold = None):

        n0 = len(grid)
        g = grid
        sol = solution
        self.fctr = 1.0
        savedir = self.direction
        
        ll = self.loglevel
        #self.loglevel = 0
        j = 1
        while j < len(g)-1:
            g0 = g
            s0 = sol
            pt = g[j]
            nn = len(g)
            
            # remove point j
            g = Numeric.take(g,range(0,j)+range(j+1,nn))

            # remove row j
            sol = Numeric.take(sol, range(0,j)+range(j+1,nn))
            np = len(g)

            self.direction = 1
            gnew, gn, snew, ok = self.refine(g, sol, threshold)
            if (len(gnew) > np):
                g = g0
                sol = s0
                j += 1
                if ll > 0:
                    print 'cannot remove point at ',pt
            else:
                if ll > 0:
                    print 'removed point at ',pt
        self.loglevel = ll
        self.fctr = 0.2
        self.direction = savedir
        return (g, sol)
    
            
    def refine(self, grid = None, solution = None, threshold = None, prune = 1):

        self.ok = 0
        # grid parameters
        n0 = len(grid)
        dz0 = grid[-1] - grid[0]

        maxpts = self.fctr*n0 + 1
        
        ncomp = Numeric.shape(solution)[1]
        
        if threshold:
            self.threshold = threshold
        else:
            self.threshold = self.eps * Numeric.ones(ncomp, 'd')
            
        if Numeric.shape(solution)[0] <> n0:
            raise RefineError('Number of solution points differs from '+
                              'number of grid points.')

        # if the solution components to examine for refinement have
        # not been specified, use all components.
        nc = Numeric.shape(solution)[1]
        if not self.components: self.components = range(nc)

        c = {}
        p = {}

        dz = delta(grid)
        for j in range(1,n0-1):
            if dz[j] > self.delta[0]*dz[j-1]:
                p[j] = 1
                c['point '+`j`] = 1
            if dz[j] < dz[j-1]/self.delta[0]:
                p[j-1] = 1
                c['point '+`j-1`] = 1                

        for i in self.components:
            try:
                name = self.names[i]
            except:
                name = 'component '+`i`
            
            # get component i at all points, and compute its slope
            v = solution[:,i]
            s = slope(grid, v)

            # compute the change in value and slope
            dv = delta(v)
            ds = delta(s)

            # find the range of values and slopes
            vmin = min(v)
            vmax = max(v)
            smin = min(s)
            smax = max(s)

            # max absolute values of v and s
            aa = max((abs(vmax), abs(vmin)))
            ss = max((abs(smax), abs(smin)))


            # refine based on component i only if the range of v is
            # greater than a fraction 'min_range' of max |v|. This
            # eliminates components that consist of small fluctuations
            # on a constant background.

            if (vmax - vmin) > self.min_range*aa:

                # maximum allowable difference in value between
                # adjacent points.
                
                dmax = self.delta[1]*(vmax - vmin) + self.threshold[i]
                for j in range(len(dv)):
                    r = abs(dv[j])/dmax
                    if r > 1.0:
                        p[j] = 1
                        c[name] = 1

                        
            # refine based on the slope of component i only if the
            # range of s is greater than a fraction 'min_range' of max
            # |s|. This eliminates components that consist of small
            # fluctuations on a constant slope background.
            
            if (smax - smin) > self.min_range*ss:

                # maximum allowable difference in slope between
                # adjacent points.
                dmax = self.delta[2]*(smax - smin)
                
                for j in range(len(ds)):
                    r = abs(ds[j]) / (dmax + self.threshold[i]/dz[j])
                    if r > 1:
                        c[name] = 1
                        p[j] = 1
                        p[j+1] = 1

        if len(p) == 0: self.ok = 1

        znew = []
        nnew = len(p)
        nadded = nnew
                
        if self.loglevel > 0:
            if nnew > 0:
                print '\nRefining grid.'
                print 'New points inserted after grid points',

        for j in range(n0 - 1):
            znew.append(grid[j])
            if p.has_key(j):
                if self.loglevel > 0: print j,
                znew.append(0.5*(grid[j] + grid[j+1]))
        if self.loglevel > 0: print
        znew.append(grid[-1])
        if self.loglevel > 0 and nnew > 0:
            print 'to resolve ',
            ck = c.keys()
            for s in ck:
                if s <> ck[-1]:
                    print s+',',
                else:
                    print s,
            print
        
        npts = len(znew)
        
        newsoln = Numeric.zeros((npts, ncomp),'d')
        for i in range(ncomp):
            for j in range(npts):
                newsoln[j,i] = interp.interp(znew[j],grid,solution[:,i])

        return (Numeric.array(znew), Numeric.array(znew), newsoln, self.ok)



def refine(grid = None, solution = None, components = [], delta = (0.8, 1.0), threshold = None):
    """Refine a grid and interpolate the solution onto the new grid."""
    r = Refiner(components = components, delta = delta)
    return r.refine(grid, solution, threshold)


def prune(grid = None, solution = None, components = [], delta = (0.8, 1.0), threshold = None):
    """Remove unneeded points from a grid and solution array."""
    r = Refiner(components = components, delta = delta)
    return r.prune(grid, solution, threshold)



# test it
if __name__ == '__main__':

    grid = Numeric.array([0.0, 0.2, 0.3, 1.0, 4.0])
    soln = Numeric.array([[100.0, 0.4, -9.0], 
                          [500.0, 0.0, -89.0], 
                          [700.0, 0.9, 99.0],
                          [-99.0, 8.0, 77.0],
                          [567.0, 8.0, 0.0]])
    grid_new, soln_new = refine(grid, soln, 
                                components = [0,2], 
                                delta = (0.5, 0.8))
    
    print 'new grid = ',grid_new
    print 'new solution = ',soln_new
    
                
            
    
