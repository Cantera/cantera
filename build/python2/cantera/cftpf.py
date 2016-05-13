# coding: utf-8

import cantera as ct
import numpy as np

class CounterflowTwinPremixedFlame(ct.FlameBase):
    """ A twin premixed counterflow flame. Two opposed jets shooting into each other.
        Tweaked from CounterFlowPremixedFlame.py that shipped with Cantera 2.2.0
    """
    __slots__ = ('reactants', 'flame', 'products')

    def __init__(self, gas, grid=None):
        """                                                                                                             
        :param gas:                                                                                                     
            `Solution` (using the IdealGas thermodynamic model) used to                                                 
            evaluate all gas properties and reaction rates.                                                             
        :param grid:                                                                                                    
            Array of initial grid points                                                                                
                                                                                                                        
        A domain of class `AxisymmetricStagnationFlow` named ``flame`` will                                             
        be created to represent the flame. The three domains comprising the                                             
        stack are stored as ``self.reactants``, ``self.flame``, and                                                     
        ``self.products``.                                                                                              
        """
        self.reactants = ct.Inlet1D(name='reactants', phase=gas)
        self.reactants.T = gas.T

        self.flame = ct.AxisymmetricStagnationFlow(gas, name='flame')

        #The right most boundary will be a symmetry plane
        self.products = ct.SymmetryPlane1D(name='products', phase=gas)
        self.products.T = gas.T

        super(CounterflowTwinPremixedFlame, self).__init__(
                (self.reactants, self.flame, self.products), gas, grid)

        # Setting X needs to be deferred until linked to the flow domain                                                
        self.reactants.X = gas.X

    def set_initial_guess(self, equilibrate=True):
        """                                                                                                             
        Set the initial guess for the solution.                                                                         

        If `equilibrate` is True, then the products composition and temperature                                         
        will be set to the equilibrium state of the reactants mixture.                                                  
        """
        super(CounterflowTwinPremixedFlame, self).set_initial_guess()

        Yu = self.reactants.Y
        Tu = self.reactants.T
        self.gas.TPY = Tu, self.flame.P, Yu
        rhou = self.gas.density
        uu = self.reactants.mdot / rhou

        self.gas.equilibrate('HP')
        Teq = self.gas.T
        Yeq = self.gas.Y

        if equilibrate:
            Tb = Teq
            Yb = Yeq
            self.products.T = Tb
        else:
            Tb = self.products.T
            Yb = self.products.Y

        self.gas.TPY = Tb, self.flame.P, Yb
        rhob = self.gas.density
        ub = self.products.mdot / rhob

        locs = np.array([0.0, 0.4, 0.6, 1.0])
        self.set_profile('T', locs, [Tu, Tu, Teq, Tb])
        for k in range(self.gas.n_species):
            self.set_profile(self.gas.species_name(k), locs,
                             [Yu[k], Yu[k], Yeq[k], Yb[k]])

        # estimate strain rate                                                                                          
        self.gas.TPY = Teq, self.flame.P, Yeq
        zz = self.flame.grid
        dz = zz[-1] - zz[0]
        a = (uu + ub)/dz
        # estimate stagnation point                                                                                     
        x0 = rhou*uu * dz / (rhou*uu + rhob*ub)

        self.set_profile('u', [0.0, 1.0], [uu, -ub])
        self.set_profile('V', [0.0, x0/dz, 1.0], [0.0, a, 0.0])


def solveOpposedFlame(oppFlame, massFlux=0.12, tol_ss = [1.0e-7, 1.0e-13], tol_ts = [1.0e-7, 1.0e-11],\
                      loglevel = 1, \
                      ratio = 3, slope = 0.1, curve = 0.2, prune = 0.02):
    """ 
    Execute this function to run the Oppposed Flow Simulation 

    This function takes a CounterFlowTwinPremixedFlame object
    """

    oppFlame.reactants.mdot = massFlux
    oppFlame.products.mdot = massFlux

    oppFlame.flame.set_steady_tolerances(default=tol_ss)
    oppFlame.flame.set_transient_tolerances(default=tol_ts)
    oppFlame.set_initial_guess()  # assume adiabatic equilibrium products
    oppFlame.show_solution()

    oppFlame.energy_enabled = False
    oppFlame.solve(loglevel, False)

    oppFlame.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=prune)
    oppFlame.energy_enabled = True
    oppFlame.solve(loglevel)

    #Compute the strain rate, just before the flame. It also turns out to the maximum. This is the strain rate that computations comprare against, like when plotting Su vs. K 
    peakStrain = np.max(np.gradient(oppFlame.u, np.gradient(oppFlame.grid)))
    return (np.max(oppFlame.T), peakStrain)
    
