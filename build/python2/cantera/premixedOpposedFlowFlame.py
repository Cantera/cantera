# coding: utf-8
import cantera as ct

if ct.__version__ == '2.3.0a2':
    print "STOP: This code does not work for the development version yet. Revert to 2.2.0"

#Import the CounterFlowPremixedFlame class and solver
import cftpf

#Use the 16 species Smooke-Giovangigli mechanism. Fast and good for CH4
gas = ct.Solution('smooke.cti')

#Create a CH4/Air premixed mixture with equivalence ratio=0.75, and at room temp and pressure.
gas.X = {'CH4':0.75, 'O2':2.0, 'N2':7.52}
gas.TP = 300, ct.one_atm

#Set initial conditions. There are other initial conditions as well that assume default values to avoid the user mulling over too many options
#To see a complete list execute (assuming you are working in an ipython environment) execute: cftpf.solveOpposedFlame?

#Set the velocity 
axial_velocity = 25 #in cm/s
                                                                                                
#Initial grid: 2.5cm, meaning the whole domain is 5cm wide
initial_grid = np.linspace(0.0, 0.025, 10)

#Done with initial conditions

#Compute the mass flux, as this is what the Flame object requires
massFlux = gas.density*axial_velocity/100 #units kg/m2/s  
#Create the flame object
oppFlame = cftpf.CounterflowTwinPremixedFlame(gas, grid=initial_grid)
#Set a guess for flame temperature
oppFlame.set_initial_guess(equilibrate=True)

#Now run the solver

#The solver returns the temperature and strain rate at the flame as default values
#However you can plot/see all state space variables by calling oppFlame.X where X is T, Y[i] or whatever
#The spatial variable (distance in metres) is in oppFlame.grid
#Thus to plot temperature vs distance, use oppFlame.grid and oppFlame.T

(T, K) = cftpf.solveOpposedFlame(oppFlame, massFlux)
print "Strain Rate: {0}".format(K)
