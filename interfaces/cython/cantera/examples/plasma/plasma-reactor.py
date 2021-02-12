"""
This example uses custom.py as the template to simulate a constant pressure
constant temperature plasma reactor.
"""

import cantera as ct
import numpy as np
import scipy.integrate

class ReactorOde:
    def __init__(self, gas):
        # Parameters of the ODE system and auxiliary data are stored in the
        # ReactorOde object.
        self.gas = gas
        self.P = gas.P

    def __call__(self, t, y):
        """the ODE function, y' = f(t,y) """

        # State vector is [T, Y_1, Y_2, ... Y_K]
        self.gas.set_unnormalized_mass_fractions(y[1:])
        self.gas.TP = y[0], self.P
        rho = self.gas.density
        h = self.gas.partial_molar_enthalpies / self.gas.molecular_weights
        Q_electron = ct.electron_charge * self.gas.electron_total_power_loss * self.gas.number_density("E")
        wdot = self.gas.net_production_rates

        dTdt = (-np.dot(self.gas.partial_molar_enthalpies, wdot) / rho + Q_electron / rho +
               1.0 / t_res * (np.dot(Y0, h0) - np.dot(Y0, h))) / self.gas.cp
        dYdt = wdot * self.gas.molecular_weights / rho + 1.0 / t_res * (Y0 - self.gas.Y)

        return np.hstack((dTdt, dYdt))

gas = ct.Plasma('oxygen_plasma.yaml')

# Initial condition
P = ct.one_atm
gas.TPX = 1000, P, 'O2:1.0, E:1e-9, O2^+:1e-9'
gas.electric_field = 2e5
gas.electric_field_freq = 1e9
gas.set_electron_energy_grid(np.linspace(0, 50, 500))
gas.set_initial_mean_electron_energy(2.0)
gas.set_reuse_EEDF(True)
y0 = np.hstack((gas.T, gas.Y))

# inlet composition
Y0 = gas.Y

# init partial mass enthalpies
h0 = gas.partial_molar_enthalpies / gas.molecular_weights

# residual time [s]
t_res = 1e-3

# Set up objects representing the ODE and the solver
ode = ReactorOde(gas)
solver = scipy.integrate.ode(ode)
solver.set_integrator('vode', method='bdf', with_jacobian=True)
solver.set_initial_value(y0, 0.0)

# Integrate the equations, keeping T(t) and Y(k,t)
t_end = 1e-2
states = ct.SolutionArray(gas, 1, extra={'t': [0.0]})
dt = 1e-5
while solver.successful() and solver.t < t_end:
    solver.integrate(solver.t + dt)
    gas.TPY = solver.y[0], P, solver.y[1:]
    states.append(gas.state, t=solver.t)

# Plot the results
try:
    import matplotlib.pyplot as plt
    L1 = plt.plot(states.t, states.T, color='r', label='T', lw=2)
    plt.xlabel('time (s)')
    plt.ylabel('Temperature (K)')
    plt.twinx()
    L2 = plt.plot(states.t, states('E').X, label='E', lw=2)
    plt.ylabel('Mole Fraction')
    plt.legend(L1+L2, [line.get_label() for line in L1+L2], loc='lower right')
    plt.show()
except ImportError:
    print('Matplotlib not found. Unable to plot results.')
