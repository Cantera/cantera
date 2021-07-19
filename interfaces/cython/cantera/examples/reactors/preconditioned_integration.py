# -*- coding: utf-8 -*-
"""
Ideal gas, constant-pressure, adiabatic kinetics simulation that compares preconditioned and non-preconditioned integration of nDodecane.

Requires: cantera >= 2.6.0
Keywords: combustion, reactor network, preconditioner
"""
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer

def integrate_reactor(n_reactors=15, preconditioner=True):
    """
        Compare the integrations of a preconditioned reactor network and a
        non-preconditioned reactor network. Increase the number of reactors in the
        network with the keyword argument `n_reactors` to see greater differences in
        performance.
    """
    # Use a reduced n-dodecane mechanism with PAH formation pathways
    gas = ct.Solution('nDodecane_Reitz.yaml', 'nDodecane_IG')
    # Create multiple reactors and set initial contents
    gas.TP = 1000, ct.one_atm
    gas.set_equivalence_ratio(0.75, 'c12h26', 'n2:3.76, o2:1.0')
    reactors = [ct.IdealGasConstPressureMoleReactor(gas) for n in range(n_reactors)]
    # set volume for reactors
    for r in reactors:
        r.volume = 0.1
    # Create reactor network
    sim = ct.ReactorNet(reactors)
    # Add preconditioner
    if preconditioner:
        sim.derivative_settings = {"skip-third-bodies":True, "skip-falloff":True}
        sim.preconditioner = ct.AdaptivePreconditioner()
    sim.initialize()
    # Advance to steady state
    integ_time = default_timer()
    # indexs to keep
    tidx = reactors[0].component_index("temperature")
    cidx = reactors[0].component_index("c12h26")
    coidx = reactors[0].component_index("co2")
    time = []
    T = []
    CO2 = []
    C12H26 = []
    # advance to steady state manually
    while (sim.time < 0.1):
        previous_state = sim.get_state()
        T.append(previous_state[tidx])
        CO2.append(previous_state[coidx])
        C12H26.append(previous_state[cidx])
        time.append(sim.time)
        sim.step()
    integ_time = default_timer() - integ_time
    # Return time to integrate
    if preconditioner:
        print("Preconditioned Integration Time: {:f}".format(integ_time))
    else:
        print("Non-preconditioned Integration Time: {:f}".format(integ_time))
    # Get and output solver stats
    lin_stats = sim.linear_solver_stats
    nonlin_stats = sim.nonlinear_solver_stats
    for key in lin_stats:
        print(key, lin_stats[key])
    for key in nonlin_stats:
        print(key, nonlin_stats[key])
    print("\n")
    return time, T, CO2, C12H26

if __name__ == "__main__":
    # integrate both to steady state
    timep, Tp, CO2p, C12H26p = integrate_reactor(preconditioner=True)
    timenp, Tnp, CO2np, C12H26np  = integrate_reactor(preconditioner=False)
    # plot selected state variables
    fig, axs = plt.subplots(1, 3)
    fig.tight_layout()
    ax1, ax2, ax3 = axs
    # temperature plot
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Temperature")
    ax1.plot(timenp, Tnp, linewidth=2)
    ax1.plot(timep, Tp, linewidth=2, linestyle=":")
    ax1.legend(["Normal", "Preconditioned"])
    # CO2 plot
    ax2.set_xlabel("Time")
    ax2.set_ylabel("CO2")
    ax2.plot(timenp, CO2np, linewidth=2)
    ax2.plot(timep, CO2p, linewidth=2, linestyle=":")
    ax2.legend(["Normal", "Preconditioned"])
    # C12H26 plot
    ax3.set_xlabel("Time")
    ax3.set_ylabel("C12H26")
    ax3.plot(timenp, C12H26np, linewidth=2)
    ax3.plot(timep, C12H26p, linewidth=2, linestyle=":")
    ax3.legend(["Normal", "Preconditioned"])
    plt.show()
