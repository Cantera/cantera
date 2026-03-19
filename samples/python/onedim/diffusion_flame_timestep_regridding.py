# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

"""
Timestep Regridding for a High-Pressure Diffusion Flame
=======================================================

This example compares two attempts to solve the same high-pressure hydrogen /
oxygen counterflow diffusion flame on a deliberately coarse initial grid. With
regrid-on-timestep-failure disabled, the solver gives up after the first failed
timestepping sequence. Allowing a few regrid retries lets the solver add points
and recover a converged solution.

Requires: cantera >= 4.0; matplotlib is optional for plotting

.. tags:: Python, combustion, 1D flow, diffusion flame, solver, grid refinement
"""

import time

import numpy as np

import cantera as ct

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


WIDTH = 30e-3
INITIAL_POINTS = 20
PRESSURE = 7e6
FUEL_MDOT = 0.3
OXIDIZER_MDOT_FACTOR = 3.0


def make_flame(regrid_max: int) -> ct.CounterflowDiffusionFlame:
    gas = ct.Solution("h2o2.yaml")
    flame = ct.CounterflowDiffusionFlame(
        gas, grid=np.linspace(0.0, WIDTH, INITIAL_POINTS))
    flame.max_time_step_count = 200
    flame.set_refine_criteria(ratio=2.0, slope=0.06, curve=0.08, prune=0.02)
    flame.time_step_regrid = regrid_max

    flame.P = PRESSURE
    flame.fuel_inlet.X = "H2:1.0"
    flame.fuel_inlet.T = 800.0
    flame.oxidizer_inlet.X = "O2:1.0"
    flame.oxidizer_inlet.T = 711.0

    rho_f = flame.fuel_inlet.phase.density
    rho_o = flame.oxidizer_inlet.phase.density
    flame.fuel_inlet.mdot = FUEL_MDOT
    flame.oxidizer_inlet.mdot = (FUEL_MDOT / rho_f) * rho_o * OXIDIZER_MDOT_FACTOR
    return flame


def solve_case(regrid_max: int) -> tuple[ct.CounterflowDiffusionFlame, bool, float, str]:
    flame = make_flame(regrid_max)
    flame.clear_stats()
    t0 = time.perf_counter()
    try:
        flame.solve(loglevel=0, auto=False)
        success = True
        message = ""
    except ct.CanteraError as err:
        success = False
        message = str(err).splitlines()[0]
    elapsed = time.perf_counter() - t0
    return flame, success, elapsed, message


def summarize(flame: ct.CounterflowDiffusionFlame, label: str, success: bool,
              elapsed: float, message: str) -> None:
    print(f"\n{label}")
    print(f"  success: {success}")
    print(f"  wall time [s]: {elapsed:.2f}")
    print(f"  final grid points: {len(flame.grid)}")
    print(f"  timesteps: {sum(flame.time_step_stats)}")
    print(f"  Jacobian evaluations: {sum(flame.jacobian_count_stats)}")
    if success:
        print(f"  peak temperature [K]: {np.max(flame.T):.1f}")
    else:
        print(f"  error: {message}")


baseline, baseline_success, baseline_time, baseline_message = solve_case(regrid_max=0)
regridded, regridded_success, regridded_time, regridded_message = solve_case(regrid_max=3)

summarize(baseline, "Regrid fallback disabled", baseline_success,
          baseline_time, baseline_message)
summarize(regridded, "Regrid fallback enabled", regridded_success,
          regridded_time, regridded_message)

if regridded_success and plt is not None:
    fig, ax = plt.subplots()
    ax.plot(regridded.grid * 1e3, regridded.T, "o-", label="Recovered solution")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("Temperature [K]")
    ax.set_title(
        "Counterflow diffusion flame at "
        f"{PRESSURE / 1e6:.1f} MPa with timestep regridding"
    )
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    plt.show()
elif regridded_success:
    print("\nInstall matplotlib to plot the recovered temperature profile.")
