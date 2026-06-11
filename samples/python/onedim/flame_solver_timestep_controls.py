"""
Timestep Controls for 1D Flame Solvers
======================================

This example demonstrates four solver controls that affect pseudo-transient
convergence in one-dimensional flame calculations:

1. ``time_step_growth_strategy``, which controls when the timestep grows after a
   successful transient step.
2. ``time_step_growth_factor``, which controls how much the timestep grows when
   growth is triggered.
3. ``time_step_regrid``, which controls how many times the solver may refine the
   grid and retry after timestepping is exhausted.
4. ``max_time_step_count``, which controls the number of transient steps allowed
   before a timestep attempt gives up.

The first comparison uses a burner-stabilized premixed flame. The following two
comparisons use a high-pressure counterflow diffusion flame, where solver convergence
can depend on both regridding and the transient-step budget.

Requires: cantera >= 4.0; matplotlib >= 2.0

.. tags:: Python, combustion, 1D flow, premixed flame, diffusion flame,
          strained flame, plotting
"""

from __future__ import annotations

import time

import matplotlib.pyplot as plt
plt.rcParams['figure.constrained_layout.use'] = True
import numpy as np

import cantera as ct


# %%
# Timestep Growth Strategies
# --------------------------
#
# The one-dimensional solver falls back to pseudo-transient timestepping when the Newton
# solver cannot converge directly to steady state. After a successful transient step,
# Cantera can either grow the timestep by a fixed factor or apply one of several
# heuristics that wait for evidence that the larger step is helpful.
#
# This first problem uses a premixed methane / air burner-stabilized flame. The physical
# solution is the same for each case, so the differences in the table are from the
# timestep-growth rule. We test several combinations of growth strategies and growth
# factors.
growth_strategies = (
    ("fixed-growth", 1.5),  # default growth strategy and factor
    ("fixed-growth", 2.0),
    ("steady-norm", 2.0),
    ("transient-residual", 2.0),
    ("residual-ratio", 2.0),
    ("newton-iterations", 2.0),
    ("steady-norm", 6.0),
    ("newton-iterations", 6.0),
)

def make_growth_flame() -> ct.BurnerFlame:
    gas = ct.Solution("gri30.yaml")
    gas.set_equivalence_ratio(0.8, "CH4", {"O2": 1.0, "N2": 3.76})
    gas.TP = 300.0, ct.one_atm

    flame = ct.BurnerFlame(gas, width=0.02)
    flame.burner.mdot = 0.09
    flame.set_refine_criteria(ratio=3.0, slope=0.25, curve=0.35, prune=0.05)
    return flame


def run_growth_strategy(strategy: str, growth_factor: float) -> dict[str, object]:
    flame = make_growth_flame()
    timesteps = []
    flame.set_time_step_callback(
        lambda dt: timesteps.append(float(dt)) or 0
    )
    flame.time_step_growth_factor = growth_factor
    flame.time_step_growth_strategy = strategy

    tic = time.perf_counter()
    print("** Running growth strategy:", strategy)
    flame.solve(loglevel=1, refine_grid=True)
    elapsed = time.perf_counter() - tic

    return {
        "strategy": f"{strategy} (growth factor = {growth_factor:.1f})",
        "n_steps": len(timesteps),
        "dt_sum": float(np.sum(timesteps)),
        "dt_min": float(np.min(timesteps)),
        "dt_max": float(np.max(timesteps)),
        "jacobians": int(sum(flame.solver_stats["jacobian_evals"])),
        "T_max": float(np.max(flame.T)),
        "wall_time": elapsed,
        "dts": np.asarray(timesteps, dtype=float),
    }


growth_results = [run_growth_strategy(strategy, growth_factor)
                  for strategy, growth_factor in growth_strategies]

# %%
# The timestep history is collected with ``set_time_step_callback``.
# The columns below report the number of timesteps, the total transient
# time covered by those steps, the largest timestep, the number of
# Jacobian evaluations, and the maximum temperature of the converged flame.

growth_header = (
    f"{'Strategy':<42} {'Steps':>8} {'sum(dt) [s]':>12} "
    f"{'max(dt) [s]':>12} {'Jac':>5} {'T_max [K]':>12} {'Runtime [s]':>9}"
)
print(growth_header)
print("-" * len(growth_header))
for row in growth_results:
    print(
        f"{row['strategy']:<42} {row['n_steps']:>8d} "
        f"{row['dt_sum']:>12.4e} {row['dt_max']:>12.4e} "
        f"{row['jacobians']:>5d} {row['T_max']:>12.4f} {row['wall_time']:>9.2f}"
    )

# %%
# In short, ``fixed-growth`` always applies ``time_step_growth_factor`` after a
# successful step. ``steady-norm`` and ``transient-residual`` require the
# corresponding residual norm to decrease. ``residual-ratio`` scales the growth
# using the transient residual improvement, capped by
# ``time_step_growth_factor``. ``newton-iterations`` grows the step only after a
# low-iteration Newton solve. Setting ``time_step_growth_factor = 1.0`` disables
# growth.

# %%
# Timestep Histories
# ------------------
#
# The timestep histories show how the strategy changes the sequence of transient
# steps even when all cases reach the same steady fixed-temperature solution.
cmap = plt.get_cmap('Dark2')

fig, ax = plt.subplots()
for i,row in enumerate(growth_results):
    steps = np.arange(1, row["n_steps"] + 1)
    ax.semilogy(
        steps + 0.1*i, row["dts"], "o-", linewidth=1.6, markersize=4,
        color=cmap(i), label=row["strategy"]
    )
ax.set_xlabel("timestep number")
ax.set_ylabel("dt [s]")
ax.grid(True, which="both", alpha=0.3)
ax.legend(fontsize=8)
plt.show()


# %%
# Regridding After Timestep Failure
# ---------------------------------
#
# The next comparison uses a high-pressure hydrogen / oxygen counterflow
# diffusion flame. The initial grid has only 20 points across a 30 mm domain, which is
# too coarse for the steady-state solver to converge. The ``time_step_regrid`` option
# lets the solver refine the grid and retry when a timestepping sequence has reached
# ``max_time_step_count``.
diffusion_width = 30e-3
diffusion_initial_points = 20
diffusion_fuel_mdot = 0.3
diffusion_oxidizer_mdot_factor = 3.0


def make_regrid_flame(
    pressure: float,
    regrid_max: int,
    max_time_step_count: int = 200,
    growth_strategy: str = "fixed-growth",
    growth_factor: float = 1.5,
) -> ct.CounterflowDiffusionFlame:
    gas = ct.Solution("h2o2.yaml")
    flame = ct.CounterflowDiffusionFlame(
        gas, grid=np.linspace(0.0, diffusion_width, diffusion_initial_points)
    )
    flame.max_time_step_count = max_time_step_count
    flame.set_refine_criteria(ratio=2.0, slope=0.06, curve=0.08, prune=0.02)
    flame.time_step_regrid = regrid_max
    flame.time_step_growth_strategy = growth_strategy
    flame.time_step_growth_factor = growth_factor

    flame.P = pressure
    flame.fuel_inlet.X = "H2:1.0"
    flame.fuel_inlet.T = 800.0
    flame.oxidizer_inlet.X = "O2:1.0"
    flame.oxidizer_inlet.T = 711.0

    rho_fuel = flame.fuel_inlet.phase.density
    rho_oxidizer = flame.oxidizer_inlet.phase.density
    flame.fuel_inlet.mdot = diffusion_fuel_mdot
    flame.oxidizer_inlet.mdot = (
        diffusion_fuel_mdot / rho_fuel * rho_oxidizer * diffusion_oxidizer_mdot_factor
    )
    return flame


def run_regrid_case(case: dict[str, object]) -> dict[str, object]:
    flame = make_regrid_flame(
        pressure=case["pressure"],
        regrid_max=case["regrid_max"],
        max_time_step_count=case["max_time_step_count"],
        growth_strategy=case["growth_strategy"],
        growth_factor=case["growth_factor"],
    )
    flame.clear_stats()

    tic = time.perf_counter()
    try:
        flame.solve(loglevel=0, auto=False)
        success = True
    except ct.CanteraError as err:
        success = False
    elapsed = time.perf_counter() - tic

    return {
        "label": case["label"],
        "plot_label": case["plot_label"],
        "pressure": case["pressure"],
        "success": success,
        "grid_points": len(flame.grid),
        "timesteps": int(sum(flame.solver_stats["steps"])),
        "jacobians": int(sum(flame.solver_stats["jacobian_evals"])),
        "peak_temperature": float(np.max(flame.T)) if success else None,
        "grid": np.array(flame.grid, copy=True),
        "temperature": np.array(flame.T, copy=True),
        "wall_time": elapsed,
        "max_time_step_count": case["max_time_step_count"],
        "regrid_max": case["regrid_max"],
        "growth_strategy": case["growth_strategy"],
        "growth_factor": case["growth_factor"],
    }


def print_regrid_rows(results: list[dict[str, object]], label_width: int = 28) -> None:
    header = (
        f"{'Mode':<{label_width}} {'Success':>7} {'Grid':>6} "
        f"{'Steps':>7} {'Jac':>5} {'Tmax [K]':>10} {'Wall [s]':>9}"
    )
    print(header)
    print("-" * len(header))
    for row in results:
        peak_temperature = (
            f"{row['peak_temperature']:>10.1f}"
            if row["peak_temperature"] is not None else f"{'--':>10}"
        )
        print(
            f"{row['label']:<{label_width}} {str(row['success']):>7} "
            f"{row['grid_points']:>6d} {row['timesteps']:>7d} "
            f"{row['jacobians']:>5d} {peak_temperature} {row['wall_time']:>9.2f}"
        )


regrid_pressure = 7e6
regrid_cases = (
    {
        "label": "time_step_regrid = 0",
        "plot_label": "0",
        "pressure": regrid_pressure,
        "regrid_max": 0,
        "max_time_step_count": 200,
        "growth_strategy": "fixed-growth",
        "growth_factor": 1.5,
    },
    {
        "label": "time_step_regrid = 1",
        "plot_label": "1",
        "pressure": regrid_pressure,
        "regrid_max": 1,
        "max_time_step_count": 200,
        "growth_strategy": "fixed-growth",
        "growth_factor": 1.5,
    },
    {
        "label": "time_step_regrid = 3",
        "plot_label": "3",
        "pressure": regrid_pressure,
        "regrid_max": 3,
        "max_time_step_count": 200,
        "growth_strategy": "fixed-growth",
        "growth_factor": 1.5,
    },
)

regrid_results = [run_regrid_case(case) for case in regrid_cases]

# %%
# Here the growth-factor and growth-strategy settings are held at their defaults. Only
# the number of regrid retries changes.
print("Regrid retries after timestep failure")
print(
    f"Case: H2/O2 counterflow diffusion flame at {regrid_pressure / 1e6:.1f} MPa "
    f"on an initial {diffusion_initial_points}-point grid"
)
print("Varied option: time_step_regrid")
print()
print_regrid_rows(regrid_results, label_width=43)

# %%
# ``time_step_regrid = 0`` disables retry-after-regrid. Positive values allow up to that
# many grid-refinement retries after timestep failure. If the refinement criteria do not
# change the grid, the retry path aborts because there is no new discretization to try.

# %%
# Regrid Retry Comparison
# ~~~~~~~~~~~~~~~~~~~~~~~
#
# The first two cases exit before reaching a steady solution. With three regrid
# retries, the solver has enough opportunities to add points and continue.
colors = ["#ae2012" if not row["success"] else "#0a9396" for row in regrid_results]
fig, ax = plt.subplots()
bars = ax.bar(
    [row["plot_label"] for row in regrid_results],
    [row["timesteps"] for row in regrid_results],
    color=colors,
)
ax.set_xlabel("time_step_regrid")
ax.set_ylabel("timesteps before exit")
ax.grid(True, axis="y", alpha=0.3)
for bar, row in zip(bars, regrid_results):
    status = "ok" if row["success"] else "fail"
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + 10,
        f"grid={row['grid_points']}\n{status}",
        ha="center",
        va="bottom",
        fontsize=8,
    )
plt.show()

# %%
# Flame Profile
# ~~~~~~~~~~~~~
#
# The successful 7 MPa case is not just a solver-status change: after regridding,
# the solver finds a steady-state solution for a resolved diffusion flame with a hot
# reaction zone in the counterflow domain.
regrid_result = next(
    row for row in regrid_results if row["success"] and row["regrid_max"] == 3
)

fig, ax = plt.subplots()
ax.plot(
    1e3 * regrid_result["grid"],
    regrid_result["temperature"],
    ".-",
)
ax.set_xlabel("distance from fuel inlet [mm]")
ax.set_ylabel("temperature [K]")
ax.grid(True, alpha=0.3)
plt.show()

# %%
# A More Difficult Pressure
# ~~~~~~~~~~~~~~~~~~~~~~~~~
#
# At 15 MPa, the same coarse counterflow problem is more demanding. Regridding is
# enabled for all cases in this section, but a 200-step timestep budget is still
# insufficient. The first comparison shows the effect of increasing
# ``max_time_step_count`` while leaving the growth strategy fixed.
difficult_pressure = 15e6
difficult_budget_cases = (
    {
        "label": "max_time_step_count = 200",
        "plot_label": "200",
        "pressure": difficult_pressure,
        "max_time_step_count": 200,
        "regrid_max": 3,
        "growth_strategy": "fixed-growth",
        "growth_factor": 1.5,
    },
    {
        "label": "max_time_step_count = 1250",
        "plot_label": "1250",
        "pressure": difficult_pressure,
        "max_time_step_count": 1250,
        "regrid_max": 3,
        "growth_strategy": "fixed-growth",
        "growth_factor": 1.5,
    },
)

difficult_budget_results = [run_regrid_case(case) for case in difficult_budget_cases]

print("Solving with a larger timestep budget")
print(
    f"Case: H2/O2 counterflow diffusion flame at "
    f"{difficult_pressure / 1e6:.1f} MPa"
)
print("Fixed settings: time_step_regrid = 3, time_step_growth_strategy = fixed-growth")
print()
print_regrid_rows(difficult_budget_results, label_width=28)

# %%
# Once the timestep budget is large enough enable solution of the steady-state problem,
# the growth strategy can still affect the work required to reach the converged steady
# solution. The fixed-growth row below reuses the successful result from the previous
# table, and the remaining rows rerun the same problem with adaptive growth strategies.
timestep_budget_result = next(
    row for row in difficult_budget_results if row["max_time_step_count"] == 1250
)
difficult_strategy_cases = (
    {
        "label": "residual-ratio",
        "plot_label": "ratio",
        "pressure": difficult_pressure,
        "max_time_step_count": 1250,
        "regrid_max": 3,
        "growth_strategy": "residual-ratio",
        "growth_factor": 1.5,
    },
    {
        "label": "newton-iterations",
        "plot_label": "newton",
        "pressure": difficult_pressure,
        "max_time_step_count": 1250,
        "regrid_max": 3,
        "growth_strategy": "newton-iterations",
        "growth_factor": 1.5,
    },
)

difficult_strategy_results = [
    {
        **timestep_budget_result,
        "label": "fixed-growth",
        "plot_label": "fixed",
    },
    *[run_regrid_case(case) for case in difficult_strategy_cases],
]

print("Growth-strategy sensitivity after recovery")
print(
    "Fixed settings: max_time_step_count = 1250, time_step_regrid = 3, "
    "time_step_growth_factor = 1.5"
)
print()
print_regrid_rows(difficult_strategy_results, label_width=18)

# %%
# ``max_time_step_count`` limits the number of transient steps before a timestep attempt
# gives up. Larger values can allow harder cases to converge, but solving these cases
# will be slow. On this 15 MPa case, adaptive growth strategies can reduce wall time
# once the timestep budget is large enough.
