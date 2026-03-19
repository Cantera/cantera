# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

"""
Timestep Controls for 1D Flame Solvers
======================================

This example compares two sets of solver controls for one-dimensional flames:

1. Named strategies for growing the timestep after a successful transient step.
2. Regridding retries after timestepping is exhausted on a difficult problem.

The growth-strategy comparison uses a fixed-temperature free flame so the
accepted timestep history can be compared directly at low cost. The regridding
comparison uses a deliberately coarse 7 MPa hydrogen / oxygen counterflow
diffusion flame that fails without recovery.

Requires: cantera >= 4.0; matplotlib is optional for plotting

.. tags:: Python, combustion, 1D flow, premixed flame, diffusion flame, strained flame,
   plotting
"""

from __future__ import annotations

import time

import numpy as np

import cantera as ct

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


GROWTH_CASE = {
    "Tin": 250.0,
    "pressure": ct.one_atm,
    "reactants": "H2:3, O2:1, AR:10",
    "width": 20e-3,
}
GROWTH_FACTOR = 2.0
GROWTH_STRATEGIES = (
    "fixed-growth",
    "steady-norm",
    "transient-residual",
    "residual-ratio",
    "newton-iterations",
)

REGRID_WIDTH = 30e-3
REGRID_INITIAL_POINTS = 20
REGRID_PRESSURE = 7e6
REGRID_REFINE = (2.0, 0.06, 0.08, 0.02)
REGRID_CASES = (
    {
        "label": "time_step_regrid = 0",
        "plot_label": "0",
        "regrid_max": 0,
        "refine_criteria": REGRID_REFINE,
        "description": "No retry after the first exhausted timestepping sequence.",
    },
    {
        "label": "time_step_regrid = 1",
        "plot_label": "1",
        "regrid_max": 1,
        "refine_criteria": REGRID_REFINE,
        "description": "One retry is enough to change the grid, but not enough to converge.",
    },
    {
        "label": "time_step_regrid = 3",
        "plot_label": "3",
        "regrid_max": 3,
        "refine_criteria": REGRID_REFINE,
        "description": "Three retries let the solver refine repeatedly and recover.",
    },
    {
        "label": "time_step_regrid = 3 with a frozen grid",
        "plot_label": "3\nfrozen",
        "regrid_max": 3,
        "refine_criteria": (1000.0, 1.0, 1.0, 0.0),
        "description": "The retry path is entered, but the grid cannot change.",
    },
)


def make_growth_flame() -> ct.FreeFlame:
    gas = ct.Solution("h2o2.yaml")
    gas.TPX = GROWTH_CASE["Tin"], GROWTH_CASE["pressure"], GROWTH_CASE["reactants"]

    flame = ct.FreeFlame(gas, width=GROWTH_CASE["width"])
    flame.flame.set_steady_tolerances(default=[1.0e-5, 1.0e-14])
    flame.flame.set_transient_tolerances(default=[1.0e-4, 1.0e-11])
    flame.inlet.T = GROWTH_CASE["Tin"]
    flame.inlet.X = GROWTH_CASE["reactants"]
    return flame


def run_growth_case(strategy: str) -> dict[str, object]:
    flame = make_growth_flame()
    dts = []
    flame.set_time_step_callback(lambda dt: dts.append(float(dt)) or 0)
    flame.time_step_growth_factor = GROWTH_FACTOR
    flame.time_step_growth_strategy = strategy
    flame.clear_stats()

    flame.energy_enabled = False
    t0 = time.perf_counter()
    flame.solve(loglevel=0, refine_grid=False)
    elapsed = time.perf_counter() - t0

    return {
        "strategy": strategy,
        "accepted_steps": len(dts),
        "dt_sum": float(np.sum(dts)),
        "dt_min": float(np.min(dts)),
        "dt_max": float(np.max(dts)),
        "jacobians": int(sum(flame.jacobian_count_stats)),
        "velocity": float(flame.velocity[0]),
        "wall_time": elapsed,
        "dts": np.asarray(dts, dtype=float),
    }


def make_regrid_flame(regrid_max: int,
                      refine_criteria: tuple[float, float, float, float]
                      ) -> ct.CounterflowDiffusionFlame:
    gas = ct.Solution("h2o2.yaml")
    flame = ct.CounterflowDiffusionFlame(
        gas, grid=np.linspace(0.0, REGRID_WIDTH, REGRID_INITIAL_POINTS))
    flame.max_time_step_count = 200
    flame.set_refine_criteria(
        ratio=refine_criteria[0], slope=refine_criteria[1],
        curve=refine_criteria[2], prune=refine_criteria[3]
    )
    flame.time_step_regrid = regrid_max

    flame.P = REGRID_PRESSURE
    flame.fuel_inlet.X = "H2:1.0"
    flame.fuel_inlet.T = 800.0
    flame.oxidizer_inlet.X = "O2:1.0"
    flame.oxidizer_inlet.T = 711.0

    rho_f = flame.fuel_inlet.phase.density
    rho_o = flame.oxidizer_inlet.phase.density
    flame.fuel_inlet.mdot = 0.3
    flame.oxidizer_inlet.mdot = (0.3 / rho_f) * rho_o * 3.0
    return flame


def run_regrid_case(case: dict[str, object]) -> dict[str, object]:
    flame = make_regrid_flame(case["regrid_max"], case["refine_criteria"])
    flame.clear_stats()

    t0 = time.perf_counter()
    try:
        flame.solve(loglevel=0, auto=False)
        success = True
    except ct.CanteraError:
        success = False
    elapsed = time.perf_counter() - t0

    return {
        "label": case["label"],
        "plot_label": case["plot_label"],
        "description": case["description"],
        "success": success,
        "grid_points": len(flame.grid),
        "timesteps": int(sum(flame.time_step_stats)),
        "jacobians": int(sum(flame.jacobian_count_stats)),
        "peak_temperature": float(np.max(flame.T)) if success else None,
        "wall_time": elapsed,
    }


def print_growth_summary(results: list[dict[str, object]]) -> None:
    print("\nNamed timestep-growth strategies")
    print(
        "  Case: fixed-temperature H2/O2/Ar free flame at "
        f"{GROWTH_CASE['Tin']:.0f} K, {GROWTH_CASE['pressure'] / ct.one_atm:.1f} atm, "
        f"width = {1e3 * GROWTH_CASE['width']:.0f} mm"
    )
    print(f"  Configured time_step_growth_factor = {GROWTH_FACTOR:.1f}")

    header = (
        f"{'Strategy':<20} {'Accepted':>8} {'sum(dt) [s]':>12} "
        f"{'max(dt) [s]':>12} {'Jac':>5} {'u0 [m/s]':>12}"
    )
    print(f"\n  {header}")
    print(f"  {'-' * len(header)}")
    for row in results:
        print(
            "  "
            f"{row['strategy']:<20} {row['accepted_steps']:>8d} "
            f"{row['dt_sum']:>12.4e} {row['dt_max']:>12.4e} "
            f"{row['jacobians']:>5d} {row['velocity']:>12.10f}"
        )

    print("\n  Interpretation:")
    print("  `fixed-growth` is the most aggressive strategy on this case.")
    print("  `steady-norm` and `residual-ratio` are more conservative.")
    print("  `newton-iterations` is the most conservative.")
    print("  `transient-residual` happens to match `fixed-growth` here.")
    print("  Set `time_step_growth_factor = 1.0` to disable successful-step growth.")


def print_regrid_summary(results: list[dict[str, object]]) -> None:
    print("\nRegrid retries after timestep failure")
    print(
        "  Case: H2/O2 counterflow diffusion flame at "
        f"{REGRID_PRESSURE / 1e6:.1f} MPa on an initial {REGRID_INITIAL_POINTS}-point grid"
    )

    header = (
        f"{'Mode':<43} {'Success':>7} {'Grid':>6} "
        f"{'Steps':>7} {'Jac':>5} {'Tmax [K]':>10} {'Wall [s]':>9}"
    )
    print(f"\n  {header}")
    print(f"  {'-' * len(header)}")
    for row in results:
        peak_temperature = (
            f"{row['peak_temperature']:>10.1f}"
            if row["peak_temperature"] is not None else f"{'--':>10}"
        )
        print(
            "  "
            f"{row['label']:<43} {str(row['success']):>7} "
            f"{row['grid_points']:>6d} {row['timesteps']:>7d} "
            f"{row['jacobians']:>5d} {peak_temperature} {row['wall_time']:>9.2f}"
        )

    print("\n  Interpretation:")
    for row in results:
        print(f"  {row['label']}: {row['description']}")


def plot_results(growth_results: list[dict[str, object]],
                 regrid_results: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

    growth_colors = {
        "fixed-growth": "#005f73",
        "steady-norm": "#bb3e03",
        "transient-residual": "#0a9396",
        "residual-ratio": "#ca6702",
        "newton-iterations": "#ae2012",
    }
    for row in growth_results:
        steps = np.arange(1, row["accepted_steps"] + 1)
        axes[0].semilogy(
            steps, row["dts"], "o-", linewidth=1.6, markersize=4,
            color=growth_colors[row["strategy"]], label=row["strategy"]
        )
    axes[0].set_title("Accepted dt histories")
    axes[0].set_xlabel("Accepted timestep number")
    axes[0].set_ylabel("dt [s]")
    axes[0].grid(True, which="both", alpha=0.3)
    axes[0].legend(fontsize=8)

    colors = ["#ae2012" if not row["success"] else "#0a9396" for row in regrid_results]
    bars = axes[1].bar(
        [row["plot_label"] for row in regrid_results],
        [row["timesteps"] for row in regrid_results],
        color=colors,
    )
    axes[1].set_title("Regrid retry modes")
    axes[1].set_xlabel("time_step_regrid")
    axes[1].set_ylabel("Accepted timesteps before exit")
    axes[1].grid(True, axis="y", alpha=0.3)
    for bar, row in zip(bars, regrid_results):
        tag = "ok" if row["success"] else "fail"
        axes[1].text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 10,
            f"grid={row['grid_points']}\n{tag}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    fig.tight_layout()
    plt.show()


def main() -> None:
    growth_results = [run_growth_case(strategy) for strategy in GROWTH_STRATEGIES]
    regrid_results = [run_regrid_case(case) for case in REGRID_CASES]

    print_growth_summary(growth_results)
    print_regrid_summary(regrid_results)

    if plt is not None:
        plot_results(growth_results, regrid_results)
    else:
        print("\nInstall matplotlib to plot the timestep histories and retry comparison.")


if __name__ == "__main__":
    main()
