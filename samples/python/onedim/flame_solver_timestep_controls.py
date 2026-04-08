# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

"""
Timestep Controls for 1D Flame Solvers
======================================

This example prints three compact comparisons for one-dimensional flame solvers:

1. Named strategies for growing the timestep after a successful transient step.
2. Regridding retries after timestepping is exhausted on a difficult problem.
3. A harder 15 MPa case that is sensitive to both timestep budget and
   timestep-growth strategy.

The first table varies the timestep-growth strategy on a fixed-temperature free
flame. The second varies the regrid retry count on a coarse high-pressure
counterflow diffusion flame. The third shows a numerically difficult
variant where both `max_time_step_count` and `time_step_growth_strategy`
matter.

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
    },
    {
        "label": "time_step_regrid = 1",
        "plot_label": "1",
        "regrid_max": 1,
        "refine_criteria": REGRID_REFINE,
    },
    {
        "label": "time_step_regrid = 3",
        "plot_label": "3",
        "regrid_max": 3,
        "refine_criteria": REGRID_REFINE,
    },
)

DIFFICULT_PRESSURE = 15e6
DIFFICULT_BUDGET_CASES = (
    {
        "label": "max_time_step_count = 200",
        "plot_label": "200",
        "max_time_step_count": 200,
        "regrid_max": 3,
        "growth_strategy": "fixed-growth",
        "growth_factor": 1.5,
    },
    {
        "label": "max_time_step_count = 1250",
        "plot_label": "1250",
        "max_time_step_count": 1250,
        "regrid_max": 3,
        "growth_strategy": "fixed-growth",
        "growth_factor": 1.5,
    },
)

DIFFICULT_STRATEGY_CASES = (
    {
        "label": "residual-ratio",
        "plot_label": "ratio",
        "max_time_step_count": 1250,
        "regrid_max": 3,
        "growth_strategy": "residual-ratio",
        "growth_factor": 1.5,
    },
    {
        "label": "newton-iterations",
        "plot_label": "newton",
        "max_time_step_count": 1250,
        "regrid_max": 3,
        "growth_strategy": "newton-iterations",
        "growth_factor": 1.5,
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
                      refine_criteria: tuple[float, float, float, float],
                      growth_strategy: str = "fixed-growth",
                      growth_factor: float = 1.5,
                      pressure: float = REGRID_PRESSURE,
                      max_time_step_count: int = 200
                      ) -> ct.CounterflowDiffusionFlame:
    gas = ct.Solution("h2o2.yaml")
    flame = ct.CounterflowDiffusionFlame(
        gas, grid=np.linspace(0.0, REGRID_WIDTH, REGRID_INITIAL_POINTS))
    flame.max_time_step_count = max_time_step_count
    flame.set_refine_criteria(
        ratio=refine_criteria[0], slope=refine_criteria[1],
        curve=refine_criteria[2], prune=refine_criteria[3]
    )
    flame.time_step_regrid = regrid_max
    flame.time_step_growth_strategy = growth_strategy
    flame.time_step_growth_factor = growth_factor

    flame.P = pressure
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
        "success": success,
        "grid_points": len(flame.grid),
        "timesteps": int(sum(flame.time_step_stats)),
        "jacobians": int(sum(flame.jacobian_count_stats)),
        "peak_temperature": float(np.max(flame.T)) if success else None,
        "wall_time": elapsed,
    }


def run_difficult_case(case: dict[str, object]) -> dict[str, object]:
    flame = make_regrid_flame(
        case["regrid_max"], REGRID_REFINE, case["growth_strategy"],
        case["growth_factor"], DIFFICULT_PRESSURE, case["max_time_step_count"])
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
        "growth_strategy": case["growth_strategy"],
        "growth_factor": case["growth_factor"],
        "max_time_step_count": case["max_time_step_count"],
        "regrid_max": case["regrid_max"],
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
    print("  Varied option: `time_step_growth_strategy`")
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

    print("\n  Option reminder:")
    print("  `fixed-growth`: always apply `time_step_growth_factor` after a successful step.")
    print("  `steady-norm`: grow only if the steady-state residual norm decreases.")
    print("  `transient-residual`: grow only if the transient residual norm decreases.")
    print("  `residual-ratio`: scale growth based on transient residual improvement.")
    print("  `newton-iterations`: grow only after small-iteration Newton solves.")
    print("  Set `time_step_growth_factor = 1.0` to disable successful-step growth.")


def print_regrid_summary(results: list[dict[str, object]]) -> None:
    print("\nRegrid retries after timestep failure")
    print(
        "  Case: H2/O2 counterflow diffusion flame at "
        f"{REGRID_PRESSURE / 1e6:.1f} MPa on an initial {REGRID_INITIAL_POINTS}-point grid"
    )
    print("  Varied option: `time_step_regrid`")
    print("  Growth-factor and growth-strategy settings are left at their defaults.")

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

    print("\n  Option reminder:")
    print("  `time_step_regrid = 0`: disable retry-after-regrid.")
    print("  `time_step_regrid = n > 0`: allow up to `n` regrid retries after timestep failure.")
    print("  If refinement criteria do not change the grid, the retry path aborts.")


def print_difficult_summary(budget_results: list[dict[str, object]],
                            strategy_results: list[dict[str, object]]) -> None:
    print("\nA numerically difficult 15 MPa case")
    print(
        "  Case: the same coarse H2/O2 counterflow diffusion flame at "
        f"{DIFFICULT_PRESSURE / 1e6:.1f} MPa"
    )
    print("  This case converges only after allowing more transient steps,")
    print("  and then shows meaningful sensitivity to the growth heuristic.")

    header = (
        f"{'Mode':<25} {'Success':>7} {'Grid':>6} "
        f"{'Steps':>7} {'Jac':>5} {'Tmax [K]':>10} {'Wall [s]':>9}"
    )
    print("\n  Recovery with a larger timestep budget")
    print("  Varied option: `max_time_step_count` with `time_step_regrid = 3`")
    print("  Growth settings are held at `fixed-growth` with factor 1.5.")
    print(f"\n  {header}")
    print(f"  {'-' * len(header)}")
    for row in budget_results:
        peak_temperature = (
            f"{row['peak_temperature']:>10.1f}"
            if row["peak_temperature"] is not None else f"{'--':>10}"
        )
        print(
            "  "
            f"{row['label']:<25} {str(row['success']):>7} "
            f"{row['grid_points']:>6d} {row['timesteps']:>7d} "
            f"{row['jacobians']:>5d} {peak_temperature} {row['wall_time']:>9.2f}"
        )

    strategy_header = (
        f"{'Strategy':<18} {'Success':>7} {'Grid':>6} "
        f"{'Steps':>7} {'Jac':>5} {'Tmax [K]':>10} {'Wall [s]':>9}"
    )
    print("\n  Growth-strategy sensitivity after recovery")
    print(
        "  Varied option: `time_step_growth_strategy` with "
        "`max_time_step_count = 1250`, `time_step_regrid = 3`, and "
        "`time_step_growth_factor = 1.5`."
    )
    print(f"\n  {strategy_header}")
    print(f"  {'-' * len(strategy_header)}")
    for row in strategy_results:
        peak_temperature = (
            f"{row['peak_temperature']:>10.1f}"
            if row["peak_temperature"] is not None else f"{'--':>10}"
        )
        print(
            "  "
            f"{row['growth_strategy']:<18} {str(row['success']):>7} "
            f"{row['grid_points']:>6d} {row['timesteps']:>7d} "
            f"{row['jacobians']:>5d} {peak_temperature} {row['wall_time']:>9.2f}"
        )

    print("\n  Option reminder:")
    print("  `max_time_step_count` limits accepted transient steps before a timestep attempt gives up.")
    print("  Larger values can recover harder cases, but they increase wall time substantially.")
    print("  On this 15 MPa case, adaptive growth strategies can also cut wall time substantially once the timestep budget is large enough.")


def plot_results(growth_results: list[dict[str, object]],
                 regrid_results: list[dict[str, object]],
                 difficult_strategy_results: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))

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

    colors = [
        growth_colors.get(row["growth_strategy"], "#0a9396")
        if row["success"] else "#ae2012"
        for row in difficult_strategy_results
    ]
    bars = axes[2].bar(
        [row["plot_label"] for row in difficult_strategy_results],
        [row["wall_time"] for row in difficult_strategy_results],
        color=colors,
    )
    axes[2].set_title("15 MPa heuristic sensitivity")
    axes[2].set_xlabel("time_step_growth_strategy")
    axes[2].set_ylabel("Wall time [s]")
    axes[2].grid(True, axis="y", alpha=0.3)
    for bar, row in zip(bars, difficult_strategy_results):
        tag = "ok" if row["success"] else "fail"
        axes[2].text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.3,
            f"jac={row['jacobians']}\n{tag}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    fig.tight_layout()
    if "agg" not in plt.get_backend().lower():
        plt.show()


def main() -> None:
    growth_results = [run_growth_case(strategy) for strategy in GROWTH_STRATEGIES]
    regrid_results = [run_regrid_case(case) for case in REGRID_CASES]
    difficult_budget_results = [run_difficult_case(case) for case in DIFFICULT_BUDGET_CASES]
    recovered_baseline = next(
        row for row in difficult_budget_results if row["max_time_step_count"] == 1250
    )
    difficult_strategy_results = [
        {
            **recovered_baseline,
            "label": "fixed-growth",
            "plot_label": "fixed",
        },
        *[run_difficult_case(case) for case in DIFFICULT_STRATEGY_CASES],
    ]

    print_growth_summary(growth_results)
    print_regrid_summary(regrid_results)
    print_difficult_summary(difficult_budget_results, difficult_strategy_results)

    if plt is not None:
        plot_results(
            growth_results, regrid_results, difficult_strategy_results)
    else:
        print("\nInstall matplotlib to plot the timestep histories and comparison tables.")


if __name__ == "__main__":
    main()
