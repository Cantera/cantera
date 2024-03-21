r"""
Flame Speed with Convergence Analysis
=====================================

Requires: cantera >= 3.0.0, matplotlib >= 2.0, pandas

In this example we simulate a freely-propagating, adiabatic, 1-D flame and

* Calculate its laminar burning velocity
* Estimate the uncertainty in the laminar burning velocity calculation due to grid size.

.. tags:: Python, combustion, 1D flow, flame speed, premixed flame, plotting

The figure below illustrates the setup, in a flame-fixed co-ordinate system. The
reactants enter with density :math:`\rho_u`, temperature :math:`T_u` and speed
:math:`S_u`. The products exit the flame at speed :math:`S_b`, density :math:`\rho_b`,
and temperature :math:`T_b`.

.. image:: /_static/images/samples/flame-speed.svg
   :width: 50%
   :alt: Freely Propagating Flame
   :align: center
"""

# %%
# Import Modules
# --------------

import cantera as ct
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib
import scipy
import scipy.optimize

# %%
# Define plotting preference
# --------------------------

plt.style.use("ggplot")
plt.style.use("seaborn-v0_8-deep")
plt.rcParams["figure.constrained_layout.use"] = True

# %%
# Estimate uncertainty from grid size and speeds
# ----------------------------------------------

def extrapolate_uncertainty(grids, speeds, plot=True):
    """
    Given a list of grid sizes and a corresponding list of flame speeds,
    extrapolate and estimate the uncertainty in the final flame speed.
    Also makes a plot, unless called with `plot=False`.
    """
    grids = list(grids)
    speeds = list(speeds)

    def speed_from_grid_size(grid_size, true_speed, error):
        """
        Given a grid size (or an array or list of grid sizes)
        return a prediction (or array of predictions)
        of the computed flame speed, based on
        the parameters `true_speed` and `error`.

        It seems, from experience, that error scales roughly with
        1/grid_size, so we assume that form.
        """
        return true_speed + error / np.array(grid_size)

    # Fit the chosen form of speed_from_grid_size, to the last four
    # speed and grid size values.
    popt, pcov = scipy.optimize.curve_fit(speed_from_grid_size, grids[-4:], speeds[-4:])

    # How bad the fit was gives you some error, `percent_error_in_true_speed`.
    perr = np.sqrt(np.diag(pcov))
    true_speed_estimate = popt[0]
    percent_error_in_true_speed = perr[0] / popt[0]
    print(
        f"Fitted true_speed is {popt[0] * 100:.4f} ± {perr[0] * 100:.4f} cm/s "
        f"({percent_error_in_true_speed:.1%})"
    )

    # How far your extrapolated infinite grid value is from your extrapolated
    # (or interpolated) final grid value, gives you some other error, `estimated_percent_error`
    estimated_percent_error = (
        speed_from_grid_size(grids[-1], *popt) - true_speed_estimate
    ) / true_speed_estimate
    print(f"Estimated error in final calculation {estimated_percent_error:.1%}")

    # The total estimated error is the sum of these two errors.
    total_percent_error_estimate = abs(percent_error_in_true_speed) + abs(
        estimated_percent_error
    )
    print(f"Estimated total error {total_percent_error_estimate:.1%}")

    if plot:
        fig, ax = plt.subplots()
        ax.semilogx(grids, speeds, "o-")
        ax.set_ylim(
            min(speeds[-5:] + [true_speed_estimate - perr[0]]) * 0.95,
            max(speeds[-5:] + [true_speed_estimate + perr[0]]) * 1.05,
        )
        ax.plot(grids[-4:], speeds[-4:], "or")
        extrapolated_grids = grids + [grids[-1] * i for i in range(2, 8)]
        ax.plot(
            extrapolated_grids, speed_from_grid_size(extrapolated_grids, *popt), ":r"
        )
        ax.set_xlim(*ax.get_xlim())  # Prevent automatic expansion of axis limits
        ax.hlines(true_speed_estimate, *ax.get_xlim(), colors="r", linestyles="dashed")
        ax.hlines(
            true_speed_estimate + perr[0],
            *ax.get_xlim(),
            colors="r",
            linestyles="dashed",
            alpha=0.3,
        )
        ax.hlines(
            true_speed_estimate - perr[0],
            *ax.get_xlim(),
            colors="r",
            linestyles="dashed",
            alpha=0.3,
        )
        ax.fill_between(
            ax.get_xlim(),
            true_speed_estimate - perr[0],
            true_speed_estimate + perr[0],
            facecolor="red",
            alpha=0.1,
        )

        above = popt[1] / abs(
            popt[1]
        )  # will be +1 if approach from above or -1 if approach from below

        ax.annotate(
            "",
            xy=(grids[-1], true_speed_estimate),
            xycoords="data",
            xytext=(grids[-1], speed_from_grid_size(grids[-1], *popt)),
            textcoords="data",
            arrowprops=dict(
                arrowstyle="|-|, widthA=0.5, widthB=0.5",
                linewidth=1,
                connectionstyle="arc3",
                color="black",
                shrinkA=0,
                shrinkB=0,
            ),
        )

        ax.annotate(
            f"{abs(estimated_percent_error):.1%}",
            xy=(grids[-1], speed_from_grid_size(grids[-1], *popt)),
            xycoords="data",
            xytext=(5, 15 * above),
            va="center",
            textcoords="offset points",
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),
        )

        ax.annotate(
            "",
            xy=(grids[-1] * 4, true_speed_estimate - (above * perr[0])),
            xycoords="data",
            xytext=(grids[-1] * 4, true_speed_estimate),
            textcoords="data",
            arrowprops=dict(
                arrowstyle="|-|, widthA=0.5, widthB=0.5",
                linewidth=1,
                connectionstyle="arc3",
                color="black",
                shrinkA=0,
                shrinkB=0,
            ),
        )
        ax.annotate(
            f"{abs(percent_error_in_true_speed):.1%}",
            xy=(grids[-1] * 4, true_speed_estimate - (above * perr[0])),
            xycoords="data",
            xytext=(5, -15 * above),
            va="center",
            textcoords="offset points",
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),
        )

        ax.set(xlabel="Grid size", ylabel="Flame speed (m/s)")

    return true_speed_estimate, total_percent_error_estimate


# %%
def make_callback(flame):
    """
    Create and return a callback function that you will attach to
    a flame solver. The reason we define a function to make the callback function,
    instead of just defining the callback function, is so that it can store
    a pair of lists that persist between function calls, to store the
    values of grid size and flame speed.

    This factory returns the callback function, and the two lists:
    (callback, speeds, grids)
    """
    speeds = []
    grids = []

    def callback(_):
        speed = flame.velocity[0]
        grid = len(flame.grid)
        speeds.append(speed)
        grids.append(grid)
        print(f"Iteration {len(grids)}")
        print(f"Current flame speed is is {speed * 100:.4f} cm/s")
        if len(grids) < 5:
            return 1.0  #
        try:
            extrapolate_uncertainty(grids, speeds)
        except Exception as e:
            print("Couldn't estimate uncertainty. " + str(e))
            return 1.0  # continue anyway
        return 1.0

    return callback, speeds, grids


# %%
# Define the reactant conditions, gas mixture and kinetic mechanism associated with the gas
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Inlet Temperature in Kelvin and Inlet Pressure in Pascals
# In this case we are setting the inlet T and P to room temperature conditions
To = 300
Po = 101325

# Define the gas mixture and kinetics
# In this case, we are choosing a GRI3.0 gas
gas = ct.Solution("gri30.yaml")

# Create a stoichiometric CH4/Air premixed mixture
gas.set_equivalence_ratio(1.0, "CH4", {"O2": 1.0, "N2": 3.76})
gas.TP = To, Po

# %%
# Define flame simulation conditions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Domain width in metres
width = 0.014

# Create the flame object
flame = ct.FreeFlame(gas, width=width)

# Define logging level
loglevel = 1

# Define tight tolerances for the solver
refine_criteria = {"ratio": 2, "slope": 0.01, "curve": 0.01}
flame.set_refine_criteria(**refine_criteria)

# Set maximum number of grid points to be very high (otherwise default is 1000)
flame.set_max_grid_points(flame.domains[flame.domain_index("flame")], 1e4)

# Set up the the callback function and lists of speeds and grids
callback, speeds, grids = make_callback(flame)
flame.set_steady_callback(callback)

# %%
# Solve
# ~~~~~
#
# After the first five iterations, it will start to estimate the uncertainty.

flame.solve(loglevel=loglevel, auto=True)

Su0 = flame.velocity[0]
print(f"Flame Speed is: {Su0 * 100:.2f} cm/s")

# %%
# Use the final lists of grid sizes and flame speeds to make one final extrapolation
# "best guess"

best_true_speed_estimate, best_total_percent_error_estimate = extrapolate_uncertainty(
    grids, speeds
)

best_true_speed_estimate


# %%
# Analyze the error predictions
# -----------------------------
#
# Now let's see how good our error estimates were, with hindsight.
#
# If we assume that the final answer, with a very fine grid, has actually converged and
# is is the "truth", then we can find out how large the errors were in the previous
# values, and compare these with our estimated errors. This will show if our estimates
# are reasonable, or conservative, or too optimistic.

def analyze_errors(grids, speeds, true_speed):
    """
    If we assume that the final answer, with a very fine grid,
    has actually converged and is is the "truth", then we can
    find out how large the errors were in the previous values,
    and compare these with our estimated errors.
    This will show if our estimates are reasonable, or conservative, or too optimistic.
    """
    true_speed_estimates = np.full_like(speeds, np.NaN)
    total_percent_error_estimates = np.full_like(speeds, np.NaN)
    actual_extrapolated_percent_errors = np.full_like(speeds, np.NaN)
    actual_raw_percent_errors = np.full_like(speeds, np.NaN)
    for i in range(3, len(grids)):
        print(grids[: i + 1])
        true_speed_estimate, total_percent_error_estimate = extrapolate_uncertainty(
            grids[: i + 1], speeds[: i + 1], plot=False
        )
        actual_extrapolated_percent_error = (
            abs(true_speed_estimate - true_speed) / true_speed
        )
        actual_raw_percent_error = abs(speeds[i] - true_speed) / true_speed
        print(
            "Actual extrapolated error (with hindsight) "
            f"{actual_extrapolated_percent_error:.1%}"
        )
        print(f"Actual raw error (with hindsight) {actual_raw_percent_error:.1%}")

        true_speed_estimates[i] = true_speed_estimate
        total_percent_error_estimates[i] = total_percent_error_estimate
        actual_extrapolated_percent_errors[i] = actual_extrapolated_percent_error
        actual_raw_percent_errors[i] = actual_raw_percent_error
        print()

    fig, ax = plt.subplots()
    ax.loglog(grids, actual_raw_percent_errors * 100, "o-", label="raw error")
    ax.loglog(
        grids,
        actual_extrapolated_percent_errors * 100,
        "o-",
        label="extrapolated error",
    )
    ax.loglog(
        grids, total_percent_error_estimates * 100, "o-", label="estimated error"
    )
    ax.set(xlabel="Grid size", ylabel="Error in flame speed (%)")
    ax.legend()
    ax.set_title(flame.get_refine_criteria())
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.PercentFormatter())
    flame.get_refine_criteria()

    data = pd.DataFrame(
        data={
            "actual error in raw value": actual_raw_percent_errors * 100,
            "actual error in extrapolated value": actual_extrapolated_percent_errors
            * 100,
            "estimated error": total_percent_error_estimates * 100,
        },
        index=grids,
    )
    return data


analyze_errors(grids, speeds, best_true_speed_estimate)

# %%
# Repeat with less tight refine criteria
# --------------------------------------

refine_criteria = {"ratio": 3, "slope": 0.1, "curve": 0.1}

# Reset the gas
gas.set_equivalence_ratio(1.0, "CH4", {"O2": 1.0, "N2": 3.76})
gas.TP = To, Po

# Create a new flame object
flame = ct.FreeFlame(gas, width=width)

flame.set_refine_criteria(**refine_criteria)
flame.set_max_grid_points(flame.domains[flame.domain_index("flame")], 1e4)

callback, speeds, grids = make_callback(flame)
flame.set_steady_callback(callback)

# Define logging level
loglevel = 1

flame.solve(loglevel=loglevel, auto=True)

Su0 = flame.velocity[0]
print(f"Flame Speed is: {Su0 * 100:.2f} cm/s")

# Use the best true speed estimate from the fine grid tight criteria above
analyze_errors(grids, speeds, best_true_speed_estimate)

# %%
# Default (loose) criteria
# ------------------------

flame = ct.FreeFlame(gas, width=width)
refine_criteria = flame.get_refine_criteria()
refine_criteria.update({"prune": 0})
refine_criteria

gas.set_equivalence_ratio(1.0, "CH4", {"O2": 1.0, "N2": 3.76})
gas.TP = To, Po

# Create a new flame object
flame = ct.FreeFlame(gas, width=width)

flame.set_refine_criteria(**refine_criteria)
flame.set_max_grid_points(flame.domains[flame.domain_index("flame")], 1e4)

callback, speeds, grids = make_callback(flame)
flame.set_steady_callback(callback)

# Define logging level
loglevel = 1

flame.solve(loglevel=loglevel, auto=True)

Su0 = flame.velocity[0]
print(f"Flame Speed is: {Su0 * 100:.2f} cm/s")

analyze_errors(grids, speeds, best_true_speed_estimate)

# %%
# Middling refine criteria
# ------------------------

refine_criteria = {"ratio": 3, "slope": 0.1, "curve": 0.1}

# Reset the gas
gas.set_equivalence_ratio(1.0, "CH4", {"O2": 1.0, "N2": 3.76})
gas.TP = To, Po

# Create a new flame object
flame = ct.FreeFlame(gas, width=width)

flame.set_refine_criteria(**refine_criteria)
flame.set_max_grid_points(flame.domains[flame.domain_index("flame")], 1e4)

callback, speeds, grids = make_callback(flame)
flame.set_steady_callback(callback)

# Define logging level
loglevel = 1

flame.solve(loglevel=loglevel, auto=True)

Su0 = flame.velocity[0]
print(f"Flame Speed is: {Su0 * 100:.2f} cm/s")

# %%
analyze_errors(grids, speeds, best_true_speed_estimate)

# %%
# Try a Hydrogen flame (still with GRI mech)
# ------------------------------------------

# Tight criteria
refine_criteria = {"ratio": 2, "slope": 0.01, "curve": 0.01}

# Reset the gas
gas.set_equivalence_ratio(1.0, "H2", {"O2": 1.0, "N2": 3.76})
gas.TP = To, Po

# Create a new flame object
flame = ct.FreeFlame(gas, width=width)

flame.set_refine_criteria(**refine_criteria)
flame.set_max_grid_points(flame.domains[flame.domain_index("flame")], 1e4)

callback, speeds, grids = make_callback(flame)
flame.set_steady_callback(callback)

# Define logging level
loglevel = 1

flame.solve(loglevel=loglevel, auto=True)

Su0 = flame.velocity[0]
print(f"Flame Speed is: {Su0 * 100:.2f} cm/s")

# %%
# get a new best true speed estimate
best_true_speed_estimate, best_total_percent_error_estimate = extrapolate_uncertainty(
    grids, speeds
)

# %%
analyze_errors(grids, speeds, best_true_speed_estimate)

# %%
# Middling refine criteria, Hydrogen flame
# ----------------------------------------

refine_criteria = {"ratio": 3, "slope": 0.1, "curve": 0.1}

# Reset the gas
gas.set_equivalence_ratio(1.0, "H2", {"O2": 1.0, "N2": 3.76})
gas.TP = To, Po

# Create a new flame object
flame = ct.FreeFlame(gas, width=width)

flame.set_refine_criteria(**refine_criteria)
flame.set_max_grid_points(flame.domains[flame.domain_index("flame")], 1e4)

callback, speeds, grids = make_callback(flame)
flame.set_steady_callback(callback)

# Define logging level
loglevel = 1

flame.solve(loglevel=loglevel, auto=True)

Su0 = flame.velocity[0]
print(f"Flame Speed is: {Su0 * 100:.2f} cm/s")

# %%
analyze_errors(grids, speeds, best_true_speed_estimate)
