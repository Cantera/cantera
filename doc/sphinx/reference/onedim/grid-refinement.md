# Grid Refinement for One-dimensional Flows

## Overview

The Cantera 1D solver uses a progressive mesh refinement strategy when solving the
governing equations. An initial coarse grid is solved, and then the solution is
analyzed to identify regions where the solution is not well resolved, such as areas
with steep gradients. Additional grid points are inserted in these regions and then
the solution is solved again. If too many points are inserted, the solver may have
difficulty converging, and so the algorithm prevents itself from adding too many points
during one refinement step. These additional points are added during subsequent
refinement steps. This process is repeated until the solution is deemed to be well
resolved. This incremental refinement procedure maximizes the probability that the
initial guess for each subsequent refined solution lies within the domain of
convergence of the Newton solver.

The algorithm tracks information about which grid points
should be kept and which ones need to have additional grid points inserted to the
right of them. This data is used to refine the grid as well as to keep too many
points from being removed during a refinement step.


## Refinement Criteria

There are 4 refinement criteria that can be set by a user to control the grid
refinement process:

  - ratio
  - slope
  - curve
  - prune

### Ratio

This refinement option controls the maximum ratio of the spacing between the grid
points that are adjacent to a point under consideration.If the ratio between two
adjacent grid points is greater than the ratio
specified by the user ($ \alpha_{max} $), then an additional grid point is inserted
between the two points. Depending on which interval adjacent to a point under
consideration is larger, two conditions are checked:

#### Larger Left Interval

For grid points that have a larger spacing on the left interval, as shown in the figure below,

```{figure} /_static/images/onedim/refinement_large_left_interval.svg
:width: 75%
:alt: Grid with a spacing on the right interval that is too large.
:align: center

Grid with a spacing on the right interval that is too large.
```

the ratio, $ \alpha $ is calculated as:

$$
  \alpha = \frac{{z_j - z_{j-1}}}{{z_{j+1} - z_j}}
$$

#### Larger Right Interval

For grid points that have a larger spacing on the right interval, as shown in the figure below,

```{figure} /_static/images/onedim/refinement_large_right_interval.svg
:width: 75%
:alt: Grid with a spacing on the left interval that is too large.
:align: center

Grid with a spacing on the left interval that is too large.
```

the ratio, $ \alpha $ is calculated as:

$$
  \alpha = \frac{{z_{j+1} - z_j}}{{z_j - z_{j-1}}}
$$


These conditions are checked for all grid points in the domain, and if the ratio
is greater than the user-specified value, then an additional grid point is inserted
between the two points. To illustrate the interplay between the insertion of points
in the domain and the marking of points to be kept, consider again the two cases above.

For the case where the spacing on the right interval is too large, the grid is refined
and a stencil of points around the point under consideration are marked to be kept.
The figure below shows an example of the stencil of points that are marked to be kept
when the spacing on the right interval is too large.

```{figure} /_static/images/onedim/refinement_right_interval_stencil.svg
:width: 75%
:alt: The points that are marked as needing new grid points inserted to the right of
      them and which ones should be marked to be kept.
:align: center

The points that are marked as needing new grid points inserted to the right of them
and which ones should be marked to be kept.
```

Likewise, for the case where the spacing on the left interval is too large, the grid is
refined and a stencil of points around the point under consideration are marked to be kept.
The figure below shows an example of the stencil of points that are marked to be kept
when the spacing on the left interval is too large.

```{figure} /_static/images/onedim/refinement_left_interval_stencil.svg
:width: 75%
:alt: The points that are marked as needing new grid points inserted to the right of
      them and which ones should be marked to be kept.
:align: center

The points that are marked as needing new grid points inserted to the right of them
and which ones should be marked to be kept.
```

Marking the surrounding grid points prevents them from being removed during the refinement
step that they were marked in. This is one of the ways that the algorithm ensures that
the grid does not change too drastically during a refinement step.



### Slope

This refinement option controls the maximum change in a solution component between
two adjacent grid points compared to the maximum change in that component across the
entire domain. If the ratio of the change in a solution component between two adjacent
grid points is greater than the user-specified value, then an additional grid point is
inserted between the two points. This criterion is used to ensure that the solution
is well-resolved in regions where the solution has steep gradients.

The specific criteria that is being enforced is:
$$
  |x_{n,j+1} - x_{n,j}| <= \beta_{max} (\displaystyle\max_j x_n - \displaystyle\min_j x_n)
$$

Where $ \beta_{max} $ is a number between 0 and 1 that is specified by the user, $ x_n $
is the solution component being considered, and $ \displaystyle\max_j x_n $ and $ \displaystyle\min_j x_n $ are the
maximum and minimum values of the solution component across the domain. Visually, these
quantities are shown in the image below.

```{figure} /_static/images/onedim/refinement_slope_metric_plot.svg
:width: 75%
:alt: Diagram of the slop refinement criterion.
:align: center

Diagram of the slope refinement criterion.
```

At each point in the domain, a ratio can be calculated for each component as:

$$
  \beta_j =  \frac{|{x_{n,j+1} - x_{n,j}} |}{\displaystyle\max_j x_n - \displaystyle\min_j x_n}
$$

Points where $ \beta_j $ is greater than $ \beta_{max} $ are marked to have an additional
grid point inserted to the right of them.

The rule-of-thumb for setting $ \beta_{max} $ is that it prevents any two adjacent points
from having a change in a solution component that is greater than $ \beta_{max} \cdot 100\% $ of
the maximum change in that component across the domain.



### Curve

This refinement option controls the maximum curvature of the solution component between
two adjacent grid points compared to the maximum curvature of that component across the
entire domain. If the ratio of the curvature between two adjacent grid points is greater
than the user-specified value, then an additional grid point is inserted between the two
points. This criterion is used to ensure that the solution is well-resolved in regions
where the solution has steep gradients.

The specific criteria that is being enforced is:

$$
  |\left( \frac{d x_n}{d z} \right)_{j+1} - \left(\frac{d x_n}{d z} \right)_j| <= \gamma_{max} (\displaystyle\max_j \frac{dx_n}{dz} - \displaystyle\min_j \frac{dx_n}{dz})
$$

Where $ \gamma_{max} $ is a number between 0 and 1 that is specified by the user, $ \frac{dx_n}{dz} $
is the derivative of the solution component being considered, and $ \displaystyle\max_j \frac{dx_n}{dz} $ and $ \displaystyle\min_j \frac{dx_n}{dz} $ are the
maximum and minimum values of the solution component derivative across the domain. Visually, these
quantities are shown in the image below.

```{figure} /_static/images/onedim/refinement_curve_metric_plot.svg
:width: 75%
:alt: Diagram of the curve refinement criterion.
:align: center

Diagram of the curve refinement criterion.
```

At each point in the domain, a ratio can be calculated for each component as:

$$
  \gamma_j =  \frac{|\left( \frac{d x_n}{d z} \right)_{j+1} - \left(\frac{d x_n}{d z} \right)_j|}{\displaystyle\max_j \frac{dx_n}{dz} - \displaystyle\min_j \frac{dx_n}{dz}}
$$

Points where $ \gamma_j $ is greater than $ \gamma_{max} $ are marked to have an additional
grid point inserted to the right of them.

The rule-of-thumb for setting $ \gamma_{max} $ is that it prevents any two adjacent points
from having a change in the slope of a solution component that is greater than $ \gamma_{max} \cdot 100\% $ of
the maximum change in that component's slope across the domain.

Note: The slope is computed using a backward difference scheme and so the estimate of the
derivative at a grid point depends on the solution at the grid point and the solution at the
grid point to the right.

$$
  \left( \frac{d x_n}{d z} \right)_j \approx \frac{x_{n,j+1} - x_{n,j}}{z_{j+1} - z_{j}}
$$

If a point is marked as needing an additional grid point, the refinement algorithm will
actually insert a point to the right of the point under consideration as well as to the
point that is to the right. This essentially splits the two intervals that are to the right
of the point under consideration.


### Prune

This refinement option works in contrast to the other refinement options. Instead of
adding grid points, this option removes grid points that are not needed. This option
lets a user specify a value of the ratio quantities described above that if the ratio
drops below this, then the point is removed. The user specifies a pruning value, $ \delta $,
and that serves as a lower-bound on the ratio quantities.

$$
  \delta = \beta_{min} = \gamma_{min}
$$

If $\beta_j < \delta$, or $\gamma_j < \delta$, then the point is removed.
