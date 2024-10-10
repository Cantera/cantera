# Grid Refinement for One-dimensional Flows

## Overview

The Cantera 1D solver uses a progressive mesh refinement strategy when solving the
governing equations. The system is initially solved on a coarse grid, and then the
solution is analyzed to identify regions where it is not well resolved, such as areas
with steep gradients. Additional grid points are inserted in these regions and then
the system is solved again. This process is repeated until the solution is deemed to be
well resolved. This incremental refinement procedure maximizes the probability that the
initial guess for each subsequent refined solution lies within the domain of
convergence of the Newton solver.

The algorithm tracks information about which grid points
should be kept and which ones need to have additional grid points inserted to the
right of them. This data is used to refine the grid as well as to keep too many
points from being removed during a refinement step.


## Refinement Criteria

There are 4 refinement criteria that can be set by a user to control the grid
refinement process: ratio, slope, curve, and prune.

### Ratio

This refinement option controls the maximum ratio of the spacing between the grid
points that are adjacent to a point under consideration. If the ratio between two
adjacent grid points is greater than the ratio
specified by the user ($ \alpha_{max} $), then an additional grid point is inserted
between the two points. Depending on which interval adjacent to a point under
consideration is larger, two conditions are checked:

#### Larger Left Interval

For grid points that have a larger spacing on the left interval, as shown in the figure below,

```{figure} /_static/images/onedim/refinement_large_left_interval.svg
:width: 75%
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
  \beta_j = \frac{\left|{x_{n,j+1} - x_{n,j}}\right|}{\displaystyle\max_j x_n - \displaystyle\min_j x_n} \le \beta_{max}
$$

Where $ \beta_j $ is computed at each domain point, $ \beta_{max} $ is a number
between 0 and 1 that is specified by the user, $ x_n $ is the solution component
being considered, and $ \displaystyle\max_j x_n $ and $ \displaystyle\min_j x_n $ are
the maximum and minimum values of the solution component across the domain. Points
where $ \beta_j $ is greater than $ \beta_{max} $ are marked to have an additional grid
point inserted to the right of them. Visually, the minimum and maximum solution
component values are shown in the image below. For a point, j, considered where
the change between the solution at j and the solution at j+1 is larger than the
global maximum change, the point under consideration is marked as needing an
additional grid point inserted to the right of it, and the point to the right of
the point under consideration is marked to be kept.

```{figure} /_static/images/onedim/refinement_slope_metric_plot.svg
:width: 75%
:align: center

Diagram of the slope refinement criterion.
```


### Curve

This refinement option controls the maximum curvature of the solution component between
two adjacent grid points compared to the maximum curvature of that component across the
entire domain. If the ratio of the curvature between two adjacent grid points is greater
than the user-specified value, then an additional grid point is inserted between the two
points. This criterion is used to ensure that the solution is well-resolved in regions
where the solution has steep gradients.

The specific criteria that is being enforced is:

$$
  \gamma_j = \frac{\left| x'_{n, j+1} - x'_{n, j} \right|}
             {\displaystyle\max_j x'_n - \displaystyle\min_j x'_n} \le \gamma_{max}
$$

Where $ \gamma_j$ is computed at each grid point, $ \gamma_{max} $ is a number between
0 and 1 that is specified by the user, $ x'_n $ is the derivative of the
solution component being considered, and $ \displaystyle\max_j x'_n $ and
$ \displaystyle\min_j x'_n $ are the maximum and minimum values of the
solution component derivative across the domain. Points where $ \gamma_j $ is greater
than $ \gamma_{max} $ are marked to have an additional grid point inserted to the
right of them. Visually, the minimum and maximum global quantities are shown in the
image below.

```{figure} /_static/images/onedim/refinement_curve_metric_plot.svg
:width: 75%
:alt: Diagram of the curve refinement criterion.
:align: center

Diagram of the curve refinement criterion.
```

**Note:** The slope is computed using a backward difference scheme and so the estimate of
the derivative at a grid point depends on the solution at the grid point and the
solution at the grid point to the right.

$$
  x'_{n, j} \approx \frac{x_{n,j+1} - x_{n,j}}{z_{j+1} - z_{j}}
$$

If a point is marked as needing an additional grid point, the refinement algorithm will
actually insert a point to the right of the point under consideration as well as to the
point that is to the right. This essentially splits the two intervals that are to the
right of the point under consideration.


### Prune

This refinement option works in contrast to the other refinement options. Instead of
adding grid points, this option removes grid points that are not needed. The pruning
parameter helps the grid refinement algorithm identify areas where the grid may be
overly refined. You specify a pruning value $ \delta $ that serves as a lower-bound on
the ratio quantities. This helps prevent unnecessary grid refinement by
removing points that are not needed, optimizing computational efficiency without
sacrificing the accuracy of the solution.

The following criteria are used to decide if a point is marked for removal.

$$
  \frac{\beta_j}{\beta_{\text{max}}} < \delta \quad \text{or}
  \quad \frac{\gamma_j}{\gamma_{\text{max}}} < \delta
$$

If $ \frac{\beta_j}{\beta_{max}} $, or $ \frac{\gamma_j}{\gamma_{max}} $ are less
than $ \delta $, then the point is marked for removal.

For example, consider a case where the slope refinement criteria is set to 0.5, we
can look at the pruning criteria as $$ \beta_j < \beta_{max} \delta $$. In this form
it is easy to see that the pruning parameter defines "What fraction of the original
$\beta_{max}$ is acceptable enough to remove a point?" Continuing with this example,
consider a case where the refinement has progressed such that these are the values of
the $$ \beta_j $$ at each point in a domain with 10 grid points.

$$
  \beta_j = [0.11, 0.16, 0.05, 0.1, 0.22, 0.09, 0.05, 0.02, 0.37, 0.29]
$$

All of the values are below the original $$ \beta_{max} $$ of 0.5, but some are much
lower than that threshold. If the pruning parameter is set to a value of 0.2, then
the effective cutoff threshold for removing a point would be 0.5*0.2 = 0.1. With this
pruning criteria, the pruned data would be:

$$
  \beta_j = [0.11, 0.16, 0.1, 0.22, 0.37, 0.29]
$$