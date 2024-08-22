// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REFINE_H
#define CT_REFINE_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

class Domain1D;

//! Refine Domain1D grids so that profiles satisfy adaptation tolerances
//! @ingroup onedUtilsGroup
 class Refiner
{
public:
    Refiner(Domain1D& domain);
    virtual ~Refiner() {}
    Refiner(const Refiner&) = delete;
    Refiner& operator=(const Refiner&) = delete;

    /**
     *  Set grid refinement criteria.
     *
     * The ratio parameter is the maximum allowed ratio between grid spacing at
     * adjacent intervals. The ratio parameter considers the situation where the
     * left interval is much larger than the right interval, or if the right
     * interval is much larger than the left interval.
     *
     * The slope parameter is the maximum fractional change in the value of each
     * solution component between adjacent grid points. This quantity is used to
     * check the criteria given below:
     *    val[j+1] - val[j]  < slope*(val_max - val_min)
     * Where val_max and val_min are the maximum and minimum values of the
     * component over the entire domain.
     *
     * The curve parameter is the maximum fractional change in the derivative of
     * each solution component between adjacent grid points. This quantity is
     * used to check the criteria given below:
     *   slope[j+1] - slope[j]  < curve*(slope_max - slope_min)
     * Where slope_max and slope_min are the maximum and minimum slopes of the
     * component over the entire domain.
     *
     * The prune parameter is a threshold for removing unnecessary grid points. Both
     * the slope and curve parameters utilize a normalized ratio when making refinement
     * decisions.
     *
     * @f[
     *   \text{ratio} = \frac{|x[j+1] - x[j]|}{\text{m_slope}*(x_{\text{max}} -
     *                  x_{\text{min}})}
     * @f]
     *
     * @f[
     *  \text{ratio\_curve} = \frac{|slope[j+1] - slope[j]|}{\text{m_curve}*(slope_{\text{max}} -
     *                 slope_{\text{min}})}
     * @f]
     *
     * If this ratio exceeds 1, refinement is needed. This normalized ratio is what
     * the prune parameter is compared against. If the normalized ratio is greater
     * than the prune parameter, the grid point is kept, otherwise it is marked for
     * removal. Too much refinement can lead to unnecessary grid points, which will
     * slow down the solution process.
     *
     * A rule of thumb with the prune parameter is that, if you have specified your
     * slope and curve parameters to a desireable level, then the ratio above would
     * be 1.0 for the case when the grid points just meet the refinement criteria that
     * would be "good enough" for the solution at hand. A ratio less than than one
     * still meets the criteria, but introduces more grid points. The prune parameter
     * allows you set the point at which the refiner decides that is has overly-refined
     * the grid and can remove points that are not needed.
     *
     *  @param ratio Maximum ratio between grid spacing at adjacent intervals.
     *      That is, `(x[j+1] - x[j]) / (x[j] - x[j-1]) < ratio`
     *  @param slope Maximum fractional change in the value of each solution
     *      component between adjacent grid points.
     *  @param curve Maximum fractional change in the derivative of each
     *      solution component between adjacent grid points.
     *  @param prune Threshold for removing unnecessary grid points. `prune`
     *      should be smaller than both `slope` and `curve`. Set `prune <= 0`
     *      to disable pruning.
     */
    void setCriteria(double ratio = 10.0,
                     double slope = 0.8,
                     double curve = 0.8,
                     double prune = -0.1);

    //! Get the grid refinement criteria. @see Refiner::setCriteria
    vector<double> getCriteria()
    {
        return {m_ratio, m_slope, m_curve, m_prune};
    }

    /**
     * Set the active state for a component.
     *
     * @param comp  Component index
     * @param state  True if the component should be considered for grid refinement
     */
    void setActive(int comp, bool state = true) {
        m_active[comp] = state;
    }

    //! Set the maximum number of points allowed in the domain
    void setMaxPoints(int npmax) {
        m_npmax = npmax;
    }

    //! Returns the maximum number of points allowed in the domain
    size_t maxPoints() const {
        return m_npmax;
    }

    //! Set the minimum allowable spacing between adjacent grid points [m].
    void setGridMin(double gridmin) {
        m_gridmin = gridmin;
    }

    //! Returns the the minimum allowable spacing between adjacent
    //! grid points [m].
    double gridMin() const {
        return m_gridmin;
    }

    /**
     * Determine locations in the grid that need additional grid points and
     * updates the internal state of the Refiner with this information.
     *
     * @param[in] n  Number of grid points
     * @param[in] z  Point to array of grid points
     * @param[in] x  Pointer to solution vector
     *
     * @returns The number of new grid points needed (size of #m_insertion_points)
     */
    int analyze(size_t n, const double* z, const double* x);

    /**
     * Constructs a new grid based on refinement locations determined by the analyze
     * method.
     *
     * This function generates a new grid by inserting additional points into the
     * current grid at locations where the `analyze` method has identified a need for
     * refinement. The new grid points are placed midway between existing points deemed
     * necessary for increased resolution. If no refinement is needed, the original
     * grid is copied directly.
     *
     * @param[in] n The number of points in the original grid array `z`.
     * @param[in] z Pointer to the array of original grid points.
     * @param[in] nn The maximum number of points that the new grid array `znew` can hold.
     * @param[out] znew Pointer to the array where the new grid points will be stored.
     *
     * @return The function returns 0 upon successful creation of the new grid. Throws an exception if the provided
     *         output array size is insufficient to hold the new grid.
     */
    int getNewGrid(int n, const double* z, int nn, double* znew);

    /**
     * Returns the number of new grid points that were needed.
     */
    int numNewPoints() {
        return static_cast<int>(m_insertion_points.size());
    }

    /**
     * Displays the results of the grid refinement analysis.
     *
     * This method logs information about where new grid points have been inserted
     * and the components that required these insertions, if any. If no new points
     * were needed, it logs that the current grid is sufficient.
     *
     * @note This method should be called after analyze() to report the outcomes of
     * the refinement analysis.
     */
    void show();

    /**
     * Returns true if a new grid point is needed to the right of grid index j.
     *
     * @param j  Grid point index
     */
    bool newPointNeeded(size_t j) {
        return m_insertion_points.find(j) != m_insertion_points.end();
    }

    /**
     * Returns true if the grid point at index j should be kept.
     *
     * @param j Grid point index
     */
    bool keepPoint(size_t j) {
        return (m_keep[j] != -1);
    }

    /**
     * Returns the value of the solution component, n, at grid point j.
     *
     * @param x Solution vector
     * @param n Solution Component index
     * @param j Grid point index
     */
    double value(const double* x, size_t n, size_t j);

    //! Returns the maximum ratio of grid spacing at adjacent intervals
    double maxRatio() {
        return m_ratio;
    }

    //! Returns the maximum allowable difference in value between adjacent points
    double maxDelta() {
        return m_slope;
    }

    //! Returns the maximum allowable difference in the derivative between adjacent points
    double maxSlope() {
        return m_curve;
    }

    //! Returns the threshold for removing unnecessary grid points
    double prune() {
        return m_prune;
    }

protected:
    //! Indices of grid points that need new grid points added after them
    set<size_t> m_insertion_points;

    //! Map of grid point indices that should be kept, 1=keep, -1=remove, 0=unset
    map<size_t, int> m_keep;

    //! Names of components that require the addition of new grid points
    set<string> m_component_name;

    //! Vector of flags indicating whether each component should be considered
    //! for grid refinement
    vector<bool> m_active;

    //! Refinement criteria
    double m_ratio = 10.0;
    double m_slope = 0.8;
    double m_curve = 0.8;
    double m_prune = -0.001;

    //! Threshold parameter that is used distinguish between small fluctuations
    //! around a constant value during the refinement process.
    double m_min_range = 0.01;

    Domain1D* m_domain; //! Pointer to the domain to be refined
    size_t m_nv; //!< Number of components in the domain
    size_t m_npmax = 1000; //!< Maximum number of grid points

    double m_thresh = std::sqrt(std::numeric_limits<double>::epsilon());
    double m_gridmin = 1e-10; //!< minimum grid spacing [m]
};

}

#endif
