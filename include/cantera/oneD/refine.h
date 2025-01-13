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

    //! Set grid refinement criteria
    /*!
     *  @param ratio Maximum ratio between grid spacing at adjacent intervals.
     *      That is, `(x[j+1] - x[j]) / (x[j] - x[j-1]) < ratio`
     *  @param slope Maximum fractional change in the value of each solution
     *      component between adjacent grid points
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

    int analyze(size_t n, const double* z, const double* x);
    int getNewGrid(int n, const double* z, int nn, double* znew);

    double remeshFromSolution(int np, const double* z, const double* x, const double dist_min=1e9, const double domain_size=50.0); //from MUTAGEN

    int indxtp(int np,double val,const double* array); //from MUTAGEN

    void getRatio(int np, const double* z, double* x); //from MUTAGEN

    int nNewPoints() {
        return static_cast<int>(m_loc.size());
    }
    void show();
    bool newPointNeeded(size_t j) {
        return m_loc.find(j) != m_loc.end();
    }
    bool keepPoint(size_t j) {
        return (m_keep[j] != -1);
    }


    double z_new(const int m) //from MUTAGEN
    {
      return m_z_new[m];
    }

    int z_new_size() //from MUTAGEN
    {
      return m_z_new.size();
    }

     double grad_max() //from MUTAGEN
    {
      return m_grad_max;
    }

    double curve_max() //from MUTAGEN
    {
      return m_curve_max;
    }


    double value(const double* x, size_t i, size_t j);

    double maxRatio() {
        return m_ratio;
    }
    double maxDelta() {
        return m_slope;
    }
    double maxSlope() {
        return m_curve;
    }
    double prune() {
        return m_prune;
    }

protected:
    //! Indices of grid points that need new grid points added after them
    set<size_t> m_loc;
    map<size_t, int> m_keep;
    //! Names of components that require the addition of new grid points
    set<string> m_c;
    vector<bool> m_active;
    double m_ratio = 10.0;
    double m_slope = 0.8;
    double m_curve = 0.8;
    double m_prune = -0.001;
    double m_min_range = 0.01;
    Domain1D* m_domain;
    size_t m_nv;
    size_t m_npmax = 1000;
    vector<double> m_z_new; //from MUTAGEN
    double m_grad_max, m_curve_max; //from MUTAGEN
    double m_thresh = std::sqrt(std::numeric_limits<double>::epsilon());
    double m_gridmin = 1e-10; //!< minimum grid spacing [m]
};

}

#endif
