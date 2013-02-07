#ifndef CT_REFINE_H
#define CT_REFINE_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

class Domain1D;

class Refiner
{
public:
    Refiner(Domain1D& domain);
    virtual ~Refiner() {}

    void setCriteria(doublereal ratio = 10.0,
                     doublereal slope = 0.8,
                     doublereal curve = 0.8,
                     doublereal prune = -0.1) {
        m_ratio = ratio;
        m_slope = slope;
        m_curve = curve;
        m_prune = prune;
    }
    void setActive(int comp, bool state = true) {
        m_active[comp] = state;
    }
    void setMaxPoints(int npmax) {
        m_npmax = npmax;
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

    int analyze(size_t n, const doublereal* z, const doublereal* x);
    int getNewGrid(int n, const doublereal* z, int nn, doublereal* znew);
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
    std::map<size_t, int> m_loc;
    std::map<size_t, int> m_keep;
    std::map<std::string, int> m_c;
    std::vector<bool> m_active;
    doublereal m_ratio, m_slope, m_curve, m_prune;
    doublereal m_min_range;
    Domain1D* m_domain;
    size_t m_nv, m_npmax;
    doublereal m_thresh;
    doublereal m_gridmin; //!< minimum grid spacing [m]
};

}

#endif
