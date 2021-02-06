//! @file refine.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/refine.h"
#include "cantera/oneD/StFlow.h"

using namespace std;

namespace Cantera
{
Refiner::Refiner(Domain1D& domain) :
    m_ratio(10.0), m_slope(0.8), m_curve(0.8), m_prune(-0.001),
    m_min_range(0.01), m_domain(&domain), m_npmax(1000),
    m_gridmin(1e-10)
{
    m_nv = m_domain->nComponents();
    m_active.resize(m_nv, true);
    m_thresh = std::sqrt(std::numeric_limits<double>::epsilon());
}

void Refiner::setCriteria(doublereal ratio, doublereal slope,
                          doublereal curve, doublereal prune)
{
    if (ratio < 2.0) {
        throw CanteraError("Refiner::setCriteria",
            "'ratio' must be greater than 2.0 ({} was specified).", ratio);
    } else if (slope < 0.0 || slope > 1.0) {
        throw CanteraError("Refiner::setCriteria",
            "'slope' must be between 0.0 and 1.0 ({} was specified).", slope);
    } else if (curve < 0.0 || curve > 1.0) {
        throw CanteraError("Refiner::setCriteria",
            "'curve' must be between 0.0 and 1.0 ({} was specified).", curve);
    } else if (prune > curve || prune > slope) {
        throw CanteraError("Refiner::setCriteria",
            "'prune' must be less than 'curve' and 'slope' ({} was specified).",
            prune);
    }
    m_ratio = ratio;
    m_slope = slope;
    m_curve = curve;
    m_prune = prune;
}

int Refiner::analyze(size_t n, const doublereal* z,
                     const doublereal* x)
{
    if (n >= m_npmax) {
        throw CanteraError("Refiner::analyze", "max number of grid points reached ({}).", m_npmax);
    }

    if (m_domain->nPoints() <= 1) {
        return 0;
    }

    m_loc.clear();
    m_c.clear();
    m_keep.clear();

    m_keep[0] = 1;
    m_keep[n-1] = 1;

    m_nv = m_domain->nComponents();

    // check consistency
    if (n != m_domain->nPoints()) {
        throw CanteraError("Refiner::analyze", "inconsistent");
    }

    // find locations where cell size ratio is too large.
    vector_fp v(n), s(n-1);

    vector_fp dz(n-1);
    for (size_t j = 0; j < n-1; j++) {
        dz[j] = z[j+1] - z[j];
    }

    for (size_t i = 0; i < m_nv; i++) {
        if (m_active[i]) {
            string name = m_domain->componentName(i);
            // get component i at all points
            for (size_t j = 0; j < n; j++) {
                v[j] = value(x, i, j);
            }

            // slope of component i
            for (size_t j = 0; j < n-1; j++) {
                s[j] = (value(x, i, j+1) - value(x, i, j))/(z[j+1] - z[j]);
            }

            // find the range of values and slopes
            doublereal vmin = *min_element(v.begin(), v.end());
            doublereal vmax = *max_element(v.begin(), v.end());
            doublereal smin = *min_element(s.begin(), s.end());
            doublereal smax = *max_element(s.begin(), s.end());

            // max absolute values of v and s
            doublereal aa = std::max(fabs(vmax), fabs(vmin));
            doublereal ss = std::max(fabs(smax), fabs(smin));

            // refine based on component i only if the range of v is
            // greater than a fraction 'min_range' of max |v|. This
            // eliminates components that consist of small fluctuations
            // on a constant background.
            if ((vmax - vmin) > m_min_range*aa) {
                // maximum allowable difference in value between adjacent
                // points.
                doublereal dmax = m_slope*(vmax - vmin) + m_thresh;
                for (size_t j = 0; j < n-1; j++) {
                    doublereal r = fabs(v[j+1] - v[j])/dmax;
                    if (r > 1.0 && dz[j] >= 2 * m_gridmin) {
                        m_loc[j] = 1;
                        m_c[name] = 1;
                    }
                    if (r >= m_prune) {
                        m_keep[j] = 1;
                        m_keep[j+1] = 1;
                    } else if (m_keep[j] == 0) {
                        m_keep[j] = -1;
                    }
                }
            }

            // refine based on the slope of component i only if the
            // range of s is greater than a fraction 'min_range' of max
            // |s|. This eliminates components that consist of small
            // fluctuations on a constant slope background.
            if ((smax - smin) > m_min_range*ss) {
                // maximum allowable difference in slope between
                // adjacent points.
                doublereal dmax = m_curve*(smax - smin);
                for (size_t j = 0; j < n-2; j++) {
                    doublereal r = fabs(s[j+1] - s[j]) / (dmax + m_thresh/dz[j]);
                    if (r > 1.0 && dz[j] >= 2 * m_gridmin &&
                            dz[j+1] >= 2 * m_gridmin) {
                        m_c[name] = 1;
                        m_loc[j] = 1;
                        m_loc[j+1] = 1;
                    }
                    if (r >= m_prune) {
                        m_keep[j+1] = 1;
                    } else if (m_keep[j+1] == 0) {
                        m_keep[j+1] = -1;
                    }
                }
            }
        }
    }

    StFlow* fflame = dynamic_cast<StFlow*>(m_domain);

    // Refine based on properties of the grid itself
    for (size_t j = 1; j < n-1; j++) {
        // Add a new point if the ratio with left interval is too large
        if (dz[j] > m_ratio*dz[j-1]) {
            m_loc[j] = 1;
            m_c[fmt::format("point {}", j)] = 1;
            m_keep[j-1] = 1;
            m_keep[j] = 1;
            m_keep[j+1] = 1;
            m_keep[j+2] = 1;
        }

        // Add a point if the ratio with right interval is too large
        if (dz[j] < dz[j-1]/m_ratio) {
            m_loc[j-1] = 1;
            m_c[fmt::format("point {}", j-1)] = 1;
            m_keep[j-2] = 1;
            m_keep[j-1] = 1;
            m_keep[j] = 1;
            m_keep[j+1] = 1;
        }

        // Keep the point if removing would make the ratio with the left
        // interval too large.
        if (j > 1 && z[j+1]-z[j-1] > m_ratio * dz[j-2]) {
            m_keep[j] = 1;
        }

        // Keep the point if removing would make the ratio with the right
        // interval too large.
        if (j < n-2 && z[j+1]-z[j-1] > m_ratio * dz[j+1]) {
            m_keep[j] = 1;
        }

        // Keep the point where the temperature is fixed
        if (fflame && fflame->domainType() == cFreeFlow && z[j] == fflame->m_zfixed) {
            m_keep[j] = 1;
        }
    }

    // Don't allow pruning to remove multiple adjacent grid points
    // in a single pass.
    for (size_t j = 2; j < n-1; j++) {
        if (m_keep[j] == -1 && m_keep[j-1] == -1) {
            m_keep[j] = 1;
        }
    }

    return int(m_loc.size());
}

double Refiner::value(const double* x, size_t i, size_t j)
{
    return x[m_domain->index(i,j)];
}

void Refiner::show()
{
    if (!m_loc.empty()) {
        writeline('#', 78);
        writelog(string("Refining grid in ") +
                 m_domain->id()+".\n"
                 +"    New points inserted after grid points ");
        for (const auto& loc : m_loc) {
            writelog("{} ", loc.first);
        }
        writelog("\n");
        writelog("    to resolve ");
        for (const auto& c : m_c) {
            writelog(string(c.first)+" ");
        }
        writelog("\n");
        writeline('#', 78);
    } else if (m_domain->nPoints() > 1) {
        writelog("no new points needed in "+m_domain->id()+"\n");
    }
}

int Refiner::getNewGrid(int n, const doublereal* z,
                        int nn, doublereal* zn)
{
    int nnew = static_cast<int>(m_loc.size());
    if (nnew + n > nn) {
        throw CanteraError("Refine::getNewGrid", "array size too small.");
    }

    if (m_loc.empty()) {
        copy(z, z + n, zn);
        return 0;
    }

    int jn = 0;
    for (int j = 0; j < n - 1; j++) {
        zn[jn] = z[j];
        jn++;
        if (m_loc.count(j)) {
            zn[jn] = 0.5*(z[j] + z[j+1]);
            jn++;
        }
    }
    zn[jn] = z[n-1];
    return 0;
}
}
