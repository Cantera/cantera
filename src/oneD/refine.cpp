//! @file refine.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/refine.h"
#include "cantera/oneD/Flow1D.h"
#include "cantera/base/global.h"

using namespace std;

namespace Cantera
{
Refiner::Refiner(Domain1D& domain) :
    m_domain(&domain)
{
    m_nv = m_domain->nComponents();
    m_active.resize(m_nv, true);
}

void Refiner::setCriteria(double ratio, double slope, double curve, double prune)
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

int Refiner::analyze(size_t n, const double* z, const double* x)
{
    if (n >= m_npmax) {
        throw CanteraError("Refiner::analyze", "max number of grid points reached ({}).", m_npmax);
    }

    if (m_domain->nPoints() <= 1) {
        return 0;
    }

    // check consistency
    if (n != m_domain->nPoints()) {
        throw CanteraError("Refiner::analyze", "number of grid points provided does not match domain size.");
    }

    // Reset the state of the refiner
    m_insertPts.clear();
    m_componentNames.clear();
    m_keep.clear();

    // Keep the first and last grid points
    m_keep[0] = KEEP;
    m_keep[n-1] = KEEP;

    m_nv = m_domain->nComponents();

    // find locations where cell size ratio is too large.
    vector<double> val(n);
    vector<double> slope(n-1);

    vector<double> dz(n-1); // Store the right-looking grid spacings
    for (size_t j = 0; j < n-1; j++) {
        dz[j] = z[j+1] - z[j];
    }

    for (size_t i = 0; i < m_nv; i++) {
        if (m_active[i]) {
            string name = m_domain->componentName(i);
            // get component i at all points
            for (size_t j = 0; j < n; j++) {
                val[j] = value(x, i, j);
            }

            // slope of component i (using forward difference)
            for (size_t j = 0; j < n-1; j++) {
                slope[j] = (val[j+1] - val[j]) / dz[j];
            }

            // find the range of values and slopes of component i over the domain
            double valMin = *min_element(val.begin(), val.end());
            double valMax = *max_element(val.begin(), val.end());
            double slopeMin = *min_element(slope.begin(), slope.end());
            double slopeMax = *max_element(slope.begin(), slope.end());

            // max absolute values of val and slope
            double valMagnitude = std::max(fabs(valMax), fabs(valMin));
            double slopeMagnitude = std::max(fabs(slopeMax), fabs(slopeMin));

            // refine based on component i only if the range of val is greater than a
            // fraction 'min_range' of max |val|. This eliminates components that
            // consist of small fluctuations around a constant value.
            if (valMax - valMin > m_min_range*valMagnitude) {
                // maximum allowable difference in value between adjacent points. Based
                // on the global min and max values of the component over the domain.
                double max_change = m_slope*(valMax - valMin);
                for (size_t j = 0; j < n-1; j++) {
                    double ratio = fabs(val[j+1] - val[j]) / (max_change + m_thresh);
                    if (ratio > 1.0 && dz[j] >= 2 * m_gridmin) {
                        m_insertPts.insert(j);
                        m_componentNames.insert(name);
                    }
                    if (ratio >= m_prune) {
                        m_keep[j] = KEEP;
                        m_keep[j+1] = KEEP;
                    } else if (m_keep[j] == UNSET) {
                        m_keep[j] = REMOVE;
                    }
                }
            }

            // refine based on the slope of component i only if the range of s is
            // greater than a fraction 'min_range' of max|s|. This eliminates
            // components that consist of small fluctuations on a constant slope
            // background.
            if (slopeMax - slopeMin > m_min_range*slopeMagnitude) {
                // maximum allowable difference in slope between adjacent points.
                double max_change = m_curve*(slopeMax - slopeMin);
                for (size_t j = 0; j < n-2; j++) {
                    // Using the solution component absolute tolerance (m_thresh),
                    // an absolute tolerance for the change in slope can be estimated
                    // for an interval dz as m_thresh/dz.
                    double ratio = fabs(slope[j+1] - slope[j]) / (max_change + m_thresh/dz[j]);
                    if (ratio > 1.0 && dz[j] >= 2*m_gridmin && dz[j+1] >= 2*m_gridmin) {
                        m_componentNames.insert(name);
                        m_insertPts.insert(j);
                        m_insertPts.insert(j+1);
                    }
                    if (ratio >= m_prune) {
                        m_keep[j+1] = KEEP;
                    } else if (m_keep[j+1] == UNSET) {
                        m_keep[j+1] = REMOVE;
                    }
                }
            }
        }
    }

    // Refine based on properties of the grid itself
    for (size_t j = 1; j < n-1; j++) {
        // Add a new point if the ratio with left interval is too large.
        // Extra points around the interval set under consideration are kept.
        if (dz[j] > m_ratio*dz[j-1]) {
            m_insertPts.insert(j);
            m_componentNames.insert(fmt::format("point {}", j));
            m_keep[j-1] = KEEP;
            m_keep[j] = KEEP;
            m_keep[j+1] = KEEP;
            m_keep[j+2] = KEEP;
        }

        // Add a point if the ratio with right interval is too large
        if (dz[j-1] > m_ratio*dz[j]) {
            m_insertPts.insert(j-1);
            m_componentNames.insert(fmt::format("point {}", j-1));
            m_keep[j-2] = KEEP;
            m_keep[j-1] = KEEP;
            m_keep[j] = KEEP;
            m_keep[j+1] = KEEP;
        }

        // Keep the point if removing would make the ratio with the left interval too
        // large.
        if (j > 1 && z[j+1]-z[j-1] > m_ratio * dz[j-2]) {
            m_keep[j] = KEEP;
        }

        // Keep the point if removing would make the ratio with the right interval too
        // large.
        if (j < n-2 && z[j+1]-z[j-1] > m_ratio * dz[j+1]) {
            m_keep[j] = KEEP;
        }

        Flow1D* fflame = dynamic_cast<Flow1D*>(m_domain);
        // Keep the point where the temperature is fixed
        if (fflame && fflame->isFree() && z[j] == fflame->m_zfixed) {
            m_keep[j] = KEEP;
        }

        // Keep the point if it is a control point used for two-point flame control
        if (fflame && fflame->twoPointControlEnabled() &&
            (z[j] == fflame->leftControlPointCoordinate() ||
             z[j] == fflame->rightControlPointCoordinate()))
        {
            m_keep[j] = KEEP;
        }
    }

    // Don't allow pruning to remove multiple adjacent grid points
    // in a single pass.
    for (size_t j = 2; j < n-1; j++) {
        if (m_keep[j] == REMOVE && m_keep[j-1] == REMOVE) {
            m_keep[j] = KEEP;
        }
    }

    return int(m_insertPts.size());
}

double Refiner::value(const double* x, size_t n, size_t j)
{
    return x[m_domain->index(n,j)];
}

void Refiner::show()
{
    if (!m_insertPts.empty()) {
        writeline('#', 78);
        writelog(string("Refining grid in ") +
                 m_domain->id()+".\n"
                 +"    New points inserted after grid points ");
        for (const auto& loc : m_insertPts) {
            writelog("{} ", loc);
        }
        writelog("\n");
        writelog("    to resolve ");
        for (const auto& c : m_componentNames) {
            writelog(c + " ");
        }
        writelog("\n");
        writeline('#', 78);
    } else if (m_domain->nPoints() > 1) {
        writelog("no new points needed in "+m_domain->id()+"\n");
    }
}

int Refiner::getNewGrid(int n, const double* z, int nn, double* zn)
{
    warn_deprecated(
        "Refiner::getNewGrid",
        "Deprecated in Cantera 3.1; unused function that will be removed.");

    int nnew = static_cast<int>(m_insertPts.size());
    if (nnew + n > nn) {
        throw CanteraError("Refine::getNewGrid", "array size too small.");
    }

    if (m_insertPts.empty()) {
        copy(z, z + n, zn);
        return 0;
    }

    int jn = 0;
    for (int j = 0; j < n - 1; j++) {
        zn[jn] = z[j];
        jn++;
        if (m_insertPts.count(j)) {
            zn[jn] = 0.5*(z[j] + z[j+1]);
            jn++;
        }
    }
    zn[jn] = z[n-1];
    return 0;
}
}
