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
    m_insertion_points.clear();
    m_component_name.clear();
    m_keep.clear();

    // Keep the first and last grid points
    m_keep[0] = 1;
    m_keep[n-1] = 1;

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

            // slope of component i
            for (size_t j = 0; j < n-1; j++) {
                slope[j] = (val[j+1] - val[j]) / dz[j];
            }

            // find the range of values and slopes of component i over the domain
            double val_min = *min_element(val.begin(), val.end());
            double val_max = *max_element(val.begin(), val.end());
            double slope_min = *min_element(slope.begin(), slope.end());
            double slope_max = *max_element(slope.begin(), slope.end());

            // max absolute values of val and slope
            double val_magnitude = std::max(fabs(val_max), fabs(val_min));
            double slope_magnitude = std::max(fabs(slope_max), fabs(slope_min));

            // refine based on component i only if the range of val is greater than a
            // fraction 'min_range' of max |val|. This eliminates components that
            // consist of small fluctuations around a constant value.
            if ((val_max - val_min) > m_min_range*val_magnitude) {
                // maximum allowable difference in value between adjacent points. Based
                // on the global min and max values of the component over the domain.
                double max_change = m_slope*(val_max - val_min) + m_thresh;
                for (size_t j = 0; j < n-1; j++) {
                    double ratio = fabs(val[j+1] - val[j]) / max_change;
                    if (ratio > 1.0 && dz[j] >= 2 * m_gridmin) {
                        m_insertion_points.insert(j);
                        m_component_name.insert(name);
                    }
                    if (ratio >= m_prune) {
                        m_keep[j] = 1;
                        m_keep[j+1] = 1;
                    } else if (m_keep[j] == 0) {
                        m_keep[j] = -1;
                    }
                }
            }

            // refine based on the slope of component i only if the range of s is
            // greater than a fraction 'min_range' of max|s|. This eliminates
            // components that consist of small fluctuations on a constant slope
            // background.
            if ((slope_max - slope_min) > m_min_range*slope_magnitude) {
                // maximum allowable difference in slope between adjacent points.
                double max_change = m_curve*(slope_max - slope_min) + m_thresh;
                for (size_t j = 0; j < n-2; j++) {
                    double ratio = fabs(slope[j+1] - slope[j]) / max_change;
                    if (ratio > 1.0 && dz[j] >= 2*m_gridmin && dz[j+1] >= 2*m_gridmin) {
                        m_component_name.insert(name);
                        m_insertion_points.insert(j);
                        m_insertion_points.insert(j+1);
                    }
                    if (ratio >= m_prune) {
                        m_keep[j+1] = 1;
                    } else if (m_keep[j+1] == 0) {
                        m_keep[j+1] = -1;
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
            m_insertion_points.insert(j);
            m_component_name.insert(fmt::format("point {}", j));
            m_keep[j-1] = 1;
            m_keep[j] = 1;
            m_keep[j+1] = 1;
            m_keep[j+2] = 1;
        }

        // Add a point if the ratio with right interval is too large
        if (dz[j-1] > m_ratio*dz[j]) {
            m_insertion_points.insert(j-1);
            m_component_name.insert(fmt::format("point {}", j-1));
            m_keep[j-2] = 1;
            m_keep[j-1] = 1;
            m_keep[j] = 1;
            m_keep[j+1] = 1;
        }

        // Keep the point if removing would make the ratio with the left interval too
        // large.
        if (j > 1 && z[j+1]-z[j-1] > m_ratio * dz[j-2]) {
            m_keep[j] = 1;
        }

        // Keep the point if removing would make the ratio with the right interval too
        // large.
        if (j < n-2 && z[j+1]-z[j-1] > m_ratio * dz[j+1]) {
            m_keep[j] = 1;
        }

        Flow1D* fflame = dynamic_cast<Flow1D*>(m_domain);
        // Keep the point where the temperature is fixed
        if (fflame && fflame->isFree() && z[j] == fflame->m_zfixed) {
            m_keep[j] = 1;
        }

        // Keep the point if it is a control point used for two-point flame control
        if (fflame && fflame->twoPointControlEnabled() &&
            (z[j] == fflame->leftControlPointCoordinate() ||
             z[j] == fflame->rightControlPointCoordinate()))
        {
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

    return int(m_insertion_points.size());
}

double Refiner::value(const double* x, size_t n, size_t j)
{
    return x[m_domain->index(n,j)];
}

void Refiner::show()
{
    if (!m_insertion_points.empty()) {
        writeline('#', 78);
        writelog(string("Refining grid in ") +
                 m_domain->id()+".\n"
                 +"    New points inserted after grid points ");
        for (const auto& loc : m_insertion_points) {
            writelog("{} ", loc);
        }
        writelog("\n");
        writelog("    to resolve ");
        for (const auto& c : m_component_name) {
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
    int nnew = static_cast<int>(m_insertion_points.size());
    if (nnew + n > nn) {
        throw CanteraError("Refine::getNewGrid", "array size too small.");
    }

    if (m_insertion_points.empty()) {
        copy(z, z + n, zn);
        return 0;
    }

    int jn = 0;
    for (int j = 0; j < n - 1; j++) {
        zn[jn] = z[j];
        jn++;
        if (m_insertion_points.count(j)) {
            zn[jn] = 0.5*(z[j] + z[j+1]);
            jn++;
        }
    }
    zn[jn] = z[n-1];
    return 0;
}
}
