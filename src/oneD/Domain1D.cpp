/**
 * @file Domain1D.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/Domain1D.h"
#include "cantera/oneD/MultiJac.h"
#include "cantera/oneD/refine.h"
#include "cantera/base/AnyMap.h"
#include "cantera/base/Solution.h"
#include "cantera/base/SolutionArray.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{

Domain1D::Domain1D(size_t nv, size_t points, double time)
{
    resize(nv, points);
}

Domain1D::~Domain1D()
{
    if (m_solution) {
        m_solution->thermo()->removeSpeciesLock();
    }
}

void Domain1D::setSolution(shared_ptr<Solution> sol)
{
    if (!sol || !(sol->thermo())) {
        throw CanteraError("Domain1D::setSolution",
            "Missing or incomplete Solution object.");
    }
    warn_deprecated("Domain1D::setSolution",
        "After Cantera 3.2, a change of domain contents after instantiation "
        "will be disabled.");
    if (m_solution) {
        m_solution->thermo()->removeSpeciesLock();
    }
    m_solution = sol;
    m_solution->thermo()->addSpeciesLock();
}

string Domain1D::info(const vector<string>& keys, int rows, int width)
{
    return toArray()->info(keys, rows, width);
}

void Domain1D::resize(size_t nv, size_t np)
{
    // if the number of components is being changed, then a
    // new grid refiner is required.
    if (nv != m_nv || !m_refiner) {
        m_nv = nv;
        m_refiner = make_unique<Refiner>(*this);
    }
    m_nv = nv;
    m_name.resize(m_nv,"");
    m_max.resize(m_nv, 0.0);
    m_min.resize(m_nv, 0.0);
    // Default error tolerances for all domains
    m_rtol_ss.resize(m_nv, 1.0e-4);
    m_atol_ss.resize(m_nv, 1.0e-9);
    m_rtol_ts.resize(m_nv, 1.0e-4);
    m_atol_ts.resize(m_nv, 1.0e-11);
    m_points = np;
    m_z.resize(np, 0.0);
    m_slast.resize(m_nv * m_points, 0.0);
    locate();
}

string Domain1D::componentName(size_t n) const
{
    if (n < m_nv) {
        if (m_name[n] != "") {
            return m_name[n];
        } else {
            return fmt::format("component {}", n);
        }
    }
    throw IndexError("Domain1D::componentName", "component", n, m_nv);
}

size_t Domain1D::componentIndex(const string& name, bool checkAlias) const
{
    for (size_t n = 0; n < m_nv; n++) {
        if (name == componentName(n)) {
            return n;
        }
    }
    throw CanteraError("Domain1D::componentIndex",
                       "Component '{}' not found", name);
}

bool Domain1D::hasComponent(const string& name, bool checkAlias) const
{
    for (size_t n = 0; n < m_nv; n++) {
        if (name == componentName(n)) {
            return true;
        }
    }
    throw CanteraError("Domain1D::hasComponent",
                       "Component '{}' not found", name);
}

void Domain1D::setTransientTolerances(double rtol, double atol, size_t n)
{
    if (n == npos) {
        for (n = 0; n < m_nv; n++) {
            m_rtol_ts[n] = rtol;
            m_atol_ts[n] = atol;
        }
    } else {
        m_rtol_ts[n] = rtol;
        m_atol_ts[n] = atol;
    }
}

void Domain1D::setSteadyTolerances(double rtol, double atol, size_t n)
{
    if (n == npos) {
        for (n = 0; n < m_nv; n++) {
            m_rtol_ss[n] = rtol;
            m_atol_ss[n] = atol;
        }
    } else {
        m_rtol_ss[n] = rtol;
        m_atol_ss[n] = atol;
    }
}

void Domain1D::needJacUpdate()
{
    if (m_container) {
        m_container->linearSolver()->setAge(10000);
        m_container->saveStats();
    }
}

AnyMap Domain1D::getMeta() const
{
    auto wrap_tols = [this](const vector<double>& tols) {
        // If all tolerances are the same, just store the scalar value.
        // Otherwise, store them by component name
        set<double> unique_tols(tols.begin(), tols.end());
        if (unique_tols.size() == 1) {
            return AnyValue(tols[0]);
        } else {
            AnyMap out;
            for (size_t i = 0; i < tols.size(); i++) {
                out[componentName(i)] = tols[i];
            }
            return AnyValue(out);
        }
    };
    AnyMap state;
    state["type"] = type();
    state["points"] = static_cast<long int>(nPoints());
    if (nComponents() && nPoints()) {
        state["tolerances"]["transient-abstol"] = wrap_tols(m_atol_ts);
        state["tolerances"]["steady-abstol"] = wrap_tols(m_atol_ss);
        state["tolerances"]["transient-reltol"] = wrap_tols(m_rtol_ts);
        state["tolerances"]["steady-reltol"] = wrap_tols(m_rtol_ss);
    }
    return state;
}

void Domain1D::setMeta(const AnyMap& meta)
{
    auto set_tols = [&](const AnyValue& tols, const string& which, vector<double>& out)
    {
        if (!tols.hasKey(which)) {
            return;
        }
        const auto& tol = tols[which];
        if (tol.isScalar()) {
            out.assign(nComponents(), tol.asDouble());
        } else {
            for (size_t i = 0; i < nComponents(); i++) {
                string name = componentName(i);
                if (tol.hasKey(name)) {
                    out[i] = tol[name].asDouble();
                } else {
                    warn_user("Domain1D::setMeta", "No {} found for component '{}'",
                              which, name);
                }
            }
        }
    };

    if (meta.hasKey("tolerances")) {
        const auto& tols = meta["tolerances"];
        set_tols(tols, "transient-abstol", m_atol_ts);
        set_tols(tols, "transient-reltol", m_rtol_ts);
        set_tols(tols, "steady-abstol", m_atol_ss);
        set_tols(tols, "steady-reltol", m_rtol_ss);
    }
}

void Domain1D::locate()
{
    if (m_left) {
        // there is a domain on the left, so the first grid point in this domain
        // is one more than the last one on the left
        m_jstart = m_left->lastPoint() + 1;

        // the starting location in the solution vector
        m_iloc = m_left->loc() + m_left->size();
    } else {
        // this is the left-most domain
        m_jstart = 0;
        m_iloc = 0;
    }
    // if there is a domain to the right of this one, then repeat this for it
    if (m_right) {
        m_right->locate();
    }
}

void Domain1D::setupGrid(const vector<double>& grid)
{
    setupGrid(grid.size(), grid.data());
}

void Domain1D::setupGrid(size_t n, const double* z)
{
    if (n > 1) {
        resize(m_nv, n);
        for (size_t j = 0; j < m_points; j++) {
            m_z[j] = z[j];
        }
    }
}

void Domain1D::setupUniformGrid(size_t points, double length, double start)
{
    vector<double> grid(points);
    double dz = length / (double)(points - 1);
    for (int iz = 0; iz < points; iz++) {
        grid[iz] = start + iz * dz;
    }
    setupGrid(grid);
}

void Domain1D::show(const double* x)
{
    size_t nn = m_nv/5;
    for (size_t i = 0; i < nn; i++) {
        writeline('-', 79, false, true);
        writelog("\n          z ");
        for (size_t n = 0; n < 5; n++) {
            writelog(" {:>10s} ", componentName(i*5 + n));
        }
        writeline('-', 79, false, true);
        for (size_t j = 0; j < m_points; j++) {
            writelog("\n {:10.4g} ", m_z[j]);
            for (size_t n = 0; n < 5; n++) {
                writelog(" {:10.4g} ", x[index(i*5 + n, j)]);
            }
        }
        writelog("\n");
    }
    size_t nrem = m_nv - 5*nn;
    writeline('-', 79, false, true);
    writelog("\n          z ");
    for (size_t n = 0; n < nrem; n++) {
        writelog(" {:>10s} ", componentName(nn*5 + n));
    }
    writeline('-', 79, false, true);
    for (size_t j = 0; j < m_points; j++) {
        writelog("\n {:10.4g} ", m_z[j]);
        for (size_t n = 0; n < nrem; n++) {
            writelog(" {:10.4g} ", x[index(nn*5 + n, j)]);
        }
    }
    writelog("\n");
}

void Domain1D::setProfile(const string& name, double* values, double* soln)
{
    warn_deprecated(
        "Domain1D::setProfile", "To be removed after Cantera 3.2. Replaceable by "
        "version using vector arguments.");
    for (size_t n = 0; n < m_nv; n++) {
        if (name == componentName(n)) {
            for (size_t j = 0; j < m_points; j++) {
                soln[index(n, j) + m_iloc] = values[j];
            }
            return;
        }
    }
    throw CanteraError("Domain1D::setProfile", "unknown component: "+name);
}

void Domain1D::setRefineCriteria(double ratio, double slope, double curve, double prune)
{
    m_refiner->setCriteria(ratio, slope, curve, prune);
}

vector<double> Domain1D::getRefineCriteria()
{
    return m_refiner->getCriteria();
}

void Domain1D::_getInitialSoln(double* x)
{
    for (size_t j = 0; j < m_points; j++) {
        for (size_t n = 0; n < m_nv; n++) {
            x[index(n,j)] = initialValue(n,j);
        }
    }
}

double Domain1D::initialValue(size_t n, size_t j)
{
    throw NotImplementedError("Domain1D::initialValue");
}

} // namespace
