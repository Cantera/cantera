/**
 * @file Domain1D.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/Domain1D.h"
#include "cantera/oneD/MultiJac.h"
#include "cantera/oneD/refine.h"
#include "cantera/base/AnyMap.h"
#include "cantera/base/SolutionArray.h"

namespace Cantera
{

Domain1D::Domain1D(size_t nv, size_t points, double time)
{
    resize(nv, points);
}

Domain1D::~Domain1D()
{
}

int Domain1D::domainType()
{
    warn_deprecated("Domain1D::domainType",
        "To be changed after Cantera 3.0; for new behavior, see 'type'.");
    return m_type;
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

std::string Domain1D::componentName(size_t n) const
{
    if (m_name[n] != "") {
        return m_name[n];
    } else {
        return fmt::format("component {}", n);
    }
}

size_t Domain1D::componentIndex(const std::string& name) const
{
    for (size_t n = 0; n < nComponents(); n++) {
        if (name == componentName(n)) {
            return n;
        }
    }
    throw CanteraError("Domain1D::componentIndex",
                       "no component named "+name);
}

void Domain1D::setTransientTolerances(doublereal rtol, doublereal atol, size_t n)
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

void Domain1D::setSteadyTolerances(doublereal rtol, doublereal atol, size_t n)
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
        m_container->jacobian().setAge(10000);
        m_container->saveStats();
    }
}

AnyMap Domain1D::getMeta() const
{
    auto wrap_tols = [this](const vector_fp& tols) {
        // If all tolerances are the same, just store the scalar value.
        // Otherwise, store them by component name
        std::set<double> unique_tols(tols.begin(), tols.end());
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

AnyMap Domain1D::serialize(const double* soln) const
{
    warn_deprecated("Domain1D::serialize",
        "To be removed after Cantera 3.0; superseded by asArray.");
    AnyMap out;
    auto arr = asArray(soln);
    arr->writeEntry(out, "", "");
    return out;
}

shared_ptr<SolutionArray> Domain1D::toArray(bool normalize) const {
    if (!m_state) {
        throw CanteraError("Domain1D::toArray",
            "Domain needs to be installed in a container before calling asArray.");
    }
    auto ret = asArray(m_state->data() + m_iloc);
    if (normalize) {
        ret->normalize();
    }
    return ret;
}

void Domain1D::fromArray(const shared_ptr<SolutionArray>& arr)
{
    if (!m_state) {
        throw CanteraError("Domain1D::fromArray",
            "Domain needs to be installed in a container before calling fromArray.");
    }
    resize(nComponents(), arr->size());
    m_container->resize();
    fromArray(*arr, m_state->data() + m_iloc);
    _finalize(m_state->data() + m_iloc);
}

void Domain1D::setMeta(const AnyMap& meta)
{
    auto set_tols = [&](const AnyValue& tols, const string& which, vector_fp& out)
    {
        if (!tols.hasKey(which)) {
            return;
        }
        const auto& tol = tols[which];
        if (tol.isScalar()) {
            out.assign(nComponents(), tol.asDouble());
        } else {
            for (size_t i = 0; i < nComponents(); i++) {
                std::string name = componentName(i);
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

void Domain1D::restore(const AnyMap& state, double* soln, int loglevel)
{
    warn_deprecated("Domain1D::restore",
        "To be removed after Cantera 3.0; restore from SolutionArray instead.");
    auto arr = SolutionArray::create(solution());
    arr->readEntry(state, "", "");
    fromArray(*arr, soln);
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

void Domain1D::setupGrid(size_t n, const doublereal* z)
{
    if (n > 1) {
        resize(m_nv, n);
        for (size_t j = 0; j < m_points; j++) {
            m_z[j] = z[j];
        }
    }
}

void Domain1D::showSolution_s(std::ostream& s, const double* x)
{
    warn_deprecated("Domain1D::showSolution_s",
        "To be removed after Cantera 3.0; replaced by 'show'.");
    show(s, x);
}

void Domain1D::showSolution(const double* x)
{
    warn_deprecated("Domain1D::showSolution",
        "To be removed after Cantera 3.0; replaced by 'show'.");
    show(x);
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
                double v = value(x, i*5+n, j);
                writelog(" {:10.4g} ", v);
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
            double v = value(x, nn*5+n, j);
            writelog(" {:10.4g} ", v);
        }
    }
    writelog("\n");
}

void Domain1D::setProfile(const std::string& name, double* values, double* soln)
{
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

void Domain1D::_getInitialSoln(doublereal* x)
{
    for (size_t j = 0; j < m_points; j++) {
        for (size_t n = 0; n < m_nv; n++) {
            x[index(n,j)] = initialValue(n,j);
        }
    }
}

doublereal Domain1D::initialValue(size_t n, size_t j)
{
    throw NotImplementedError("Domain1D::initialValue");
}

} // namespace
