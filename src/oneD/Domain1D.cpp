/**
 * @file Domain1D.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/Domain1D.h"
#include "cantera/oneD/MultiJac.h"
#include "cantera/base/ctml.h"

using namespace std;

namespace Cantera
{

Domain1D::Domain1D(size_t nv, size_t points, double time) :
    m_rdt(0.0),
    m_nv(0),
    m_container(0),
    m_index(npos),
    m_type(0),
    m_iloc(0),
    m_jstart(0),
    m_left(0),
    m_right(0),
    m_bw(-1),
    m_force_full_update(false)
{
    resize(nv, points);
}

void Domain1D::resize(size_t nv, size_t np)
{
    // if the number of components is being changed, then a
    // new grid refiner is required.
    if (nv != m_nv || !m_refiner) {
        m_nv = nv;
        m_refiner.reset(new Refiner(*this));
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

XML_Node& Domain1D::save(XML_Node& o, const doublereal* const sol)
{
    XML_Node& d = o.addChild("domain");
    d.addAttribute("points", nPoints());
    d.addAttribute("components", nComponents());
    d.addAttribute("id", id());
    addFloatArray(d, "abstol_transient", nComponents(), &m_atol_ts[0]);
    addFloatArray(d, "reltol_transient", nComponents(), &m_rtol_ts[0]);
    addFloatArray(d, "abstol_steady", nComponents(), &m_atol_ss[0]);
    addFloatArray(d, "reltol_steady", nComponents(), &m_rtol_ss[0]);
    return d;
}

void Domain1D::restore(const XML_Node& dom, doublereal* soln, int loglevel)
{
    vector_fp values;
    vector<XML_Node*> nodes = dom.getChildren("floatArray");
    for (size_t i = 0; i < nodes.size(); i++) {
        string title = nodes[i]->attrib("title");
        getFloatArray(*nodes[i], values, false);
        if (values.size() != nComponents()) {
            if (loglevel > 0) {
                warn_user("Domain1D::restore", "Received an array of length "
                    "{} when one of length {} was expected. Tolerances for "
                    "individual species may not be preserved.",
                    values.size(), nComponents());
            }
            // The number of components will differ when restoring from a
            // mechanism with a different number of species. Assuming that
            // tolerances are the same for all species, we can just copy the
            // tolerance from the last species.
            if (!values.empty()) {
                values.resize(nComponents(), values[values.size()-1]);
            } else {
                // If the tolerance vector is empty, just leave the defaults
                // in place.
                continue;
            }
        }
        if (title == "abstol_transient") {
            m_atol_ts = values;
        } else if (title == "reltol_transient") {
            m_rtol_ts = values;
        } else if (title == "abstol_steady") {
            m_atol_ss = values;
        } else if (title == "reltol_steady") {
            m_rtol_ss = values;
        } else {
            throw CanteraError("Domain1D::restore",
                               "Got an unexpected array, '" + title + "'");
        }
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

void Domain1D::setupGrid(size_t n, const doublereal* z)
{
    if (n > 1) {
        resize(m_nv, n);
        for (size_t j = 0; j < m_points; j++) {
            m_z[j] = z[j];
        }
    }
}

void Domain1D::showSolution(const doublereal* x)
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
