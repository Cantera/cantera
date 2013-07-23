/**
 * @file Domain1D.cpp
 */

#include "cantera/oneD/Domain1D.h"
 #include "cantera/base/ctml.h"

#include <cstdio>

using namespace std;
using namespace ctml;

namespace Cantera
{

void Domain1D::
setTolerances(size_t nr, const doublereal* rtol,
              size_t na, const doublereal* atol, int ts)
{
    warn_deprecated("Domain1D::setTolerances",
                    "Use setTransientTolerances or setSteadyTolerances.");
    if (nr < m_nv || na < m_nv)
        throw CanteraError("Domain1D::setTolerances",
                           "wrong array size for solution error tolerances. "
                           "Size should be at least "+int2str(m_nv));
    if (ts >= 0) {
        copy(rtol, rtol + m_nv, m_rtol_ss.begin());
        copy(atol, atol + m_nv, m_atol_ss.begin());
    }
    if (ts <= 0) {
        copy(rtol, rtol + m_nv, m_rtol_ts.begin());
        copy(atol, atol + m_nv, m_atol_ts.begin());
    }
}

void Domain1D::
setTolerances(size_t n, doublereal rtol, doublereal atol, int ts)
{
    warn_deprecated("Domain1D::setTolerances",
                    "Use setTransientTolerances or setSteadyTolerances.");
    if (ts >= 0) {
        m_rtol_ss[n] = rtol;
        m_atol_ss[n] = atol;
    }
    if (ts <= 0) {
        m_rtol_ts[n] = rtol;
        m_atol_ts[n] = atol;
    }
}

void Domain1D::
setTolerances(doublereal rtol, doublereal atol,int ts)
{
    warn_deprecated("Domain1D::setTolerances",
                    "Use setTransientTolerances or setSteadyTolerances.");
    for (size_t n = 0; n < m_nv; n++) {
        if (ts >= 0) {
            m_rtol_ss[n] = rtol;
            m_atol_ss[n] = atol;
        }
        if (ts <= 0) {
            m_rtol_ts[n] = rtol;
            m_atol_ts[n] = atol;
        }
    }
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

void Domain1D::setTolerancesTS(doublereal rtol, doublereal atol, size_t n)
{
    warn_deprecated("Domain1D::setTolerancesTS",
                    "Use setTransientTolerances");
    setTransientTolerances(rtol, atol, n);
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

void Domain1D::setTolerancesSS(doublereal rtol, doublereal atol, size_t n)
{
    warn_deprecated("Domain1D::setTolerancesSS",
                    "Use setSteadyTolerances");
    setSteadyTolerances(rtol, atol, n);
}

void Domain1D::
eval(size_t jg, doublereal* xg, doublereal* rg,
     integer* mask, doublereal rdt)
{

    if (jg != npos && (jg + 1 < firstPoint() || jg > lastPoint() + 1)) {
        return;
    }

    // if evaluating a Jacobian, compute the steady-state residual
    if (jg != npos) {
        rdt = 0.0;
    }

    // start of local part of global arrays
    doublereal* x = xg + loc();
    doublereal* rsd = rg + loc();
    integer* diag = mask + loc();

    size_t jmin, jmax, jpt,  j, i;
    jpt = jg - firstPoint();

    if (jg == npos) {      // evaluate all points
        jmin = 0;
        jmax = m_points - 1;
    } else {          // evaluate points for Jacobian
        jmin = std::max<size_t>(jpt, 1) - 1;
        jmax = std::min(jpt+1,m_points-1);
    }

    for (j = jmin; j <= jmax; j++) {
        if (j == 0 || j == m_points - 1) {
            for (i = 0; i < m_nv; i++) {
                rsd[index(i,j)] = residual(x,i,j);
                diag[index(i,j)] = 0;
            }
        } else {
            for (i = 0; i < m_nv; i++) {
                rsd[index(i,j)] = residual(x,i,j)
                                  - timeDerivativeFlag(i)*rdt*(value(x,i,j) - prevSoln(i,j));
                diag[index(i,j)] = timeDerivativeFlag(i);
            }
        }
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
    vector<XML_Node*> nodes;
    dom.getChildren("floatArray", nodes);
    for (size_t i = 0; i < nodes.size(); i++) {
        string title = nodes[i]->attrib("title");
        getFloatArray(*nodes[i], values, false);
        if (values.size() != nComponents()) {
            if (loglevel > 0) {
                writelog("Warning: Domain1D::restore: Got an array of length " +
                         int2str(values.size()) + " when one of length " +
                         int2str(nComponents()) + " was expected. " +
                         "Tolerances for individual species may not be preserved.\n");
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

void Domain1D::setupGrid(size_t n, const doublereal* z)
{
    if (n > 1) {
        resize(m_nv, n);
        for (size_t j = 0; j < m_points; j++) {
            m_z[j] = z[j];
        }
    }
}

void drawline()
{
    writelog("\n-------------------------------------"
             "------------------------------------------");
}

void Domain1D::showSolution(const doublereal* x)
{
    size_t nn = m_nv/5;
    size_t i, j, n;
    //char* buf = new char[100];
    char buf[100];
    doublereal v;
    for (i = 0; i < nn; i++) {
        drawline();
        sprintf(buf, "\n        z   ");
        writelog(buf);
        for (n = 0; n < 5; n++) {
            sprintf(buf, " %10s ",componentName(i*5 + n).c_str());
            writelog(buf);
        }
        drawline();
        for (j = 0; j < m_points; j++) {
            sprintf(buf, "\n %10.4g ",m_z[j]);
            writelog(buf);
            for (n = 0; n < 5; n++) {
                v = value(x, i*5+n, j);
                sprintf(buf, " %10.4g ",v);
                writelog(buf);
            }
        }
        writelog("\n");
    }
    size_t nrem = m_nv - 5*nn;
    drawline();
    sprintf(buf, "\n        z   ");
    writelog(buf);
    for (n = 0; n < nrem; n++) {
        sprintf(buf, " %10s ", componentName(nn*5 + n).c_str());
        writelog(buf);
    }
    drawline();
    for (j = 0; j < m_points; j++) {
        sprintf(buf, "\n %10.4g ",m_z[j]);
        writelog(buf);
        for (n = 0; n < nrem; n++) {
            v = value(x, nn*5+n, j);
            sprintf(buf, " %10.4g ", v);
            writelog(buf);
        }
    }
    writelog("\n");
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
    throw CanteraError("Domain1D::initialValue",
                       "base class method called!");
    return 0.0;
}

} // namespace
