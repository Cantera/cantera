/**
 * @file Domain1D.cpp
 *
 */

#include "Domain1D.h"

using namespace std;

namespace Cantera {

    void Domain1D::
    setTolerances(int nr, const doublereal* rtol, 
        int na, const doublereal* atol, int ts) {
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
    setTolerances(int n, doublereal rtol, doublereal atol, int ts) {
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
    setTolerances(doublereal rtol, doublereal atol,int ts) {
        for (int n = 0; n < m_nv; n++){
            if(ts >= 0) {
                m_rtol_ss[n] = rtol;
                m_atol_ss[n] = atol;
            }
            if (ts <= 0) {
                m_rtol_ts[n] = rtol;
                m_atol_ts[n] = atol;
            }
        }
    }

    void Domain1D::
    setTolerancesTS(doublereal rtol, doublereal atol) {
        for (int n = 0; n < m_nv; n++){
            m_rtol_ts[n] = rtol;
            m_atol_ts[n] = atol;
        }
    }

    void Domain1D::
    setTolerancesSS(doublereal rtol, doublereal atol) {
        for (int n = 0; n < m_nv; n++){
            m_rtol_ss[n] = rtol;
            m_atol_ss[n] = atol;
        }
    }
    
    void Domain1D::
    eval(int jg, doublereal* xg, doublereal* rg, 
        integer* mask, doublereal rdt) {

        if (jg >=0 && (jg < firstPoint() - 1 || jg > lastPoint() + 1)) return;

        // if evaluating a Jacobian, compute the steady-state residual
        if (jg >= 0) rdt = 0.0;

        // start of local part of global arrays
        doublereal* x = xg + loc();
        doublereal* rsd = rg + loc();
        integer* diag = mask + loc();
        
        int jmin, jmax, jpt,  j, i;
        jpt = jg - firstPoint();
        
        if (jg < 0) {      // evaluate all points
            jmin = 0;
            jmax = m_points - 1;
        }
        else {            // evaluate points for Jacobian
            jmin = std::max(jpt-1, 0);
            jmax = std::min(jpt+1,m_points-1);
        }

        for (j = jmin; j <= jmax; j++) {
            if (j == 0 || j == m_points - 1) {
                for (i = 0; i < m_nv; i++) {
                    rsd[index(i,j)] = residual(x,i,j);
                    diag[index(i,j)] = 0;
                }
            }
            else {
                for (i = 0; i < m_nv; i++) {
                    rsd[index(i,j)] = residual(x,i,j) 
                                      - timeDerivativeFlag(i)*rdt*(value(x,i,j) - prevSoln(i,j));
                    diag[index(i,j)] = timeDerivativeFlag(i);
                }
            }
        }
    }
    

    // called to set up initial grid, and after grid refinement
    void Domain1D::setupGrid(int n, const doublereal* z) {
        if (n > 1) {
            resize(m_nv, n);
            int j;
            for (j = 0; j < m_points; j++) m_z[j] = z[j];
        }
    }
    

    void drawline() {
        writelog("\n-------------------------------------"
            "------------------------------------------");
    }
    
    
    /**
     * Print the solution.
     */
    void Domain1D::showSolution(const doublereal* x) {
        int nn = m_nv/5;
        int i, j, n;
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
        int nrem = m_nv - 5*nn;
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

    
    // initial solution
    void Domain1D::_getInitialSoln(doublereal* x) {
        for (int j = 0; j < m_points; j++) {
            for (int n = 0; n < m_nv; n++) {
                x[index(n,j)] = initialValue(n,j);
            }
        }
    }

    doublereal Domain1D::initialValue(int n, int j) { 
        throw CanteraError("Domain1D::initialValue",
            "base class method called!");
    }
    

} // namespace
