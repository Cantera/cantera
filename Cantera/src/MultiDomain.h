deprecated
/**
 *
 *  @file Resid1D.h
 *
 *  >>>>> Under construction! <<<<<
 *
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2002 California Institute of Technology
 *
 */

#ifndef CT_MULTIDOM_H
#define CT_MULTIDOM_H


#include "stringUtils.h"
//#include "Array.h"
#include "Resid1D.h"

namespace Cantera {


    /**
     * Residual function evaluator for a one-dimensional problem.
     */
    class MultiDomain {
    public:

        /**
         * Constructor.
         * @param nv Number of variables at each grid point.
         * @param points Number of grid points.
         */
        MultiDomain() : m_bw(0), m_nd(0), m_rdt(0.0), m_jac_ok(false) {}

        /// Destructor.
        virtual ~MultiDomain(){}

        int nDomains() const { return m_nd; }

        Resid1D& domain(int i) const { return *m_dom[i]; }
        int start(int i) const { return m_start[i]; }
        int size() const { return m_size; }
        int bandwidth() const { return m_bw; }

        virtual void addDomain(Resid1D* d) {
            m_dom.push_back(d);
            int sz = d->nPoints() 
                     * d->nComponents();
            if (m_nd > 0) {
                m_start.push_back(m_start.back() + m_states.back());
            }
            else
                m_start.push_back(0);
            m_states.push_back(sz);
            m_comp.push_back(d->nComponents());
            m_points.push_back(d->nPoints());
            int bw1, bw2 = 0;
            bw1 = 2*d->nComponents() - 1;
            if (m_nd > 0) {
                bw2 = d->nComponents() + m_dom[m_nd-1]->nComponents() - 1;
            }
            if (bw1 > m_bw) m_bw = bw1;
            if (bw2 > m_bw) m_bw = bw2;
            m_nd = m_states.size();
            m_size = m_start.back() + m_states.back();
        }

        void evalDomain(int i, int j, doublereal* x, doublereal* r, 
            doublereal rdt) {
            int jpt;
            if (j < 0) 
                jpt = j;
            else
                jpt = (j - m_start[i])/m_comp[i];
            //cout << "calling Resid1D::eval. jpt, start = " << jpt << "  " << m_start[i] << endl;
            m_dom[i]->eval(jpt, x + m_start[i], r + m_start[i], rdt);
        } 

        
        doublereal ssnorm(doublereal* x, doublereal* r) {
            //cout << " calling eval to get ss norm " << endl;
            eval(-1, x, r, 0.0);
            doublereal ss = 0.0;
            for (int i = 0; i < m_size; i++) 
                ss = fmax(fabs(r[i]),ss);
            return ss;
        }
            
        doublereal rdt() const { return m_rdt; }
        
        void initTimeInteg(doublereal dt, doublereal* x) {
            doublereal rdt_old = m_rdt;
            m_rdt = 1.0/dt;
            if (fabs(rdt_old - m_rdt) > Tiny) {
                m_jac_ok = false;
            }
            int i;
            for (i = 0; i < m_nd; i++) 
                m_dom[i]->initTimeInteg(dt, x + m_start[i]);
        }

        bool transient() const { return (m_rdt != 0.0);}
        bool steady() const { return (m_rdt == 0.0); }

        void setSteadyMode() {
            if (m_rdt > 0)
                m_jac_ok = false;
            m_rdt = 0.0;
        }
 
        void eval(int j, double* x, double* r, doublereal rdt=-1.0) {
            int i, jpt;
            if (rdt < 0.0) rdt = m_rdt;
            if (j < 0) {
                for (i = 0; i < m_nd; i++) {
                    evalDomain(i, j, x, r, rdt);
                }
            }
            else {
                for (i = 0; i < m_nd; i++) {
                    if (j >= m_start[i]) {
                        //cout << "calling evalDomain with j = " << j << endl;
                        evalDomain(i, j, x, r, rdt);
                        jpt = j - m_start[i];
                        if (jpt == 0 && i > 0) 
                            evalDomain(i-1, j-1, x, r, rdt);
                        else if (jpt == m_states[i] - 1 && i < m_nd - 1) 
                            evalDomain(i+1, j+1, x, r, rdt);
                        break;
                    }
                }
            }
        }

    protected:
        doublereal m_rdt;
        bool m_jac_ok;
        int m_nd, m_bw, m_size;
        vector_int m_states;
        vector_int m_start;
        vector_int m_comp, m_points;
        vector<Resid1D*> m_dom;

    private:

    };
}

#endif


