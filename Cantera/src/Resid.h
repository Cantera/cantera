/**
 *
 *  @file Resid.h
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

#ifndef CT_RESID_H
#define CT_RESID_H

#include <vector>
#include "ctexceptions.h"
#include "stringUtils.h"

namespace Cantera {


    /**
     * Residual function evaluator for a zero-dimensional problem.
     */
    class Resid {
    public:

        /**
         * Constructor.
         * @param nv Number of variables at each grid point.
         * @param points Number of grid points.
         */
        Resid(int nv=1, doublereal time = 0.0) {
            m_nv = nv;
            m_max.resize(m_nv, 0.0);
            m_min.resize(m_nv, 0.0);
            m_rtol.resize(m_nv, 0.0);
            m_atol.resize(m_nv, 0.0);
            m_time = time;
            m_slast.resize(m_nv);
            setSteadyMode();
        }

        void resize(int nv) {
            m_nv = nv;
            m_max.resize(m_nv, 0.0);
            m_min.resize(m_nv, 0.0);
            m_rtol.resize(m_nv, 0.0);
            m_atol.resize(m_nv, 0.0);
            m_slast.resize(m_nv);
            setSteadyMode();
        }

        /// Destructor.
        virtual ~Resid(){}

        /// Number of components
        int nComponents() const { return m_nv; }

        /// Name of the nth component.
        virtual string componentName(int n) const { 
            return "component " + int2str(n); }

        void setBounds(int nl, const doublereal* lower, 
            int nu, const doublereal* upper) {
            if (nl != m_nv || nu != m_nv)
                throw CanteraError("Resid::setBounds",
                    "wrong array size for solution bounds");
            copy(upper, upper + m_nv, m_max.begin());
            copy(lower, lower + m_nv, m_min.begin());
        }

        void setTolerances(int nr, const doublereal* rtol, 
            int na, const doublereal* atol) {
            if (nr != m_nv || na != m_nv)
                throw CanteraError("Resid::setTolerances",
                    "wrong array size for solution error tolerances");
            copy(rtol, rtol + m_nv, m_rtol.begin());
            copy(atol, atol + m_nv, m_atol.begin());
        }

        doublereal rtol(int n) { return m_rtol[n]; }
        doublereal atol(int n) { return m_atol[n]; }

        doublereal upperBound(int n) const { return m_max[n]; }
        doublereal lowerBound(int n) const { return m_min[n]; }

        void initTimeInteg(doublereal dt, const doublereal* x0) {
            copy(x0, x0 + m_nv, m_slast.begin());
            m_rdt = 1.0/dt;
        } 
        
        void setSteadyMode() { m_rdt = 0.0; }

        bool steady() { return (m_rdt == 0.0); }
        bool transient() { return (m_rdt != 0.0); }

        
        /**
         * Evaluate the residual function.
         */         
        virtual void eval(doublereal* x, doublereal* r) {
            throw CanteraError("Resid::eval",
                "residual function not defined.");
        }

        virtual void update(doublereal* x) {}

        void evalss(doublereal* x, doublereal* r) {
            doublereal rdt_save = m_rdt;
            m_rdt = 0.0;
            eval(x, r);
            m_rdt = rdt_save;
        }

        doublereal time() { return m_time;}
        doublereal rdt() { return m_rdt; }
        void incrementTime(doublereal dt) { m_time += dt; }

    protected:
        int m_nv;
        int m_points;
        vector_fp m_slast;
        doublereal m_rdt;
        doublereal m_time;
        vector_fp m_max;
        vector_fp m_min;
        vector_fp m_rtol;
        vector_fp m_atol;

    private:

    };
}

#endif


