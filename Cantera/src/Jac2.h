/**
 *
 *  @file Jac2.h
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

#ifndef CT_JAC2_H
#define CT_JAC2_H


#include "BandMatrix.h"
#include "ctlapack.h"
#include "../ext/math/gmres.h"
#include "stringUtils.h"
#include "Array.h"
#include "time.h"

namespace Cantera {


    /**
     * Residual function evaluator for a one-dimensional problem.
     */
    class ResidFunc2 {
    public:

        /**
         * Constructor.
         * @param nv Number of variables at each grid point.
         * @param points Number of grid points.
         */
        ResidFunc2(int nv, int points) {
            m_nv = nv;
            m_points = points;
            m_max.resize(m_nv, 0.0);
            m_min.resize(m_nv, 0.0);
            m_rtol.resize(m_nv, 0.0);
            m_atol.resize(m_nv, 0.0);
            m_slast.resize(m_nv, m_points);
            m_rdt = 0.0;

            //            m_soln.resize(m_nv, m_points, 0.0);
            //m_resid.resize(m_nv, m_points, 0.0);
        }

        /// Destructor.
        virtual ~ResidFunc2(){}

        /// Number of components at each grid point.
        int nComponents() const { return m_nv; }

        /// Number of grid points.
        int nPoints() const { return m_points; }

        /// Name of the nth component.
        virtual string componentName(int n) const { 
            return "component " + int2str(n); }

        void setBounds(const vector_fp& lower, const vector_fp& upper) {
            if (lower.size() != m_nv || upper.size() != m_nv)
                throw CanteraError("ResidFunc2::setBounds",
                    "wrong array size for solution bounds");
            m_max = upper;
            m_min = lower;
        }

        void setTolerances(vector_fp& rtol, vector_fp& atol) {
            m_rtol = rtol;
            m_atol = atol;
        }

        doublereal rtol(int n) { return m_rtol[n]; }
        doublereal atol(int n) { return m_atol[n]; }

        doublereal upperBound(int n) const {
            return m_max[n];
        }

        doublereal lowerBound(int n) const {
            return m_min[n];
        }


        void initTimeInteg(doublereal dt, const Array2D& x0) {
            m_slast = x0;
            m_rdt = 1.0/dt;
        } 
        
        void setSteadyMode() {
            m_rdt = 0.0;
        }

        bool steady() { return (m_rdt == 0.0); }
        bool transient() { return (m_rdt != 0.0); }

        /// Evaluate the residual function at point j.         
        virtual void eval(int j, Array2D& x, Array2D& r) {
            throw CanteraError("ResidFunc2::eval",
                "residual function not defined.");
        }

        void evalss(Array2D& x, Array2D& r) {
            doublereal rdt_save = m_rdt;
            m_rdt = 0.0;
            eval(-1, x, r);
            m_rdt = rdt_save;
        }

    protected:
        int m_nv;
        int m_points;
        //        Array2D m_soln;
        //Array2D m_resid;
        Array2D m_slast;
        doublereal m_rdt;
        vector_fp m_max;
        vector_fp m_min;
        vector_fp m_rtol;
        vector_fp m_atol;

    private:

    };


    ///////////////////////////////////////////////////////////////


    /**
     * Class Jac2 evaluates the Jacobian of a system of equations
     * defined by a residual function of class ResidFunc2. It is
     * assumed that the Jacobian is banded.
     */
    class Jac2 : public BandMatrix {

    public:

        /** 
         * Constructor. The residual function defining the system of
         * equations must be supplied.
         */
        Jac2(ResidFunc2& r);

        /// Destructor. Does nothing.
        virtual ~Jac2(){}

        /**
         * Evaluate the Jacobian.
         */
        void eval(Array2D& x0, Array2D& resid0);

        /**
         * Returns the matrix element describing the influence of the
         * nth component at point j on the mth equation at point
         * i. Due to the assumption of a banded Jacobian, this will be
         * zero unless |i - j| <= 1.
         */
        doublereal& v(int m, int i, int n, int j) {
            return value(i*m_nv + m, j*m_nv + n);
        }

        doublereal elapsedTime() const {
            return m_elapsed;
        }

        int nEvals() const { return m_nevals; }

        int age() const { return m_age; }

        void incrementAge() { m_age++; }
        void setAge(int age) { m_age = age; }

    protected:

        ResidFunc2* m_resid;
        Array2D m_r1;
        int m_nv, m_points;
        doublereal m_atol;
        doublereal m_elapsed;
        int m_nevals;
        int m_age;

    private:

    };
}

#endif


