/**
 *
 *  @file Jac1D.h
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

#ifndef CT_JAC1D_H
#define CT_JAC1D_H

#include "Domain1D.h"
#include "BandMatrix.h"
//#include "ArrayViewer.h"
#include "Array.h"
#include "time.h"

namespace Cantera {

    /**
     * Class Jac1D evaluates the Jacobian of a system of equations
     * defined by a residual function of class Domain1D. It is
     * assumed that the Jacobian is banded.
     */
    class Jac1D : public BandMatrix {

    public:

        /** 
         * Constructor. The residual function defining the system of
         * equations must be supplied.
         */
        Jac1D(Domain1D& r);

        /// Destructor. Does nothing.
        virtual ~Jac1D(){}

        /**
         * Evaluate the Jacobian.
         */
        void eval(doublereal* x0, doublereal* resid0);

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

        Domain1D* m_resid;
        Array2D m_r1;
        //        ArrayViewer m_x0, m_r0;
        int m_nv, m_points;
        doublereal m_atol;
        doublereal m_elapsed;
        int m_nevals;
        int m_age;

    private:
        size_t index(int m, int j) { return m_nv*j + m; }
    };
}

#endif


