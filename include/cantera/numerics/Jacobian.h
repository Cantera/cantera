//! @file MultiJac.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_JACOBIAN_H
#define CT_JACOBIAN_H

// #include "BandMatrix.h"
#include "DenseMatrix.h"
#include "FuncEval.h"

namespace Cantera
{

/**
 * Class MultiJac evaluates the Jacobian of a system of equations defined by a
 * residual function supplied by an instance of class OneDim. The residual
 * function may consist of several linked 1D domains, with different variables
 * in each domain.
 * @ingroup onedim
 */
class Jacobian : public DenseMatrix
// class Jacobian : public BandMatrix
{
public:
    Jacobian(FuncEval& func);

    /**
     * Evaluate the Jacobian at x0. The unperturbed residual function is resid0,
     * which must be supplied on input. The third parameter 'rdt' is the
     * reciprocal of the time step. If zero, the steady-state Jacobian is
     * evaluated.
     */
    void eval(doublereal* x0, doublereal* resid0, double rdt);

    //wrapper for DenseMatrix solve
    void compute(doublereal* step){
        solve(*this, step);
    }

    //wrapper for DenseMatrix multiply
    void mult(doublereal* b, doublereal* result){
        multiply(*this, b, result);
    }

    //! Elapsed CPU time spent computing the Jacobian.
    doublereal elapsedTime() const {
        return m_elapsed;
    }

    //! Number of Jacobian evaluations.
    int nEvals() const {
        return m_nevals;
    }

    //! Number of times 'incrementAge' has been called since the last evaluation
    int age() const {
        return m_age;
    }

    //! Increment the Jacobian age.
    void incrementAge() {
        m_age++;
    }

    void updateTransient(doublereal rdt, integer* mask);

    //! Set the Jacobian age.
    void setAge(int age) {
        m_age = age;
    }

    vector_int& transientMask() {
        return m_mask;
    }

    void incrementDiagonal(int j, doublereal d);

protected:
    //! Residual evaluator for this Jacobian
    /*!
     * This is a pointer to the residual evaluator. This object isn't owned by
     * this Jacobian object.
     */
    FuncEval* m_residfunc;

    vector_fp m_r1;
    doublereal m_rtol, m_atol;
    doublereal m_elapsed;
    vector_fp m_ssdiag;
    vector_int m_mask;
    int m_nevals;
    int m_age;
    size_t m_size;
    size_t m_points;
};
}

#endif
