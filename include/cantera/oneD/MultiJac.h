//! @file MultiJac.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MULTIJAC_H
#define CT_MULTIJAC_H

#include "cantera/numerics/BandMatrix.h"
#include "OneDim.h"

namespace Cantera
{

/**
 * Class MultiJac evaluates the Jacobian of a system of equations defined by a
 * residual function supplied by an instance of class OneDim. The residual
 * function may consist of several linked 1D domains, with different variables
 * in each domain.
 * @ingroup onedUtilsGroup
 * @ingroup derivGroup
 */
class MultiJac : public BandMatrix
{
public:
    //! Constructor.
    //! @param r  The nonlinear system for which to compute the Jacobian.
    MultiJac(OneDim& r);

    /**
     * Evaluates the Jacobian at x0 using finite differences. The unperturbed residual
     * function is resid0, which must be supplied on input. The parameter 'rdt' is the
     * reciprocal of the time step. If zero, the steady-state Jacobian is evaluated,
     * otherwise the transient Jacobian is evaluated.
     *
     * Each variable in the input vector `x0` is perturbed to compute the
     * corresponding column of the Jacobian matrix. The Jacobian is computed by
     * perturbing each variable, evaluating the residual function, and then
     * estimating the partial derivatives numerically using finite differences.
     *
     * @param x0  Solution vector at which to evaluate the Jacobian
     * @param resid0  Residual vector at x0
     * @param rdt  Reciprocal of the time step
     */
    void eval(double* x0, double* resid0, double rdt);

    //! Elapsed CPU time spent computing the Jacobian.
    double elapsedTime() const {
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

    //! Update the transient terms in the Jacobian by using the transient mask.
    void updateTransient(double rdt, integer* mask);

    //! Set the Jacobian age.
    void setAge(int age) {
        m_age = age;
    }

    //! Return the transient mask.
    vector<int>& transientMask() {
        return m_mask;
    }

    //! @deprecated To be removed after Cantera 3.1.
    void incrementDiagonal(int j, double d);

protected:
    //! Residual evaluator for this Jacobian
    /*!
     * This is a pointer to the residual evaluator. This object isn't owned by
     * this Jacobian object.
     */
    OneDim* m_resid;

    vector<double> m_r1; //!< Perturbed residual vector
    double m_rtol = 1e-5; //!< Relative tolerance for perturbing solution components

    //! Absolute tolerance for perturbing solution components
    double m_atol = sqrt(std::numeric_limits<double>::epsilon());

    double m_elapsed = 0.0; //!< Elapsed CPU time taken to compute the Jacobian
    vector<double> m_ssdiag; //!< Diagonal of the steady-state Jacobian

    //! Transient mask for transient terms, 1 if transient, 0 if steady-state
    vector<int> m_mask;

    int m_nevals = 0; //!< Number of Jacobian evaluations.
    int m_age = 100000; //!< Age of the Jacobian (times incrementAge() has been called)
};
}

#endif
