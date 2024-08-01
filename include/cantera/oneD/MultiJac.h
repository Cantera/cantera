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
    MultiJac(OneDim& r);

    /**
     * Evaluate the Jacobian at x0. The unperturbed residual function is resid0,
     * which must be supplied on input. The third parameter 'rdt' is the
     * reciprocal of the time step. If zero, the steady-state Jacobian is
     * evaluated.
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

    void incrementDiagonal(int j, double d);

protected:
    //! Residual evaluator for this Jacobian
    /*!
     * This is a pointer to the residual evaluator. This object isn't owned by
     * this Jacobian object.
     */
    OneDim* m_resid;

    vector<double> m_r1;
    double m_rtol = 1e-5;
    double m_atol = sqrt(std::numeric_limits<double>::epsilon());
    double m_elapsed = 0.0;
    vector<double> m_ssdiag;

    //! Transient mask for transient terms, 1 if transient, 0 if steady-state
    vector<int> m_mask;
    int m_nevals = 0;
    int m_age = 100000;
    size_t m_size;
    size_t m_points;
};
}

#endif
