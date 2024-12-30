//! @file MultiJac.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MULTIJAC_H
#define CT_MULTIJAC_H

#include "cantera/numerics/BandMatrix.h"
#include "cantera/numerics/PreconditionerBase.h"
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
class MultiJac : public PreconditionerBase
{
public:
    //! Constructor.
    //! @param r  The nonlinear system for which to compute the Jacobian.
    //! @deprecated To be removed after %Cantera 3.2. Use default constructor instead.
    MultiJac(OneDim& r);

    MultiJac() = default;

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
     * @deprecated To be removed after %Cantera 3.2. Jacobian evaluation moved to
     *     OneDim::evalJacobian().
     */
    void eval(double* x0, double* resid0, double rdt);

    void reset() override;
    void setValue(size_t row, size_t col, double value) override;
    void initialize(size_t nVars) override;
    void setBandwidth(size_t bw) override;
    void updateTransient(double rdt, integer* mask) override;

    double& value(size_t i, size_t j) {
        return m_mat.value(i, j);
    }

    double value(size_t i, size_t j) const {
        return m_mat.value(i, j);
    }

    void solve(const double* const b, double* const x) {
        m_mat.solve(b, x);
    }

    void solve(const size_t stateSize, double* b, double* x) override {
        m_mat.solve(b, x);
    }

    int info() const override {
        return m_mat.info();
    }

    //! Return the transient mask.
    //! @deprecated Unused. To be removed after %Cantera 3.2.
    vector<int>& transientMask() {
        warn_deprecated("MultiJac::transientMask", "To be removed after Cantera 3.2");
        return m_mask;
    }

protected:
    //! Residual evaluator for this Jacobian
    /*!
     * This is a pointer to the residual evaluator. This object isn't owned by
     * this Jacobian object.
     *
     * @deprecated Unused. To be removed after %Cantera 3.2.
     */
    OneDim* m_resid = nullptr;

    BandMatrix m_mat; //!< Underlying matrix storage
    vector<double> m_ssdiag; //!< Diagonal of the steady-state Jacobian

    //! Transient mask for transient terms, 1 if transient, 0 if steady-state.
    //! @deprecated Unused. To be removed after %Cantera 3.2.
    vector<int> m_mask;
};

}

#endif
