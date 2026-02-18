//! @file MultiJac.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MULTIJAC_H
#define CT_MULTIJAC_H

#include "cantera/numerics/BandMatrix.h"
#include "cantera/numerics/SystemJacobian.h"
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
class MultiJac : public SystemJacobian
{
public:
    MultiJac() = default;
    void reset() override;
    const string type() const override { return "banded-direct"; }
    void setValue(size_t row, size_t col, double value) override;
    void initialize(size_t nVars) override;
    void setBandwidth(size_t bw) override;
    void updateTransient(double rdt, span<const int> mask) override;
    void factorize() override {
        m_mat.factor();
    }

    double& value(size_t i, size_t j) {
        return m_mat.value(i, j);
    }

    double value(size_t i, size_t j) const {
        return m_mat.value(i, j);
    }

    void solve(span<const double> b, span<double> x) override {
        m_mat.solve(b, x);
    }

    int info() const override {
        return m_mat.info();
    }

protected:
    BandMatrix m_mat; //!< Underlying matrix storage
    vector<double> m_ssdiag; //!< Diagonal of the steady-state Jacobian
};

}

#endif
