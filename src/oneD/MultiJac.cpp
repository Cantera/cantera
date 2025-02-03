//! @file MultiJac.cpp Implementation file for class MultiJac

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/MultiJac.h"
#include <ctime>

namespace Cantera
{

MultiJac::MultiJac(OneDim& r)
{
    warn_deprecated("MultiJac::MultiJac(OneDim&)",
                    "Non-default constructor to be removed after Cantera 3.2.");
    m_resid = &r;
    initialize(m_dim);
    setBandwidth(r.bandwidth());
}

void MultiJac::reset()
{
    m_mat.bfill(0.0);
    m_age = 10000;
}

void MultiJac::initialize(size_t nVars)
{
    m_dim = nVars;
    m_mat.resize(m_dim, m_mat.nSubDiagonals(), m_mat.nSuperDiagonals());
    m_ssdiag.resize(m_dim);
    m_mask.resize(m_dim);
}

void MultiJac::setBandwidth(size_t bw)
{
    m_mat.resize(m_dim, bw, bw);
}

void MultiJac::setValue(size_t row, size_t col, double value)
{
    m_mat(row, col) = value;
    if (row == col) {
        m_ssdiag[row] = value;
    }
}

void MultiJac::updateTransient(double rdt, integer* mask)
{
    for (size_t n = 0; n < m_dim; n++) {
        m_mat.value(n,n) = m_ssdiag[n] - mask[n]*rdt;
    }
    factorize();
}

void MultiJac::eval(double* x0, double* resid0, double rdt)
{
    warn_deprecated("MultiJac::eval", "To be removed after Cantera 3.2. "
        "Jacobian evaluation moved to OneDim::evalJacobian().");
    if (!m_resid) {
        throw CanteraError("MultiJac::eval", "Can only be used in combination with "
            "(deprecated) constructor that takes a OneDim object.");
    }
    m_resid->evalJacobian(x0);
}

} // namespace
