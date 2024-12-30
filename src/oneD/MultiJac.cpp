//! @file MultiJac.cpp Implementation file for class MultiJac

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/MultiJac.h"
#include <ctime>

namespace Cantera
{

MultiJac::MultiJac(OneDim& r)
    : m_mat(r.size(),r.bandwidth(),r.bandwidth())
    , m_n(r.size())
{
    m_resid = &r;
    m_r1.resize(m_n);
    m_ssdiag.resize(m_n);
    m_mask.resize(m_n);
}

void MultiJac::reset() {
    m_mat.bfill(0.0);
    m_age = 10000;
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
    for (size_t n = 0; n < m_n; n++) {
        m_mat.value(n,n) = m_ssdiag[n] - mask[n]*rdt;
    }
}

void MultiJac::eval(double* x0, double* resid0, double rdt)
{
    warn_deprecated("MultiJac::eval", "To be removed after Cantera 3.2. "
        "Jacobian evaluation moved to OneDim::evalJacobian().");
    m_resid->evalJacobian(x0);
}

} // namespace
