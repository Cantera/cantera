//! @file MultiJac.cpp Implementation file for class MultiJac

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/Jacobian.h"
#include "cantera/base/utilities.h"

#include <ctime>

using namespace std;

namespace Cantera
{

// MultiJac::MultiJac(OneDim& r)
Jacobian::Jacobian(FuncEval& func)
    : DenseMatrix((&func)->neq(),(&func)->neq())
    // : BandMatrix((&func)->neq(), 2*(&func)->neq()-1,2*(&func)->neq()-1)
{
    m_residfunc = &func;
    m_size = m_residfunc->neq();
    m_r1.resize(m_size);
    m_ssdiag.resize(m_size);
    m_mask.resize(m_size);
    m_elapsed = 0.0;
    m_nevals = 0;
    m_age = 100000;
    m_atol = sqrt(std::numeric_limits<double>::epsilon());
    m_rtol = 1.0e-5;
}

void Jacobian::updateTransient(doublereal rdt, integer* mask)
{
    for (size_t n = 0; n < m_size; n++) {
        value(n,n) = m_ssdiag[n] - mask[n]*rdt;
    }
}

void Jacobian::incrementDiagonal(int j, doublereal d)
{
    m_ssdiag[j] += d;
    value(j,j) = m_ssdiag[j];
}

void Jacobian::eval(doublereal* x0, doublereal* resid0, doublereal rdt)
{
    m_nevals++;
    clock_t t0 = clock();
    //bfill(0.0);

    for (size_t n = 0; n < m_size; n++) {
        // perturb x(n); preserve sign(x(n))
        double xsave = x0[n];
        double dx;
        if (xsave >= 0) {
            dx = xsave*m_rtol + m_atol;
        } else {
            dx = xsave*m_rtol - m_atol;
        }
        x0[n] = xsave + dx;
        dx = x0[n] - xsave;
        double rdx = 1.0/dx;

        // calculate perturbed residual
        fill(m_r1.begin(), m_r1.end(), 0.0);
        m_residfunc->eval(0, x0, m_r1.data(), 0);

        // compute nth column of Jacobian
        for (size_t m = 0; m < m_size; m++) {
            value(m,n) = (m_r1[m] - resid0[m])*rdx;
        }

        x0[n] = xsave;
    }

    for (size_t n = 0; n < m_size; n++) {
        m_ssdiag[n] = value(n,n);
    }

    m_elapsed += double(clock() - t0)/CLOCKS_PER_SEC;
    m_age = 0;
}

} // namespace
