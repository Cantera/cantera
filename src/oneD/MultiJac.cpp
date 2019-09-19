//! @file MultiJac.cpp Implementation file for class MultiJac

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/MultiJac.h"
#include <ctime>

using namespace std;

namespace Cantera
{

MultiJac::MultiJac(OneDim& r)
    : BandMatrix(r.size(),r.bandwidth(),r.bandwidth())
{
    m_size = r.size();
    m_points = r.points();
    m_resid = &r;
    m_r1.resize(m_size);
    m_ssdiag.resize(m_size);
    m_mask.resize(m_size);
    m_elapsed = 0.0;
    m_nevals = 0;
    m_age = 100000;
    m_atol = sqrt(std::numeric_limits<double>::epsilon());
    m_rtol = 1.0e-5;
}

void MultiJac::updateTransient(doublereal rdt, integer* mask)
{
    for (size_t n = 0; n < m_size; n++) {
        value(n,n) = m_ssdiag[n] - mask[n]*rdt;
    }
}

void MultiJac::incrementDiagonal(int j, doublereal d)
{
    m_ssdiag[j] += d;
    value(j,j) = m_ssdiag[j];
}

void MultiJac::eval(doublereal* x0, doublereal* resid0, doublereal rdt)
{
    m_nevals++;
    clock_t t0 = clock();
    bfill(0.0);
    size_t ipt=0;

    for (size_t j = 0; j < m_points; j++) {
        size_t nv = m_resid->nVars(j);
        for (size_t n = 0; n < nv; n++) {
            // perturb x(n); preserve sign(x(n))
            double xsave = x0[ipt];
            double dx;
            if (xsave >= 0) {
                dx = xsave*m_rtol + m_atol;
            } else {
                dx = xsave*m_rtol - m_atol;
            }
            x0[ipt] = xsave + dx;
            dx = x0[ipt] - xsave;
            double rdx = 1.0/dx;

            // calculate perturbed residual
            m_resid->eval(j, x0, m_r1.data(), rdt, 0);

            // compute nth column of Jacobian
            for (size_t i = j - 1; i != j+2; i++) {
                if (i != npos && i < m_points) {
                    size_t mv = m_resid->nVars(i);
                    size_t iloc = m_resid->loc(i);
                    for (size_t m = 0; m < mv; m++) {
                        value(m+iloc,ipt) = (m_r1[m+iloc] - resid0[m+iloc])*rdx;
                    }
                }
            }
            x0[ipt] = xsave;
            ipt++;
        }
    }

    for (size_t n = 0; n < m_size; n++) {
        m_ssdiag[n] = value(n,n);
    }

    m_elapsed += double(clock() - t0)/CLOCKS_PER_SEC;
    m_age = 0;
}

} // namespace
