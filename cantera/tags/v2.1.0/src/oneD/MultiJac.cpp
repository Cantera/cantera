/**
 *  @file MultiJac.cpp Implementation file for class MultiJac
 */

/*
 *  Copyright 2002 California Institute of Technology
 */

#include "cantera/oneD/MultiJac.h"
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
    doublereal ff = 1.0;
    while (1.0 + ff != 1.0) {
        ff *= 0.5;
    }
    m_atol = sqrt(ff);
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
    size_t n, m, ipt=0, j, nv, mv, iloc;
    doublereal rdx, dx, xsave;

    for (j = 0; j < m_points; j++) {
        nv = m_resid->nVars(j);
        for (n = 0; n < nv; n++) {

            // perturb x(n)
            xsave = x0[ipt];
            dx = m_atol + fabs(xsave)*m_rtol;
            x0[ipt] = xsave + dx;
            dx = x0[ipt] - xsave;
            rdx = 1.0/dx;

            // calculate perturbed residual
            m_resid->eval(j, x0, DATA_PTR(m_r1), rdt, 0);

            // compute nth column of Jacobian
            for (size_t i = j - 1; i != j+2; i++) {
                if (i != npos && i < m_points) {
                    mv = m_resid->nVars(i);
                    iloc = m_resid->loc(i);
                    for (m = 0; m < mv; m++) {
                        value(m+iloc,ipt) = (m_r1[m+iloc]
                                             - resid0[m+iloc])*rdx;
                    }
                }
            }
            x0[ipt] = xsave;
            ipt++;
        }

    }

    for (n = 0; n < m_size; n++) {
        m_ssdiag[n] = value(n,n);
    }

    m_elapsed += double(clock() - t0)/CLOCKS_PER_SEC;
    m_age = 0;
}

} // namespace
