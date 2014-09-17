/**
 * @file L_matrix.h
 * Functions to evaluate portions of the L matrix needed for
 * multicomponent transport properties.
 */

#ifndef CT_LMATRIX_H
#define CT_LMATRIX_H

#include "cantera/transport/MultiTransport.h"

namespace Cantera
{

//! Constant to compare dimensionless heat capacities against zero
const doublereal Min_C_Internal = 0.001;

bool MultiTransport::hasInternalModes(size_t j)
{
    return (m_cinternal[j] > Min_C_Internal);
}

void MultiTransport::eval_L0000(const doublereal* const x)
{
    doublereal prefactor = 16.0*m_temp/25.0;
    doublereal sum;
    for (size_t i = 0; i < m_nsp; i++)  {
        //  subtract-off the k=i term to account for the first delta
        //  function in Eq. (12.121)

        sum = -x[i]/m_bdiff(i,i);
        for (size_t k = 0; k < m_nsp; k++) {
            sum += x[k]/m_bdiff(i,k);
        }

        sum /= m_mw[i];
        for (size_t j = 0; j != m_nsp; ++j) {
            m_Lmatrix(i,j) = prefactor * x[j]
                             * (m_mw[j] * sum + x[i]/m_bdiff(i,j));
        }
        // diagonal term is zero
        m_Lmatrix(i,i) = 0.0;
    }
}

void MultiTransport::eval_L0010(const doublereal* const x)
{
    doublereal prefactor = 1.6*m_temp;

    doublereal sum, wj, xj;
    for (size_t j = 0; j < m_nsp; j++) {
        //constant = prefactor * x[j];
        xj = x[j];
        wj = m_mw[j];
        sum = 0.0;
        for (size_t i = 0; i < m_nsp; i++) {
            m_Lmatrix(i,j + m_nsp) = - prefactor * x[i] * xj * m_mw[i] *
                                     (1.2 * m_cstar(j,i) - 1.0) /
                                     ((wj + m_mw[i]) * m_bdiff(j,i));

            //  the next term is independent of "j";
            //  need to do it for the "j,j" term
            sum -= m_Lmatrix(i,j+m_nsp);
        }
        m_Lmatrix(j,j+m_nsp) += sum;
    }
}

void MultiTransport::eval_L1000()
{
    for (size_t j = 0; j < m_nsp; j++) {
        for (size_t i = 0; i < m_nsp; i++) {
            m_Lmatrix(i+m_nsp,j) = m_Lmatrix(j,i+m_nsp);
        }
    }
}

void MultiTransport::eval_L1010(const doublereal* x)
{
    const doublereal fiveover3pi = 5.0/(3.0*Pi);
    doublereal prefactor = (16.0*m_temp)/25.0;

    doublereal constant1, wjsq, constant2, constant3, constant4,
               fourmj, threemjsq, sum, sumwij;;
    doublereal term1, term2;

    for (size_t j = 0; j < m_nsp; j++) {

        // get constant terms that depend on just species "j"

        constant1 =   prefactor*x[j];
        wjsq      =   m_mw[j]*m_mw[j];
        constant2 =   13.75*wjsq;
        constant3 =   m_crot[j]/m_rotrelax[j];
        constant4 =   7.5*wjsq;
        fourmj    =   4.0*m_mw[j];
        threemjsq =   3.0*m_mw[j]*m_mw[j];
        sum =         0.0;
        for (size_t i = 0; i < m_nsp; i++) {

            sumwij = m_mw[i] + m_mw[j];
            term1 = m_bdiff(i,j) * sumwij*sumwij;
            term2 = fourmj*m_astar(i,j)*(1.0 + fiveover3pi*
                                         (constant3 +
                                          (m_crot[i]/m_rotrelax[i])));   //  see Eq. (12.125)

            m_Lmatrix(i+m_nsp,j+m_nsp) = constant1*x[i]*m_mw[i] /(m_mw[j]*term1) *
                                         (constant2 - threemjsq*m_bstar(i,j)
                                          - term2*m_mw[j]);

            sum += x[i] /(term1) *
                   (constant4 + m_mw[i]*m_mw[i]*
                    (6.25 - 3.0*m_bstar(i,j)) + term2*m_mw[i]);
        }

        m_Lmatrix(j+m_nsp,j+m_nsp) -= sum*constant1;
    }
}

void MultiTransport::eval_L1001(const doublereal* x)
{
    doublereal prefactor = 32.00*m_temp/(5.00*Pi);
    doublereal constant, sum;
    size_t n2 = 2*m_nsp;
    int npoly = 0;
    for (size_t j = 0; j < m_nsp; j++) {
        //        collect terms that depend only on "j"
        if (hasInternalModes(j)) {
            constant = prefactor*m_mw[j]*x[j]*m_crot[j]/(m_cinternal[j]*m_rotrelax[j]);
            sum = 0.0;
            for (size_t i = 0; i < m_nsp; i++) {
                //           see Eq. (12.127)
                m_Lmatrix(i+m_nsp,j+n2) = constant * m_astar(j,i) * x[i] /
                                          ((m_mw[j] + m_mw[i]) * m_bdiff(j,i));
                sum += m_Lmatrix(i+m_nsp,j+n2);
            }
            npoly++;
            m_Lmatrix(j+m_nsp,j+n2) += sum;
        } else {
            for (size_t i = 0; i < m_nsp; i++) {
                m_Lmatrix(i+m_nsp,j+n2) = 0.0;
            }
        }
    }
}

void MultiTransport::eval_L0001()
{
    size_t n2 = 2*m_nsp;
    for (size_t j = 0; j < m_nsp; j++) {
        for (size_t i = 0; i < m_nsp; i++) {
            m_Lmatrix(i,j+n2) = 0.0;
        }
    }
}

void MultiTransport::eval_L0100()
{
    size_t n2 = 2*m_nsp;
    for (size_t j = 0; j < m_nsp; j++)
        for (size_t i = 0; i < m_nsp; i++) {
            m_Lmatrix(i+n2,j) = 0.0;    //  see Eq. (12.123)
        }
}

void MultiTransport::eval_L0110()
{
    size_t n2 = 2*m_nsp;
    for (size_t j = 0; j < m_nsp; j++)
        for (size_t i = 0; i < m_nsp; i++) {
            m_Lmatrix(i+n2,j+m_nsp) = m_Lmatrix(j+m_nsp,i+n2);    //  see Eq. (12.123)
        }
}

void MultiTransport::eval_L0101(const doublereal* x)
{
    const doublereal fivepi = 5.00*Pi;
    const doublereal eightoverpi = 8.0 / Pi;

    doublereal prefactor = 4.00*m_temp;
    size_t n2 = 2*m_nsp;
    doublereal constant1, constant2, diff_int, sum;
    for (size_t i = 0; i < m_nsp; i++) {
        if (hasInternalModes(i)) {
            //        collect terms that depend only on "i"
            constant1 = prefactor*x[i]/m_cinternal[i];
            constant2 = 12.00*m_mw[i]*m_crot[i] /
                        (fivepi*m_cinternal[i]*m_rotrelax[i]);
            sum = 0.0;
            for (size_t k = 0; k < m_nsp; k++) {
                //           see Eq. (12.131)
                diff_int = m_bdiff(i,k);
                m_Lmatrix(k+n2,i+n2) = 0.0;
                sum += x[k]/diff_int;
                if (k != i) sum += x[k]*m_astar(i,k)*constant2 /
                                       (m_mw[k]*diff_int);
            }
            //        see Eq. (12.130)
            m_Lmatrix(i+n2,i+n2) =
                - eightoverpi*m_mw[i]*x[i]*x[i]*m_crot[i] /
                (m_cinternal[i]*m_cinternal[i]*GasConstant*m_visc[i]*m_rotrelax[i])
                - constant1*sum;
        } else {
            for (size_t k = 0; k < m_nsp; k++) {
                m_Lmatrix(i+n2,i+n2) = 1.0;
            }
        }
    }
}
}

#endif
