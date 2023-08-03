//! @file ResidJacEval.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#define CT_SKIP_DEPRECATION_WARNINGS
#include "cantera/numerics/ResidJacEval.h"
#include "cantera/base/global.h"

namespace Cantera
{
ResidJacEval::ResidJacEval(double atol) :
    m_atol(atol)
{
    warn_deprecated("class ResidJacEval", "To be removed after Cantera 3.0");
}

int ResidJacEval::nEquations() const
{
    return neq_;
}

void ResidJacEval::setAtol(double atol)
{
    m_atol = atol;
    if (m_atol <= 0.0) {
        throw CanteraError("ResidJacEval::setAtol",
                           "atol must be greater than zero");
    }
}

int ResidJacEval::getInitialConditions(double t0, double* const y,
                                       double* const ydot)
{
    for (int i = 0; i < neq_; i++) {
        y[i] = 0.0;
    }
    if (ydot) {
        for (int i = 0; i < neq_; i++) {
            ydot[i] = 0.0;
        }
    }
    return 1;
}

void ResidJacEval::user_out2(const int ifunc, const double t,
                             const double deltaT, const double* y,
                             const double* ydot)
{
}

void ResidJacEval::user_out(const int ifunc, const double t,
                            const double* y, const double* ydot)
{
    user_out2(ifunc, t, 0.0, y, ydot);
}

int ResidJacEval::evalTimeTrackingEqns(const double t,
                                       const double delta_t,
                                       const double* y,
                                       const double* ydot)
{
    return 1;
}

int ResidJacEval::calcDeltaSolnVariables(const double t,
                                         const double* const ySoln,
                                         const double* const ySolnDot,
                                         double* const deltaYSoln,
                                         const double* const solnWeights)
{
    if (!solnWeights) {
        for (int i = 0; i < neq_; i++) {
            deltaYSoln[i] = m_atol + fabs(1.0E-6 * ySoln[i]);
        }
    } else {
        for (int i = 0; i < neq_; i++) {
            deltaYSoln[i] = std::max(1.0E-2 * solnWeights[i], 1.0E-6 * fabs(ySoln[i]));
        }
    }
    return 1;
}

void ResidJacEval::calcSolnScales(const double t,
                                  const double* const ysoln,
                                  const double* const ysolnOld,
                                  double* const ysolnScales)
{
    if (ysolnScales && ysolnScales[0] == 0.0) {
        for (int i = 0; i < neq_; i++) {
            ysolnScales[i] = 1.0;
        }
    }
}

double ResidJacEval::filterNewStep(double t, const double* const ybase, double* const step)
{
    return 0.0;
}

double ResidJacEval::filterSolnPrediction(double t, double* const y)
{
    return 0.0;
}

bool ResidJacEval::evalStoppingCritera(const double t,
                                       const double delta_t,
                                       const double* const y,
                                       const double* const ydot)
{
    return false;
}

int ResidJacEval::matrixConditioning(double* const matrix, const int nrows,
                                     double* const rhs)
{
    return 1;
}

int ResidJacEval::evalResidNJ(const double t, const double deltaT,
                              const double* y, const double* ydot,
                              double* const resid,
                              const ResidEval_Type_Enum evalType,
                              const int id_x, const double delta_x)
{
    throw NotImplementedError("ResidJacEval::evalResidNJ");
}

int ResidJacEval::eval(const double t, const double* const y, const double* const ydot,
                       double* const r)
{
    double deltaT = -1.0;
    return evalResidNJ(t, deltaT, y, ydot, r);
}

int ResidJacEval::evalJacobian(const double t, const double delta_t,
                               double cj, const double* const y,
                               const double* const ydot, DenseMatrix& J,
                               double* const resid)
{
    double* const* jac_colPts = J.colPts();
    return evalJacobianDP(t, delta_t, cj, y, ydot, jac_colPts, resid);
}

int ResidJacEval::evalJacobianDP(const double t, const double delta_t,
                                 const double c_j,
                                 const double* const y,
                                 const double* const ydot,
                                 double* const* jac_colPts,
                                 double* const resid)
{
    throw NotImplementedError("ResidJacEval::evalJacobianDP");
}

}
