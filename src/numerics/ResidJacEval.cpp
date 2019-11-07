//! @file ResidJacEval.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ct_defs.h"
#include "cantera/numerics/ResidJacEval.h"

using namespace std;

namespace Cantera
{
ResidJacEval::ResidJacEval(doublereal atol) :
    m_atol(atol)
{
}

int ResidJacEval::nEquations() const
{
    return neq_;
}

void ResidJacEval::setAtol(doublereal atol)
{
    m_atol = atol;
    if (m_atol <= 0.0) {
        throw CanteraError("ResidJacEval::setAtol",
                           "atol must be greater than zero");
    }
}

int ResidJacEval::getInitialConditions(doublereal t0, doublereal* const y,
                                       doublereal* const ydot)
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

void ResidJacEval::user_out2(const int ifunc, const doublereal t,
                             const doublereal deltaT, const doublereal* y,
                             const doublereal* ydot)
{
}

void ResidJacEval::user_out(const int ifunc, const doublereal t,
                            const doublereal* y, const doublereal* ydot)
{
    user_out2(ifunc, t, 0.0, y, ydot);
}

int ResidJacEval::evalTimeTrackingEqns(const doublereal t,
                                       const doublereal delta_t,
                                       const doublereal* y,
                                       const doublereal* ydot)
{
    return 1;
}

int ResidJacEval::calcDeltaSolnVariables(const doublereal t,
                                         const doublereal* const ySoln,
                                         const doublereal* const ySolnDot,
                                         doublereal* const deltaYSoln,
                                         const doublereal* const solnWeights)
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

void ResidJacEval::calcSolnScales(const doublereal t,
                                  const doublereal* const ysoln,
                                  const doublereal* const ysolnOld,
                                  doublereal* const ysolnScales)
{
    if (ysolnScales && ysolnScales[0] == 0.0) {
        for (int i = 0; i < neq_; i++) {
            ysolnScales[i] = 1.0;
        }
    }
}

doublereal ResidJacEval::filterNewStep(doublereal t, const doublereal* const ybase, doublereal* const step)
{
    return 0.0;
}

doublereal ResidJacEval::filterSolnPrediction(doublereal t, doublereal* const y)
{
    return 0.0;
}

bool ResidJacEval::evalStoppingCritera(const doublereal t,
                                       const doublereal delta_t,
                                       const doublereal* const y,
                                       const doublereal* const ydot)
{
    return false;
}

int ResidJacEval::matrixConditioning(doublereal* const matrix, const int nrows,
                                     doublereal* const rhs)
{
    return 1;
}

int ResidJacEval::evalResidNJ(const doublereal t, const doublereal deltaT,
                              const doublereal* y, const doublereal* ydot,
                              doublereal* const resid,
                              const ResidEval_Type_Enum evalType,
                              const int id_x, const doublereal delta_x)
{
    throw NotImplementedError("ResidJacEval::evalResidNJ");
}

int ResidJacEval::eval(const doublereal t, const doublereal* const y, const doublereal* const ydot,
                       doublereal* const r)
{
    double deltaT = -1.0;
    return evalResidNJ(t, deltaT, y, ydot, r);
}

int ResidJacEval::evalJacobian(const doublereal t, const doublereal delta_t,
                               doublereal cj, const doublereal* const y,
                               const doublereal* const ydot, DenseMatrix& J,
                               doublereal* const resid)
{
    doublereal* const* jac_colPts = J.colPts();
    return evalJacobianDP(t, delta_t, cj, y, ydot, jac_colPts, resid);
}

int ResidJacEval::evalJacobianDP(const doublereal t, const doublereal delta_t,
                                 const doublereal c_j,
                                 const doublereal* const y,
                                 const doublereal* const ydot,
                                 doublereal* const* jac_colPts,
                                 doublereal* const resid)
{
    throw NotImplementedError("ResidJacEval::evalJacobianDP");
}

}
