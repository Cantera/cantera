//! @file ResidEval.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_RESIDEVAL_H
#define CT_RESIDEVAL_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/utilities.h"

#include <cstdio>

namespace Cantera
{

const int c_NONE = 0;
const int c_GE_ZERO = 1;
const int c_GT_ZERO = 2;
const int c_LE_ZERO = -1;
const int c_LT_ZERO = -2;

/**
 *  Virtual base class for DAE residual function evaluators.
 *  Classes derived from ResidEval evaluate the residual function
 * \f[
 *             \vec{F}(t,\vec{y}, \vec{y^\prime})
 * \f]
 * The DAE solver attempts to find a solution y(t) such that F = 0.
 *  @ingroup DAE_Group
 */
class ResidEval
{
public:
    ResidEval() {}
    virtual ~ResidEval() {}

    /**
     * Constrain solution component k. Possible values for
     * 'flag' are:
     *   - c_NONE       no constraint
     *   - c_GE_ZERO    >= 0
     *   - c_GT_ZERO    >  0
     *   - c_LE_ZERO    <= 0
     *   - c_LT_ZERO    <  0
     */
    virtual void constrain(const int k, const int flag) {
        m_constrain[k] = flag;
    }
    int constraint(const int k) const {
        return getValue(m_constrain, k, c_NONE);
    }

    //! Initialization function
    virtual void initSizes() {
        int neq = nEquations();
        m_alg.resize(neq, 0);
    }

    /**
     * Specify that solution component k is purely algebraic - that is, the
     * derivative of this component does not appear in the residual function.
     */
    virtual void setAlgebraic(const int k) {
        if ((int) m_alg.size() < (k+1)) {
            initSizes();
        }
        m_alg[k] = 1;
    }

    virtual bool isAlgebraic(const int k) {
        return (m_alg[k] == 1);
    }

    /**
     * Evaluate the residual function. Called by the integrator.
     * @param t time. (input)
     * @param y solution vector. (input)
     * @param ydot rate of change of solution vector. (input)
     * @param r residual vector (output)
     */
    virtual int eval(const doublereal t, const doublereal* const y,
                     const doublereal* const ydot,
                     doublereal* const r) {
        throw NotImplementedError("ResidEval::eval");
    }

    virtual int evalSS(const doublereal t, const doublereal* const y,
                       doublereal* const r) {
        return eval(t, y, 0, r);
    }

    virtual int evalSimpleTD(const doublereal t, const doublereal* const y,
                             const doublereal* const yold, doublereal deltaT,
                             doublereal* const r) {
        int nn = nEquations();
        vector_fp ydot(nn);
        for (int i = 0; i < nn; i++) {
            ydot[i] = (y[i] - yold[i]) / deltaT;
        }
        return eval(t, y, ydot.data(), r);
    }

    //! Fill in the initial conditions
    /*!
     * Values for both the solution and the value of ydot may be provided.
     *
     * @param[in] t0     Time
     * @param[out] y     Solution vector
     * @param[out] ydot  Rate of change of solution vector.
     *
     * @returns a flag to indicate that operation is successful.
     *            1  Means a successful operation
     *           -0 or neg value Means an unsuccessful operation
     */
    virtual int getInitialConditions(const doublereal t0, doublereal* const y,
                                     doublereal* const ydot) {
        initSizes();
        throw NotImplementedError("ResidEval::GetInitialConditions");
        return 1;
    }

    //! Return the number of equations in the equation system
    virtual int nEquations() const = 0;

    //! Write out to a file or to standard output the current solution
    /*!
     * ievent is a description of the event that caused this function to be
     * called.
     */
    virtual void writeSolution(int ievent, const double time,
                               const double deltaT,
                               const int time_step_num,
                               const double* y, const double* ydot) {
        int k;
        writelog("ResidEval::writeSolution\n");
        writelogf("     Time = %g, ievent = %d, deltaT = %g\n", time, ievent, deltaT);
        if (ydot) {
            writelogf(" k    y[]  ydot[]\n");
            for (k = 0; k < nEquations(); k++) {
                writelogf("%d %g %g\n", k, y[k], ydot[k]);
            }
        } else {
            writelogf(" k    y[]\n");
            for (k = 0; k < nEquations(); k++) {
                writelogf("%d %g \n", k, y[k]);
            }
        }
    }

    //! Return the number of parameters in the calculation
    /*!
     * This is the number of parameters in the sensitivity calculation. We have
     * set this to zero and have included it for later expansion
     */
    int nparams() const {
        return 0;
    }

protected:
    //! Mapping vector that stores whether a degree of freedom is a DAE or not
    /*!
     * The first index is the equation number. The second index is 1 if it is a
     * DAE, and zero if it is not.
     */
    vector_int m_alg;
    std::map<int, int> m_constrain;
};

}

#endif
