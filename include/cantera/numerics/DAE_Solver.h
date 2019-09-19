/**
 *  @file DAE_Solver.h
 *  Header file for class DAE_Solver
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_DAE_Solver_H
#define CT_DAE_Solver_H

#include "ResidJacEval.h"
#include "cantera/base/global.h"

namespace Cantera
{

class Jacobian
{
public:
    Jacobian() {}
    virtual ~Jacobian() {}
    virtual bool supplied() {
        return false;
    }
    virtual bool isBanded() {
        return false;
    }
    virtual int lowerBandWidth() {
        return 0;
    }
    virtual int upperBandWidth() {
        return 0;
    }
};

class BandedJacobian : public Jacobian
{
public:
    BandedJacobian(int ml, int mu) {
        m_ml = ml;
        m_mu = mu;
    }
    virtual bool supplied() {
        return false;
    }
    virtual bool isBanded() {
        return true;
    }
    virtual int lowerBandWidth() {
        return m_ml;
    }
    virtual int upperBandWidth() {
        return m_mu;
    }
protected:
    int m_ml, m_mu;
};

const int cDirect = 0;
const int cKrylov = 1;


/**
 * Wrapper for DAE solvers
 *
 * @attention This class currently does not have any test cases or examples. Its
 *     implementation may be incomplete, and future changes to Cantera may
 *     unexpectedly cause this class to stop working. If you use this class,
 *     please consider contributing examples or test cases. In the absence of
 *     new tests or examples, this class may be deprecated and removed in a
 *     future version of Cantera. See
 *     https://github.com/Cantera/cantera/issues/267 for additional information.
 */
class DAE_Solver
{
public:
    DAE_Solver(ResidJacEval& f) :
        m_resid(f),
        m_neq(f.nEquations()),
        m_time(0.0) {
    }

    virtual ~DAE_Solver() {}

    /**
     * Set error tolerances. This version specifies a scalar
     * relative tolerance, and a vector absolute tolerance.
     */
    virtual void setTolerances(doublereal reltol,
                               doublereal* abstol) {
        warn("setTolerances");
    }

    /**
     * Set error tolerances. This version specifies a scalar
     * relative tolerance, and a scalar absolute tolerance.
     */
    virtual void setTolerances(doublereal reltol, doublereal abstol) {
        warn("setTolerances");
    }

    /**
     * Specify a Jacobian evaluator. If this method is not called,
     * the Jacobian will be computed by finite difference.
     */
    void setJacobian(Jacobian& jac) {
        warn("setJacobian");
    }

    virtual void setLinearSolverType(int solverType) {
        warn("setLinearSolverType");
    }

    virtual void setDenseLinearSolver() {
        warn("setDenseLinearSolver");
    }

    virtual void setBandedLinearSolver(int m_upper, int m_lower) {
        warn("setBandedLinearSolver");
    }
    virtual void setMaxStepSize(doublereal dtmax) {
        warn("setMaxStepSize");
    }
    virtual void setMaxOrder(int n) {
        warn("setMaxOrder");
    }
    virtual void setMaxNumSteps(int n) {
        warn("setMaxNumSteps");
    }
    virtual void setInitialStepSize(doublereal h0) {
        warn("setInitialStepSize");
    }
    virtual void setStopTime(doublereal tstop) {
        warn("setStopTime");
    }
    virtual void setMaxErrTestFailures(int n) {
        warn("setMaxErrTestFailures");
    }
    virtual void setMaxNonlinIterations(int n) {
        warn("setMaxNonlinIterations");
    }
    virtual void setMaxNonlinConvFailures(int n) {
        warn("setMaxNonlinConvFailures");
    }
    virtual void inclAlgebraicInErrorTest(bool yesno) {
        warn("inclAlgebraicInErrorTest");
    }

    //! Calculate consistent value of the starting solution given the starting
    //! solution derivatives
    /**
     * This method may be called if the initial conditions do not satisfy the
     * residual equation F = 0. Given the derivatives of all variables, this
     * method computes the initial y values.
     */
    virtual void correctInitial_Y_given_Yp(doublereal* y, doublereal* yp,
                                           doublereal tout) {
        warn("correctInitial_Y_given_Yp");
    }

    //! Calculate consistent value of the algebraic constraints and derivatives
    //! at the start of the problem
    /**
     * This method may be called if the initial conditions do not satisfy the
     * residual equation F = 0. Given the initial values of all differential
     * variables, it computes the initial values of all algebraic variables and
     * the initial derivatives of all differential variables.
     *  @param y      Calculated value of the solution vector after the procedure ends
     *  @param yp     Calculated value of the solution derivative after the procedure
     *  @param tout   The first value of t at which a soluton will be
     *                requested (from IDASolve).  (This is needed here to
     *                determine the direction of integration and rough scale
     *                in the independent variable t.
     */
    virtual void correctInitial_YaYp_given_Yd(doublereal* y, doublereal* yp,
            doublereal tout) {
        warn("correctInitial_YaYp_given_Yd");
    }

    /**
     * Solve the system of equations up to time tout.
     */
    virtual int solve(doublereal tout) {
        warn("solve");
        return 0;
    }

    /**
     * Take one internal step.
     */
    virtual doublereal step(doublereal tout) {
        warn("step");
        return 0;
    }

    /// Number of equations.
    int nEquations() const {
        return m_resid.nEquations();
    }

    /**
     * initialize. Base class method does nothing.
     */
    virtual void init(doublereal t0) {}

    /**
     * Set a solver-specific input parameter.
     */
    virtual void setInputParameter(int flag, doublereal value) {
        warn("setInputParameter");
    }

    /**
     * Get the value of a solver-specific output parameter.
     */
    virtual doublereal getOutputParameter(int flag) const {
        warn("getOutputParameter");
        return 0.0;
    }

    /// the current value of solution component k.
    virtual doublereal solution(int k) const {
        warn("solution");
        return 0.0;
    }

    virtual const doublereal* solutionVector() const {
        warn("solutionVector");
        return &m_dummy;
    }

    /// the current value of the derivative of solution component k.
    virtual doublereal derivative(int k) const {
        warn("derivative");
        return 0.0;
    }

    virtual const doublereal* derivativeVector() const {
        warn("derivativeVector");
        return &m_dummy;
    }

protected:
    doublereal m_dummy;
    ResidJacEval& m_resid;

    //! Number of total equations in the system
    integer m_neq;
    doublereal m_time;

private:
    void warn(const std::string& msg) const {
        writelog(">>>> Warning: method "+msg+" of base class "
                 +"DAE_Solver called. Nothing done.\n");
    }
};


//! Factor method for choosing a DAE solver
/*!
 * @param itype  String identifying the type (IDA is the only option)
 * @param f      Residual function to be solved by the DAE algorithm
 * @returns a point to the instantiated DAE_Solver object
 */
DAE_Solver* newDAE_Solver(const std::string& itype, ResidJacEval& f);

}

#endif
