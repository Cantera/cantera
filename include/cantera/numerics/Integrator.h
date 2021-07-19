/**
 *  @file Integrator.h
 */

/**
 * @defgroup odeGroup ODE Integrators
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_INTEGRATOR_H
#define CT_INTEGRATOR_H
#include "FuncEval.h"
#include "PreconditionerBase.h"

#include "cantera/base/global.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

const int DIAG = 1;
const int DENSE = 2;
const int NOJAC = 4;
const int JAC = 8;
const int GMRES = 16;
const int BAND = 32;

/**
 * Specifies the method used to integrate the system of equations.
 * Not all methods are supported by all integrators.
 */
enum MethodType {
    BDF_Method, //!< Backward Differentiation
    Adams_Method //! Adams
};

//! Specifies the method used for iteration.
/*!
 * Not all methods are supported by all integrators.
 */
enum IterType {
    //!  Newton Iteration
    Newton_Iter,
    //! Functional Iteration
    Functional_Iter
};

//!  Abstract base class for ODE system integrators.
/*!
 *  @ingroup odeGroup
 */
class Integrator
{
public:
    //! Default Constructor
    Integrator() {
    }

    //!  Destructor
    virtual ~Integrator() {
    }

    //! Set error tolerances.
    /*!
     * @param reltol scalar relative tolerance
     * @param n      Number of equations
     * @param abstol array of N absolute tolerance values
     */
    virtual void setTolerances(doublereal reltol, size_t n,
                               doublereal* abstol) {
        warn("setTolerances");
    }

    //! Set error tolerances.
    /*!
     * @param reltol scalar relative tolerance
     * @param abstol scalar absolute tolerance
     */
    virtual void setTolerances(doublereal reltol, doublereal abstol) {
        warn("setTolerances");
    }

    //! Set the sensitivity error tolerances
    /*!
     * @param reltol scalar relative tolerance
     * @param abstol scalar absolute tolerance
     */
    virtual void setSensitivityTolerances(doublereal reltol, doublereal abstol)
    { }

    //! Set the problem type.
    /*!
     * @param probtype    Type of the problem
     *
     * @deprecated This funciton is to be removed along with the integer constants used
     * in conditionals to set the problem type currently. This includes DENSE, JAC,
     * NOJAC, BAND, and DIAG
     */
    virtual void setProblemType(int probtype) {
        warn_deprecated("Integrator::setProblemType()",
            "To be removed. Set linear solver type with setLinearSolverType");
        warn("setProblemType");
    }

    //! Set the linear solver type.
    /*!
     * @param linSolverType    Type of the linear solver
     */
    virtual void setLinearSolverType(std::string linSolverType) {
        warn("setLinearSolverType");
    }

    //! Set the preconditioner type.
    /*!
     * @param prectype    Type of the problem
     */
    virtual void setPreconditionerType(PreconditionerType prectype) {
        warn("setPreconditionerType");
    }

    //! Set preconditioner used by the linear solver
    /*!
     * @param preconditioner preconditioner object used for the linear solver
     */
    virtual void setPreconditioner(shared_ptr<PreconditionerBase> preconditioner) {
        m_preconditioner = preconditioner;
        m_prec_type = m_preconditioner->preconditionerType();
    }

    //! Solve a linear system Ax=b where A is the preconditioner
    /*!
     * @param[in] stateSize length of the rhs and output vectors
     * @param[in] rhs right hand side vector used in linear system
     * @param[out] output output vector for solution
     */
    virtual void preconditionerSolve(size_t stateSize, double* rhs, double* output) {
        m_preconditioner->solve(stateSize, rhs, output);
    }

    //! Return the preconditioner type
    virtual PreconditionerType preconditionerType() {
        return m_prec_type;
    }

    //! Return preconditioner reference to object
    virtual shared_ptr<PreconditionerBase> preconditioner() {
        return m_preconditioner;
    }

    //! Return the integrator problem type
    virtual std::string linearSolverType() {
        warn("linearSolverType");
        return "";
    }

    /**
     * Initialize the integrator for a new problem. Call after all options have
     * been set.
     * @param t0   initial time
     * @param func RHS evaluator object for system of equations.
     */
    virtual void initialize(doublereal t0, FuncEval& func) {
        warn("initialize");
    }

    virtual void reinitialize(doublereal t0, FuncEval& func) {
        warn("reinitialize");
    }

    //! Integrate the system of equations.
    /*!
     * @param tout Integrate to this time. Note that this is the
     *             absolute time value, not a time interval.
     */
    virtual void integrate(doublereal tout) {
        warn("integrate");
    }

    /**
     * Integrate the system of equations.
     * @param tout integrate to this time. Note that this is the
     * absolute time value, not a time interval.
     */
    virtual doublereal step(doublereal tout) {
        warn("step");
        return 0.0;
    }

    //! The current value of the solution of equation k.
    virtual doublereal& solution(size_t k) {
        warn("solution");
        return m_dummy;
    }

    //! The current value of the solution of the system of equations.
    virtual doublereal* solution() {
        warn("solution");
        return 0;
    }

    //! n-th derivative of the output function at time tout.
    virtual double* derivative(double tout, int n) {
        warn("derivative");
        return 0;
    }

    //! Order used during the last solution step
    virtual int lastOrder() const {
        warn("lastOrder");
        return 0;
    }

    //! The number of equations.
    virtual int nEquations() const {
        warn("nEquations");
        return 0;
    }

    //! The number of function evaluations.
    virtual int nEvals() const {
        warn("nEvals");
        return 0;
    }

    //! Set the maximum integration order that will be used.
    virtual void setMaxOrder(int n) {
        warn("setMaxorder");
    }

    //! Set the solution method
    virtual void setMethod(MethodType t) {
        warn("setMethodType");
    }

    //! Set the maximum step size
    virtual void setMaxStepSize(double hmax) {
        warn("setMaxStepSize");
    }

    //! Set the minimum step size
    virtual void setMinStepSize(double hmin) {
        warn("setMinStepSize");
    }

    //! Set the maximum permissible number of error test failures
    virtual void setMaxErrTestFails(int n) {
        warn("setMaxErrTestFails");
    }

    //! Set the maximum number of time-steps the integrator can take
    //!  before reaching the next output time
    //!  @param nmax The maximum number of steps, setting this value
    //!              to zero disables this option.
    virtual void setMaxSteps(int nmax) {
        warn("setMaxStep");
    }

    //! Returns the maximum number of time-steps the integrator can take
    //!  before reaching the next output time
    virtual int maxSteps() {
        warn("maxSteps");
        return 0;
    }

    virtual void setBandwidth(int N_Upper, int N_Lower) {
        warn("setBandwidth");
    }

    virtual int nSensParams() {
        warn("nSensParams()");
        return 0;
    }

    virtual double sensitivity(size_t k, size_t p) {
        warn("sensitivity");
        return 0.0;
    }

    //! Get nonlinear solver stats from integrator
    virtual AnyMap nonlinearSolverStats() const {
        AnyMap stats;
        warn("nonlinearSolverStats");
        return stats;
    }

    //! Get linear solver stats from integrator
    virtual AnyMap linearSolverStats() const {
        AnyMap stats;
        warn("linearSolverStats");
        return stats;
    }

protected:
    //! Pointer to preconditioner object used in integration which is
    //! set by setPreconditioner and initialized inside of
    //! ReactorNet::initialize()
    shared_ptr<PreconditionerBase> m_preconditioner;
    //! Type of preconditioning used in applyOptions
    PreconditionerType m_prec_type = PreconditionerType::NO_PRECONDITION;

private:
    doublereal m_dummy;
    void warn(const std::string& msg) const {
        writelog(">>>> Warning: method "+msg+" of base class "
                 +"Integrator called. Nothing done.\n");
    }
};

// defined in ODE_integrators.cpp
Integrator* newIntegrator(const std::string& itype);

} // namespace

#endif
