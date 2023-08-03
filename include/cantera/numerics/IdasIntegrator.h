/**
 *  @file IdasIntegrator.h
 *  Header file for class IdasIntegrator
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_IDASINTEGRATOR_H
#define CT_IDASINTEGRATOR_H

#include "cantera/numerics/Integrator.h"
#include "cantera/base/ctexceptions.h"
#include "sundials/sundials_nvector.h"
#include "cantera/numerics/SundialsContext.h"

namespace Cantera
{

/**
 * Wrapper for Sundials IDAS solver
 * @see FuncEval.h. Classes that use IdasIntegrator:
 * FlowReactor
 */
class IdasIntegrator : public Integrator
{
public:
    /**
     *  Constructor. Default settings: dense Jacobian, no user-supplied
     *  Jacobian function, Newton iteration.
     */
    IdasIntegrator();
    ~IdasIntegrator() override;
    void setTolerances(double reltol, size_t n, double* abstol) override;
    void setTolerances(double reltol, double abstol) override;
    void setSensitivityTolerances(double reltol, double abstol) override;
    void setLinearSolverType(const string& linearSolverType) override;
    void initialize(double t0, FuncEval& func) override;
    void reinitialize(double t0, FuncEval& func) override;
    void integrate(double tout) override;
    double step(double tout) override;
    double& solution(size_t k) override;
    double* solution() override;
    int nEquations() const override {
        return static_cast<int>(m_neq);
    }
    int maxOrder() const override {
        return m_maxord;
    }
    void setMaxOrder(int n) override;
    void setMaxStepSize(double hmax) override;
    void setMaxSteps(int nmax) override;
    int maxSteps() override;
    void setMaxErrTestFails(int n) override;
    AnyMap solverStats() const override;
    int nSensParams() override {
        return static_cast<int>(m_np);
    }
    double sensitivity(size_t k, size_t p) override;

    //! Returns a string listing the weighted error estimates associated
    //! with each solution component.
    //! This information can be used to identify which variables are
    //! responsible for integrator failures or unexpected small timesteps.
    string getErrorInfo(int N);

    //! Error message information provide by IDAS
    string m_error_message;

    int maxNonlinIterations() const override {
        return m_maxNonlinIters;
    }
    void setMaxNonlinIterations(int n) override;

    int maxNonlinConvFailures() const override {
        return m_maxNonlinConvFails;
    }
    void setMaxNonlinConvFailures(int n) override;

    bool algebraicInErrorTest() const override {
        return !m_setSuppressAlg;
    }
    void includeAlgebraicInErrorTest(bool yesno) override;

    void setMethod(MethodType t) override;

protected:
    protected:
    //! Applies user-specified options to the underlying IDAS solver. Called
    //! during integrator initialization or reinitialization.
    void applyOptions();

private:
    void sensInit(double t0, FuncEval& func);

    //! Check whether an IDAS method indicated an error. If so, throw an exception
    //! containing the method name and the error code stashed by the ida_err() function.
    void checkError(long flag, const string& ctMethod, const string& idaMethod) const;

    size_t m_neq; //!< Number of equations/variables in the system
    void* m_ida_mem = nullptr; //!< Pointer to the IDA memory for the problem
    void* m_linsol = nullptr; //!< Sundials linear solver object
    void* m_linsol_matrix = nullptr; //!< matrix used by Sundials
    SundialsContext m_sundials_ctx; //!< SUNContext object for Sundials>=6.0

    //! Object implementing the DAE residual function @f$ f(t, y, \dot{y}) = 0 @f$
    FuncEval* m_func = nullptr;

    double m_t0 = 0.0; //!< The start time for the integrator
    double m_time; //!< The current integrator time, corresponding to #m_y

    //! The latest time reached by the integrator. May be greater than #m_time
    double m_tInteg;

    N_Vector m_y = nullptr; //!< The current system state
    N_Vector m_ydot = nullptr; //!< The time derivatives of the system state
    N_Vector m_abstol = nullptr; //!< Absolute tolerances for each state variable
    string m_type = "DENSE"; //!< The linear solver type. @see setLinearSolverType()

    //! Flag indicating whether scalar (`IDA_SS`) or vector (`IDA_SV`) absolute
    //! tolerances are being used.
    int m_itol;

    int m_maxord = 0; //!< Maximum order allowed for the BDF method
    double m_reltol = 1.0e-9; //!< Relative tolerance for all state variables
    double m_abstols = 1.0e-15; //!< Scalar absolute tolerance
    double m_reltolsens; //!< Scalar relative tolerance for sensitivities
    double m_abstolsens; //!< Scalar absolute tolerance for sensitivities

    //!! Number of variables for which absolute tolerances were provided
    size_t m_nabs = 0;

    double m_hmax = 0.0; //!< Maximum time step size. Zero means infinity.

    //! Maximum number of internal steps taken in a call to integrate()
    int m_maxsteps = 20000;

    //! Maximum number of error test failures in attempting one step
    int m_maxErrTestFails = -1;

    size_t m_np; //!< Number of sensitivity parameters
    N_Vector* m_yS = nullptr; //!< Sensitivities of y, size #m_np by #m_neq.
    N_Vector* m_ySdot = nullptr; //!< Sensitivities of ydot, size #m_np by #m_neq.
    N_Vector m_constraints = nullptr; //!<

    //! Indicates whether the sensitivities stored in #m_yS and #m_ySdot have been
    //! updated for the current integrator time.
    bool m_sens_ok;

    //! Maximum number of nonlinear solver iterations at one solution
    /*!
     *  If zero, this is the default of 4.
     */
    int m_maxNonlinIters = 0;

    //! Maximum number of nonlinear convergence failures
    int m_maxNonlinConvFails = -1;

    //! If true, the algebraic variables don't contribute to error tolerances
    bool m_setSuppressAlg = false;

    //! Initial IDA step size
    double m_init_step = 1e-14;
};

}

#endif
