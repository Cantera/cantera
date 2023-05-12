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
 * Wrapper for Sundials IDA solver
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
    int nEvals() const override;
    void setMaxOrder(int n) override {
        m_maxord = n;
    }
    void setMaxStepSize(double hmax) override;
    void setMaxSteps(int nmax) override;
    int maxSteps() override;
    void setMaxErrTestFails(int n) override;
    void setBandwidth(int N_Upper, int N_Lower) override {
        m_mupper = N_Upper;
        m_mlower = N_Lower;
    }
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

    void setMaxNonlinIterations(int n) override;
    void setMaxNonlinConvFailures(int n) override;
    void inclAlgebraicInErrorTest(bool yesno);
    void setMethod(MethodType t) override;

protected:
    protected:
    //! Applies user-specified options to the underlying IDAS solver. Called
    //! during integrator initialization or reinitialization.
    void applyOptions();

private:
    void sensInit(double t0, FuncEval& func);
    size_t m_neq;
    void* m_ida_mem = nullptr; //!< Pointer to the IDA memory for the problem
    void* m_linsol = nullptr; //!< Sundials linear solver object
    void* m_linsol_matrix = nullptr; //!< matrix used by Sundials
    SundialsContext m_sundials_ctx; //!< SUNContext object for Sundials>=6.0
    FuncEval* m_func = nullptr;
    double m_t0 = 0.0;
    double m_time; //!< The current integrator time
    N_Vector m_y = nullptr;
    N_Vector m_ydot = nullptr;
    N_Vector m_abstol = nullptr;
    string m_type = "DENSE";
    int m_itol;
    int m_maxord = 0;
    double m_reltol = 1.0e-9;
    double m_abstols = 1.0e-15;
    double m_reltolsens, m_abstolsens;
    size_t m_nabs = 0;
    double m_hmax = 0.0;
    int m_maxsteps = 20000;
    int m_maxErrTestFails = -1;
    N_Vector* m_yS = nullptr;
    N_Vector* m_ySdot = nullptr;
    size_t m_np;
    int m_mupper = 0;
    int m_mlower = 0;
    N_Vector m_constraints = nullptr;

    //! Indicates whether the sensitivities stored in m_yS have been updated
    //! for at the current integrator time.
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

    //! Initial IDA stepsize
    double m_init_step = 1e-14;
};

}

#endif
