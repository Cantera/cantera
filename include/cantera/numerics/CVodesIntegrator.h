/**
 *  @file CVodesIntegrator.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_CVODESWRAPPER_H
#define CT_CVODESWRAPPER_H

#include "cantera/numerics/Integrator.h"
#include "cantera/numerics/SundialsContext.h"
#include "cantera/base/ctexceptions.h"

#include "sundials/sundials_nvector.h"

namespace Cantera
{

/**
 * Wrapper class for 'cvodes' integrator from LLNL.
 *
 * See FuncEval.h. Classes that use CVodesIntegrator:
 * ImplicitSurfChem, ReactorNet
 */
class CVodesIntegrator : public Integrator
{
public:
    /**
     *  Constructor. Default settings: dense Jacobian, no user-supplied
     *  Jacobian function, Newton iteration.
     */
    CVodesIntegrator();
    virtual ~CVodesIntegrator();
    virtual void setTolerances(double reltol, size_t n, double* abstol);
    virtual void setTolerances(double reltol, double abstol);
    virtual void setSensitivityTolerances(double reltol, double abstol);
    virtual void initialize(double t0, FuncEval& func);
    virtual void reinitialize(double t0, FuncEval& func);
    virtual void integrate(double tout);
    virtual doublereal step(double tout);
    virtual double& solution(size_t k);
    virtual double* solution();
    virtual double* derivative(double tout, int n);
    virtual int lastOrder() const;
    virtual int nEquations() const {
        return static_cast<int>(m_neq);
    }
    virtual int nEvals() const;
    virtual void setMaxOrder(int n) {
        m_maxord = n;
    }
    virtual void setMethod(MethodType t);
    virtual void setMaxStepSize(double hmax);
    virtual void setMinStepSize(double hmin);
    virtual void setMaxSteps(int nmax);
    virtual int maxSteps();
    virtual void setMaxErrTestFails(int n);
    virtual AnyMap solverStats() const;
    void setLinearSolverType(const std::string& linSolverType) {
        m_type = linSolverType;
    }
    virtual std::string linearSolverType() const {
        return m_type;
    }
    virtual void setBandwidth(int N_Upper, int N_Lower) {
        m_mupper = N_Upper;
        m_mlower = N_Lower;
    }
    virtual int nSensParams() {
        return static_cast<int>(m_np);
    }
    virtual double sensitivity(size_t k, size_t p);
    virtual void setProblemType(int probtype);

    //! Returns a string listing the weighted error estimates associated
    //! with each solution component.
    //! This information can be used to identify which variables are
    //! responsible for integrator failures or unexpected small timesteps.
    virtual std::string getErrorInfo(int N);

    //! Error message information provide by CVodes
    std::string m_error_message;

protected:
    //! Applies user-specified options to the underlying CVODES solver. Called
    //! during integrator initialization or reinitialization.
    void applyOptions();

private:
    void sensInit(double t0, FuncEval& func);

    size_t m_neq = 0;
    void* m_cvode_mem = nullptr;
    SundialsContext m_sundials_ctx; //!< SUNContext object for Sundials>=6.0
    void* m_linsol = nullptr; //!< Sundials linear solver object
    void* m_linsol_matrix = nullptr; //!< matrix used by Sundials
    FuncEval* m_func = nullptr;
    double m_t0 = 0.0;
    double m_time; //!< The current integrator time
    N_Vector m_y = nullptr;
    N_Vector m_abstol = nullptr;
    N_Vector m_dky = nullptr;
    string m_type = "DENSE";
    int m_itol;
    int m_method;
    int m_maxord = 0;
    double m_reltol = 1e-9;
    double m_abstols = 1e-15;
    double m_reltolsens = 1e-5;
    double m_abstolsens = 1e-4;
    size_t m_nabs = 0;
    double m_hmax = 0.0;
    double m_hmin = 0.0;
    int m_maxsteps = 20000;
    int m_maxErrTestFails = 0;
    N_Vector* m_yS = nullptr;
    size_t m_np = 0;
    int m_mupper = 0;
    int m_mlower = 0;
    //! Indicates whether the sensitivities stored in m_yS have been updated
    //! for at the current integrator time.
    bool m_sens_ok = false;
};

} // namespace

#endif
