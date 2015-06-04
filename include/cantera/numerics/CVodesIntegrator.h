/**
 *  @file CVodesIntegrator.h
 */
// Copyright 2005  California Institute of Technology

#ifndef CT_CVODESWRAPPER_H
#define CT_CVODESWRAPPER_H

#include "cantera/numerics/Integrator.h"
#include "cantera/base/ctexceptions.h"

#ifdef HAS_SUNDIALS

#include "sundials/sundials_nvector.h"

namespace Cantera
{

class FuncData;

/**
 * Exception thrown when a CVODES error is encountered.
 */
class CVodesErr : public CanteraError
{
public:
    explicit CVodesErr(const std::string& msg) : CanteraError("CVodesIntegrator", msg) {}
};

/**
 * Wrapper class for 'cvodes' integrator from LLNL.
 *
 * @see FuncEval.h. Classes that use CVodeInt:
 * ImplicitChem, ImplicitSurfChem, Reactor
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
    virtual void setProblemType(int probtype);
    virtual void initialize(double t0, FuncEval& func);
    virtual void reinitialize(double t0, FuncEval& func);
    virtual void integrate(double tout);
    virtual doublereal step(double tout);
    virtual double& solution(size_t k);
    virtual double* solution();
    virtual int nEquations() const {
        return static_cast<int>(m_neq);
    }
    virtual int nEvals() const;
    virtual void setMaxOrder(int n) {
        m_maxord = n;
    }
    virtual void setMethod(MethodType t);
    virtual void setIterator(IterType t);
    virtual void setMaxStepSize(double hmax);
    virtual void setMinStepSize(double hmin);
    virtual void setMaxSteps(int nmax);
    virtual void setMaxErrTestFails(int n);
    virtual void setBandwidth(int N_Upper, int N_Lower) {
        m_mupper = N_Upper;
        m_mlower = N_Lower;
    }
    virtual int nSensParams() {
        return static_cast<int>(m_np);
    }
    virtual double sensitivity(size_t k, size_t p);

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

    size_t m_neq;
    void* m_cvode_mem;
    double m_t0;
    double m_time; //!< The current integrator time
    N_Vector m_y, m_abstol;
    int m_type;
    int m_itol;
    int m_method;
    int m_iter;
    int m_maxord;
    double m_reltol;
    double m_abstols;
    double m_reltolsens, m_abstolsens;
    size_t m_nabs;
    double m_hmax, m_hmin;
    int m_maxsteps;
    int m_maxErrTestFails;
    FuncData* m_fdata;
    N_Vector*  m_yS;
    size_t m_np;
    int m_mupper, m_mlower;

    //! Indicates whether the sensitivities stored in m_yS have been updated
    //! for at the current integrator time.
    bool m_sens_ok;

};

}    // namespace

#else

#error No sundials!

#endif

#endif
