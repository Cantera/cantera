/**
 *  @file CVodeInt.h
 */
// Copyright 2001  California Institute of Technology

#ifndef CT_CVODEINT_H
#define CT_CVODEINT_H

#include "cantera/numerics/Integrator.h"
#include "cantera/numerics/FuncEval.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/ct_defs.h"
#include "../../ext/cvode/include/nvector.h"

namespace Cantera
{

/**
 * Exception thrown when a CVODE error is encountered.
 */
class CVodeErr : public CanteraError
{
public:
    explicit CVodeErr(const std::string& msg) : CanteraError("CVodeInt", msg) {}
};

/**
 *  Wrapper class for 'cvode' integrator from LLNL.
 *  The unmodified cvode code is in directory ext/cvode.
 *
 * @see FuncEval.h. Classes that use CVodeInt:
 * ImplicitChem, ImplicitSurfChem, Reactor
 */
class CVodeInt : public Integrator
{
public:
    /*!
     *  Constructor. Default settings: dense jacobian, no user-supplied
     *  Jacobian function, Newton iteration.
     */
    CVodeInt();
    virtual ~CVodeInt();
    virtual void setTolerances(double reltol, size_t n, double* abstol);
    virtual void setTolerances(double reltol, double abstol);
    virtual void setProblemType(int probtype);
    virtual void initialize(double t0, FuncEval& func);
    virtual void reinitialize(double t0, FuncEval& func);
    virtual void integrate(double tout);
    virtual doublereal step(double tout);
    virtual double& solution(size_t k);
    virtual double* solution();
    virtual int nEquations() const {
        return m_neq;
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
    virtual void setMaxErrTestFails(int nmax) {}

private:
    int m_neq;
    void* m_cvode_mem;
    double m_t0;
    N_Vector m_y, m_abstol;
    int m_type;
    int m_itol;
    int m_method;
    int m_iter;
    int m_maxord;
    double m_reltol;
    double m_abstols;
    int m_nabs;
    double m_hmax, m_hmin;
    int m_maxsteps;

    vector_fp m_ropt;
    long int* m_iopt;
    void* m_data;
};

}    // namespace

#endif // CT_CVODE
