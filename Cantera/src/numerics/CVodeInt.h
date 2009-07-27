/**
 *  @file CVodeInt.h
 */

/* $Author: hkmoffa $
 * $Date: 2009/07/17 15:32:51 $
 * $Revision: 1.1 $
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_CVODEINT_H
#define CT_CVODEINT_H

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "Integrator.h"
#include "FuncEval.h"
#include "ctexceptions.h"
#include "ct_defs.h"

namespace Cantera {

  /**
   * Exception thrown when a CVODE error is encountered.
   */
  class CVodeErr : public CanteraError {
  public:
    CVodeErr(std::string msg) : CanteraError("CVodeInt", msg){}
  };


  /**
   *  Wrapper class for 'cvode' integrator from LLNL.
   *  The unmodified cvode code is in directory ext/cvode.
   *
   * @see FuncEval.h. Classes that use CVodeInt:
   * ImplicitChem, ImplicitSurfChem, Reactor
   *
   */
  class CVodeInt : public Integrator {

  public:

    CVodeInt();
    virtual ~CVodeInt();
    virtual void setTolerances(double reltol, int n, double* abstol);
    virtual void setTolerances(double reltol, double abstol);
    virtual void setProblemType(int probtype);
    virtual void initialize(double t0, FuncEval& func);
    virtual void reinitialize(double t0, FuncEval& func);
    virtual void integrate(double tout);
    virtual doublereal step(double tout);
    virtual double& solution(int k);
    virtual double* solution();
    virtual int nEquations() const { return m_neq;}
    virtual int nEvals() const;
    virtual void setMaxOrder(int n) { m_maxord = n; }
    virtual void setMethod(MethodType t);
    virtual void setIterator(IterType t);
    virtual void setMaxStepSize(double hmax);
    virtual void setMinStepSize(double hmin);
    virtual void setMaxSteps(int nmax);

  private:

    int m_neq;
    void* m_cvode_mem;
    double m_t0;
    void *m_y, *m_abstol;
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
