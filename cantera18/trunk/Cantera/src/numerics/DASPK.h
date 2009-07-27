/**
 *
 *  @file DASPK.h
 *
 *  Header file for class DASPK
 */

/* 
 *  $Date: 2009/03/27 23:39:26 $
 *  $Revision: 1.3 $
 *
 *  Copyright 2001 California Institute of Technology
 *
 */

#ifndef CT_DASPK_H
#define CT_DASPK_H

#include <vector>

#include "ct_defs.h"
#include "ResidEval.h"

namespace Cantera {

  class Jacobian {
  public:
    Jacobian(){}
    virtual ~Jacobian(){}
    virtual bool supplied() { return false; }
    virtual bool isBanded() { return false; }
    virtual int lowerBandWidth() { return 0; }
    virtual int upperBandWidth() { return 0; }
  };

  class BandedJac : public Jacobian {
  public:
    BandedJac(int ml, int mu) {
      m_ml = ml; m_mu = mu;
    }
    virtual bool supplied() { return false; }
    virtual bool isBanded() { return true; }
    virtual int lowerBandWidth() { return m_ml; }
    virtual int upperBandWidth() { return m_mu; }
  protected:
    int m_ml, m_mu;
  };

  class ResidEval;

  const int cDirect = 0;
  const int cKrylov = 1;

  /**
   * Wrapper for DASPK 2.0 DAE solver of Petzold et al.
   */
  class DASPK {
  public:

    DASPK(ResidEval& f);
    virtual ~DASPK();

    integer iwork(int n) const { return m_iwork[n-1];}
    doublereal rwork(int n) const { return m_rwork[n-1];}
    void setIwork(int n, integer m) { m_iwork[n-1] = m; }
    void setRwork(int n, doublereal v) { m_rwork[n-1] = v; }
    integer info(int n) { return m_info[n-1]; }
    void setInfo(int n, integer m) { m_info[n-1] = m; }

    void setTolerances(int nr, double* reltol, int na, double* abstol);
    void setTolerances(double reltol, double abstol);
    void setJacobian(Jacobian& jac);
    void setMethod(int methodType);
    void setMaxTime(doublereal tmax);
    void setMaxStepSize(doublereal dtmax);
    void setMaxOrder(int n);
    void setInitialIntStepSize(doublereal h0);
    void estimateInitial_Y_given_Yp();
    void estimateInitial_YaYp_given_Yd(const std::vector<int>& vartypes);
    void sizeRwork();
    void sizeIwork();
    int integrate(doublereal tout);
    void step(doublereal tout);
    int nEvals() const;
    int neq() { return m_resid.nEquations(); }
    void init(doublereal t0);

  protected:

    ResidEval& m_resid;
    vector_int m_info;
    vector_int m_iwork;
    vector_int m_ipar;

    vector_fp  m_rwork;
    vector_fp  m_atol;
    vector_fp  m_rtol;
    vector_fp  m_rpar;

    integer    m_idid;
    integer    m_neq;
    integer    m_lrw;
    integer    m_liw;
    integer    m_ml;
    integer    m_mu;
    integer    m_lenwp;

    bool       m_ok;
    bool       m_init;
    doublereal m_time;
  };

}

#endif

