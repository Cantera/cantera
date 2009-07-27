/**
 *  @file DASPK.cpp
 *
 */
/*
 *  $Date: 2009/03/27 23:39:26 $
 *  $Revision: 1.2 $
 */

// Copyright 2001  California Institute of Technology

// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "DASPK.h"
#include "ctexceptions.h"
#include "stringUtils.h"

#include <iostream>
using namespace std;

extern "C" {

  typedef void (*ResidFunc)(const doublereal* t, 
			    const doublereal* y, const doublereal* yprime,
			    const doublereal* cj, doublereal* delta, 
			    integer* ires, doublereal* rpar, integer* ipar);

  typedef void (*JacFunc)();
  typedef void (*PsolFunc)();

  extern void ddaspk_(ResidFunc res, integer* neq, doublereal* t, 
		      doublereal* y, doublereal* yprime, doublereal* tout, integer* info,
		      doublereal* rtol, doublereal* atol, integer* idid, doublereal* rwork,
		      integer* lrw, integer* iwork, integer* liw, doublereal* rpar, 
		      integer* ipar, JacFunc jac, PsolFunc psol);


  /**
   *  Function called by DASPK to evaluate the residual. 
   */
  static void ddaspk_res(const doublereal* t, 
			 const doublereal* y, const doublereal* yprime,
			 const doublereal* cj, doublereal* delta, 
			 integer* ires, doublereal* rpar, integer* ipar) {
    void **iddres_res = reinterpret_cast<void **>(&(ipar[0]));
    void *hndl = *iddres_res;
    Cantera::ResidEval* f = (Cantera::ResidEval*)hndl;
    double delta_t = 0.0;
    f->evalResid(*t, delta_t, y, yprime, delta);
  }

  static void ddaspk_jac() {}
  static void ddaspk_psol() {}

}


namespace Cantera {

  class DASPKErr : public CanteraError {
  public:
    DASPKErr(string proc, string msg) 
      : CanteraError("DASPK::"+proc,msg) {}
    virtual ~DASPKErr(){}
  };

  /**
   *  Constructor. Default settings: dense jacobian, no user-supplied
   *  Jacobian function, Newton iteration.
   */
  DASPK::DASPK(ResidEval& f) : 
    m_resid(f), 
    m_idid(0),
    m_lrw(0), 
    m_liw(0), 
    m_ml(0), 
    m_mu(0), 
    m_lenwp(0),
    m_ok(false), 
    m_init(false), 
    m_time(0.0)
  {
    m_info.resize(20);
    m_neq = f.neq();
    m_rwork.resize(20);  // will be reset later
    m_iwork.resize(20);  //  "
    m_ipar.resize(2);
    m_rpar.resize(2);
    void *iddr = static_cast<void *>(&m_resid);
    void **iddr_ipar = reinterpret_cast<void **>(&(m_ipar[0]));
    *iddr_ipar = iddr;
    setTolerances(1.e-7, 1.e-15);
  }


  /// Destructor.
  DASPK::~DASPK(){}
    
  void DASPK::setTolerances(int nr, double* reltol, int na, double* abstol) {
    // scalar tolerances
    if (nr == 1 && na == 1) {
      setInfo(2,0);
      m_rtol.resize(1);
      m_rtol[0] = reltol[0];
      m_atol.resize(1);
      m_atol[0] = abstol[0];
    }
    // vector tolerances
    else {
      setInfo(2,1);
      m_rtol.resize(neq());
      m_atol.resize(neq());
      copy(reltol, reltol + nr, m_rtol.begin());
      copy(abstol, abstol + na, m_atol.begin());
    }
  }

  void DASPK::setTolerances(double reltol, double abstol) {
    doublereal rtol = reltol;
    doublereal atol = abstol;
    setTolerances(1, &rtol, 1, &atol);
  }

  void DASPK::setJacobian(Jacobian& jac) {
        
    // No Jacobian evaluation function is supplied, so let DASPK
    // compute the Jacobian by numerical finite-difference
    if (!jac.supplied()) setInfo(5,0);
    else {
      setInfo(5,1);
    }

    if (jac.isBanded()) { 
      setInfo(6,1);
      setIwork(1, jac.lowerBandWidth());
      setIwork(2, jac.upperBandWidth());
    }
    else setInfo(6,0);
  }

  void DASPK::setMethod(int methodType) {
    if (methodType == cDirect) 
      setInfo(12,0);
    else if (methodType == cKrylov) 
      setInfo(12,1);
    else 
      throw DASPKErr("setMethod",
		     "method must be either cDirect "
		     "or cKrylov");
  }

  void DASPK::setMaxTime(doublereal tmax) {
    setInfo(4,1);
    setRwork(1,tmax);
  }

  void DASPK::setMaxStepSize(doublereal dtmax) {
    setInfo(7,1);
    setRwork(2,dtmax);
  }

  void DASPK::setInitialIntStepSize(doublereal h0) {
    setInfo(8,1);
    setRwork(3,h0);
  }

  void DASPK::setMaxOrder(int n) {
    setInfo(9,1);
    setIwork(3,n);
  }

  void DASPK::estimateInitial_Y_given_Yp() {
    setInfo(11,2);
  }

  void DASPK::estimateInitial_YaYp_given_Yd(
					    const vector<int>& vartypes) {
    setInfo(11,2);
    int m, n = neq();
    int lid = ((info(10) == 0 || info(10) == 2) ? 
	       41 : 41 + neq());
    if (int(m_iwork.size()) < lid + neq())
      m_iwork.resize(lid + neq());
    for (m = 0; m < n; m++) {
      setIwork(lid + m, vartypes[m]);
    }
  }

  void DASPK::sizeRwork() {
    int base;
    if (info(12) == 0) {
      base = 50 + 9*neq();
      if (info(6) == 0) 
	base += neq()*neq();
      else {
	base += (2*m_ml + m_mu + 1)*neq();
	if (info(5) == 0)
	  base += 2*(neq()/(m_ml + m_mu + 1) + 1);
      }
    }
    else {
      base = 91 + 18*neq() + m_lenwp;
    } 
    if (info(16) == 1) base += neq();

    /// @todo fix this!
    base = 2000000; // tmp

    m_rwork.resize(base, 0.0);
    m_lrw = base;
  }

  void DASPK::sizeIwork() {
    int base;
    if (info(12) == 0) {
      base = 40 + neq();
    }
    else {
      base = 40 + m_lenwp;
    }
    if (info(10) == 1 || info(10) == 3) base += neq();
    if (info(11) == 1 || info(16) == 1) base += neq();
    m_iwork.resize(base);
    m_liw = base;
  }


  void DASPK::init(doublereal t0) 
  {
    m_init = true;
    m_time = t0;
    setInfo(1,0);  // tells DASPK to initialize
    sizeRwork();
    sizeIwork();
    //m_resid.init(t0);
  }

  int DASPK::integrate(doublereal tout) {
    if (!m_init) init(0.0);

    doublereal tfinal = tout;
    setInfo(3,0);  // don't want intermediate output

    ddaspk_(ddaspk_res, &m_neq, &m_time, m_resid.solution(), 
            m_resid.solution_dot(), &tfinal, m_info.begin(), 
            m_rtol.begin(), m_atol.begin(), &m_idid, 
            m_rwork.begin(), &m_lrw, m_iwork.begin(), &m_liw,
            m_rpar.begin(), m_ipar.begin(), ddaspk_jac, ddaspk_psol);

    return m_idid;
  }

  void DASPK::step(double tout)
  {
    setInfo(3,1);  // do want intermediate output
    doublereal tfinal = tout;
    //        setInfo(3,0);  // don't want intermediate output

    ddaspk_(ddaspk_res, &m_neq, &m_time, m_resid.solution(), 
            m_resid.solution_dot(), &tfinal, m_info.begin(), 
            m_rtol.begin(), m_atol.begin(), &m_idid, 
            m_rwork.begin(), &m_lrw, m_iwork.begin(), &m_liw,
            m_rpar.begin(), m_ipar.begin(), ddaspk_jac, ddaspk_psol);
    if (m_idid < 0) {
      throw DASPKErr("step",
		     "DASPK returned IDID = "+int2str(m_idid));
      m_ok = false;
    }
    else if (m_idid == 1 || m_idid == 2 || m_idid == 3) {
      m_ok = true;
    }
    else {
      m_ok = false;
    }
    return;
  }

  int DASPK::nEvals() const { return iwork(12); }
}



