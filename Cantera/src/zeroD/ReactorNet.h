/**
 *  @file ReactorNet.h
 */

/*
 * $Author: hkmoffa $
 * $Revision: 1.16 $
 * $Date: 2009/07/11 17:35:05 $
 */

// Copyright 2004  California Institute of Technology

#ifndef CT_REACTORNET_H
#define CT_REACTORNET_H

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "Reactor.h"
#include "FuncEval.h"
#include "Integrator.h"
#include "Array.h"

namespace CanteraZeroD {


  class ReactorNet : public Cantera::FuncEval {

  public:

    //! Constructor
    ReactorNet();

    //! Destructor
    virtual ~ReactorNet();

    //-----------------------------------------------------

    /** @name Methods to set up a simulation. */
    //@{

    /**
     * Set initial time. Default = 0.0 s. Restarts integration
     * from this time using the current mixture state as the
     * initial condition.
     */
    void setInitialTime(doublereal time) {
      m_time = time;
      m_init = false;
    }

    /// Set the maximum time step.
    void setMaxTimeStep(double maxstep) {
      m_maxstep = maxstep;
      m_init = false;
    }
            
    void setTolerances(doublereal rtol, doublereal atol) {
      if (rtol >= 0.0) m_rtol = rtol;
      if (atol >= 0.0) m_atols = atol;
      m_init = false;
    }

    void setSensitivityTolerances(doublereal rtol, doublereal atol) {
      if (rtol >= 0.0) m_rtolsens = rtol;
      if (atol >= 0.0) m_atolsens = atol;
      m_init = false;
    }

    /// Current value of the simulation time.
    doublereal time() { return m_time; }

    /// Relative tolerance.
    doublereal rtol() { return m_rtol; }
    doublereal atol() { return m_atols; }
        
    /**
     * Initialize the reactor network. 
     */
    void initialize(doublereal t0 = 0.0);

    /**
     * Advance the state of all reactors in time.
     * @param time Time to advance to (s). 
     */
    void advance(doublereal time);

    double step(doublereal time);

    //@}

    void addReactor(ReactorBase* r, bool iown = false);

    ReactorBase& reactor(int n) {
      return *m_r[n];
    }

    bool verbose() const { return m_verbose; }
    void setVerbose(bool v = true) { m_verbose = v; }

    /// Return a reference to the integrator.
    Integrator& integrator() { return *m_integ; }

    void updateState(doublereal* y);

    double sensitivity(int k, int p) {
      return m_integ->sensitivity(k, p)/m_integ->solution(k);
    }

    double sensitivity(std::string species, int p, int reactor=0) {
      int k = globalComponentIndex(species, reactor);
      return sensitivity(k, p);
    }

    void evalJacobian(doublereal t, doublereal* y, 
		      doublereal* ydot, doublereal* p, Array2D* j);

    //-----------------------------------------------------

    // overloaded methods of class FuncEval
    virtual int neq() { return m_nv; }
    virtual void eval(doublereal t, doublereal* y, 
		      doublereal* ydot, doublereal* p);
    virtual void getInitialConditions(doublereal t0, size_t leny, 
				      doublereal* y);
    virtual int nparams() { return m_ntotpar; }

    int globalComponentIndex(std::string species, int reactor=0);

    void connect(int i, int j) {
      m_connect[j*m_nr + i] = 1;
      m_connect[i*m_nr + j] = 1;
    }

    bool connected(int i, int j) {
      return (m_connect[m_nr*i + j] == 1);
    }

  protected:

    std::vector<ReactorBase*> m_r;
    std::vector<Reactor*> m_reactors;
    int m_nr;
    int m_nreactors;
    Integrator* m_integ;
    doublereal m_time;
    bool m_init;
    int m_nv;
    vector_int m_size;
    vector_fp m_atol;
    doublereal m_rtol, m_rtolsens;
    doublereal m_atols, m_atolsens;
    doublereal m_maxstep;
    bool m_verbose;
    int m_ntotpar;
    vector_int m_nparams;
    vector_int m_connect;
    vector_fp m_ydot;

    std::vector<bool> m_iown;

  private:

  };
}

#endif

