#include "ReactorNet.h"
#include "Integrator.h"
#include "FlowDevice.h"
#include "Wall.h"

using namespace std;

namespace CanteraZeroD {

  ReactorNet::ReactorNet() : Cantera::FuncEval(), m_nr(0), m_nreactors(0),
			     m_integ(0), m_time(0.0), m_init(false), 
			     m_nv(0), m_rtol(1.0e-9), m_rtolsens(1.0e-4), 
			     m_atols(1.0e-15), m_atolsens(1.0e-4),
			     m_maxstep(-1.0),
			     m_verbose(false), m_ntotpar(0)
  {
#ifdef DEBUG_MODE
    m_verbose = true;
#endif
    m_integ = newIntegrator("CVODE");// CVodeInt;

    // use backward differencing, with a full Jacobian computed
    // numerically, and use a Newton linear iterator

    m_integ->setMethod(BDF_Method);
    m_integ->setProblemType(DENSE + NOJAC);
    m_integ->setIterator(Newton_Iter);        
  }

  ReactorNet::~ReactorNet() {
    for (int n = 0; n < m_nr; n++) {
      if (m_iown[n]) {
	delete m_r[n];
      }
      m_r[n] = 0;
    }
    m_r.clear();
    m_reactors.clear();
    deleteIntegrator(m_integ);
  }

  void ReactorNet::initialize(doublereal t0) {
    int n, nv;
    char buf[100];
    m_nv = 0;
    m_reactors.clear();
    m_nreactors = 0;
    if (m_verbose) {
      writelog("Initializing reactor network.\n");
    }
    if (m_nr == 0) 
      throw CanteraError("ReactorNet::initialize",
			 "no reactors in network!");
    for (n = 0; n < m_nr; n++) {
      if (m_r[n]->type() >= ReactorType) {
	m_r[n]->initialize(t0);
	Reactor* r = (Reactor*)m_r[n];
	m_reactors.push_back(r);
	nv = r->neq();
	m_size.push_back(nv);
	m_nparams.push_back(r->nSensParams());
	m_ntotpar += r->nSensParams();
	m_nv += nv;
	m_nreactors++;

	if (m_verbose) {
	  sprintf(buf,"Reactor %d: %d variables.\n",n,nv);
	  writelog(buf);
	  sprintf(buf,"            %d sensitivity params.\n",
		  r->nSensParams());
	  writelog(buf);
	}
	if (m_r[n]->type() == FlowReactorType && m_nr > 1) {
	  throw CanteraError("ReactorNet::initialize",
			     "FlowReactors must be used alone.");
	}
      }
    }

    m_connect.resize(m_nr*m_nr,0);
    m_ydot.resize(m_nv,0.0);
    int i, j, nin, nout, nw;
    ReactorBase *r, *rj;
    for (i = 0; i < m_nr; i++) {
      r = m_reactors[i];
      for (j = 0; j < m_nr; j++) {
	if (i == j) connect(i,j);
	else {
	  rj = m_reactors[j];
	  nin = rj->nInlets();
	  for (n = 0; n < nin; n++) {
	    if (&rj->inlet(n).out() == r) connect(i,j);
	  }
	  nout = rj->nOutlets();
	  for (n = 0; n < nout; n++) {
	    if (&rj->outlet(n).in() == r) connect(i,j);
	  }
	  nw = rj->nWalls();
	  for (n = 0; n < nw; n++) {
	    if (&rj->wall(n).left() == rj  
		&& &rj->wall(n).right() == r) connect(i,j);
	    else if (&rj->wall(n).left() == r  
		     && &rj->wall(n).right() == rj) connect(i,j);
	  }
	}
      }
    }

    m_atol.resize(neq());
    fill(m_atol.begin(), m_atol.end(), m_atols);
    m_integ->setTolerances(m_rtol, neq(), DATA_PTR(m_atol));
    m_integ->setSensitivityTolerances(m_rtolsens, m_atolsens);
    m_integ->setMaxStepSize(m_maxstep);
    if (m_verbose) {
      sprintf(buf, "Number of equations: %d\n", neq());
      writelog(buf);
      sprintf(buf, "Maximum time step:   %14.6g\n", m_maxstep);
      writelog(buf);
    }
    m_integ->initialize(t0, *this);
    m_init = true;
  }

  void ReactorNet::advance(doublereal time) {
    if (!m_init) {
      if (m_maxstep < 0.0)
	m_maxstep = time - m_time;
      initialize();
    }
    m_integ->integrate(time);
    m_time = time;
    updateState(m_integ->solution());
  }

  double ReactorNet::step(doublereal time) {
    if (!m_init) {
      if (m_maxstep < 0.0)
	m_maxstep = time - m_time;
      initialize();
    }
    m_time = m_integ->step(time);
    updateState(m_integ->solution());
    return m_time;
  }


 void ReactorNet::addReactor(ReactorBase* r, bool iown) {
   if (r->type() >= ReactorType) {
     m_r.push_back(r);
     m_iown.push_back(iown);
     m_nr++;
     if (m_verbose) {
       writelog("Adding reactor "+r->name()+"\n");
     }
   }
   else {
     if (m_verbose) {
       writelog("Not adding reactor "+r->name()+
		", since type = "+int2str(r->type())+"\n");
     }
   }
 }

  //     void ReactorNet::addSensitivityParam(int n, int stype, int i) {
  //         m_reactors[n]->addSensitivityParam(int stype, int i);
  //         m_sensreactor.push_back(n);
  //         m_nSenseParams++;
  //     }

  //     void ReactorNet::setParameters(int np, double* p) {
  //         int n, nr;
  //         for (n = 0; n < np; n++) {
  //             if (n < m_nSenseParams) {
  //                 nr = m_sensreactor[n];
  //                 m_reactors[nr]->setParameter(n, p[n]);
  //             }
  //         }
  //     }

        
  void ReactorNet::eval(doublereal t, doublereal* y, 
			doublereal* ydot, doublereal* p) {
    int n;
    int start = 0;
    int pstart = 0;
    // use a try... catch block, since exceptions are not passed
    // through CVODE, since it is C code
    try {
      updateState(y);
      for (n = 0; n < m_nreactors; n++) {
	m_reactors[n]->evalEqs(t, y + start, 
			       ydot + start, p + pstart);
	start += m_size[n];
	pstart += m_nparams[n];
      }
    }
    catch (...) {
      showErrors();
      error("Terminating execution.");
    }
  }


        
  void ReactorNet::evalJacobian(doublereal t, doublereal* y, 
				doublereal* ydot, doublereal* p, Array2D* j) {
    int n, m;
    doublereal ysave, dy;
    Array2D& jac = *j;

    // use a try... catch block, since exceptions are not passed
    // through CVODE, since it is C code
    try {
      //evaluate the unperturbed ydot
      eval(t, y, ydot, p);
      for (n = 0; n < m_nv; n++) {
             
	// perturb x(n)
	ysave = y[n];
	dy = m_atol[n] + fabs(ysave)*m_rtol;
	y[n] = ysave + dy;
	dy = y[n] - ysave;

	// calculate perturbed residual
	eval(t, y, DATA_PTR(m_ydot), p);

	// compute nth column of Jacobian
	for (m = 0; m < m_nv; m++) {
	  jac(m,n) = (m_ydot[m] - ydot[m])/dy;
	}
	y[n] = ysave;
      }
    }
    catch (...) {
      showErrors();
      error("Terminating execution.");
    }
  }

  void ReactorNet::updateState(doublereal* y) {
    int n;
    int start = 0;
    for (n = 0; n < m_nreactors; n++) {
      m_reactors[n]->updateState(y + start);
      start += m_size[n];
    }
  }

  void ReactorNet::getInitialConditions(doublereal t0, 
					size_t leny, doublereal* y) {
    int n;
    int start = 0;
    for (n = 0; n < m_nreactors; n++) {
      m_reactors[n]->getInitialConditions(t0, m_size[n], y + start);
      start += m_size[n];
    }
  }

  int ReactorNet::globalComponentIndex(string species, int reactor) {
    int start = 0;
    int n;
    for (n = 0; n < reactor; n++) start += m_size[n];
    return start + m_reactors[n]->componentIndex(species);
  }

}

