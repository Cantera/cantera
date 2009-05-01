/**
 *  @file IDA_Solver.cpp
 *
 */

// Copyright 2006-2009  California Institute of Technology

#include "IDA_Solver.h"

using namespace std;

#include <sundials_types.h>
#include <sundials_math.h>
#include <ida.h>
#include <ida_dense.h>
#include <ida_spgmr.h>
#include <ida_band.h>
#include <nvector_serial.h>

inline static N_Vector nv(void* x) {
    return reinterpret_cast<N_Vector>(x);
}

namespace Cantera {
  
  ResidData::ResidData(ResidEval* f, int npar = 0) {
    func_ = f;
  }
  
  
}


extern "C" {
  
  /**
   *  Function called by IDA to evaluate the residual, given y and
   *  ydot.  IDA allows passing in a void* pointer to access
   *  external data. Instead of requiring the user to provide a
   *  residual function directly to IDA (which would require using
   *  the sundials data types N_Vector, etc.), we define this
   *  function as the single function that IDA always calls. The
   *  real evaluation of the residual is done by an instance of a
   *  subclass of ResidEval, passed in to this function as a pointer
   *  in the parameters.
   */
  static int ida_resid(realtype t, N_Vector y, N_Vector ydot, 
		       N_Vector r, void *f_data) {
    double* ydata = NV_DATA_S(y);
    double* ydotdata = NV_DATA_S(ydot);
    double* rdata = NV_DATA_S(r);
    Cantera::ResidData* d = (Cantera::ResidData*)f_data;
    Cantera::ResidEval* f = d->func;
    f->eval(t, ydata, ydotdata, rdata);
    return 0;
  }
  
}

namespace Cantera {
  
  /**
   *  Constructor. Default settings: dense jacobian, no user-supplied
   *  Jacobian function, Newton iteration.
   */
  IDA_Solver::IDA_Solver(ResidEval& f) : neq_(0), 
					 ida_mem_(0), t0_(0.0), 
					 y_(0),
					 ydot_(0),
					 abstol_(0), 
					 type_(0), 
					 itol_(IDA_SS), 
					 iter_(0), 
					 maxord_(0),
					 reltol_(1.0e-9), 
					 abstols_(1.0e-15), 
					 nabs_(0), 
					 hmax_(0.0),
					 maxsteps_(20000), 
					 mupper_(0), 
					 mlower_(0) {}
  
  
  /// Destructor.
  IDA_Solver::~IDA_Solver() {   
    if (ida_mem_) {
      IDAFree(&ida_mem_);
    }
    if (y_) N_VDestroy_Serial(nv(y_));
    if (ydot_) N_VDestroy_Serial(nv(ydot_));
    if (abstol_) N_VDestroy_Serial(nv(abstol_));
    delete fdata_;
  }
  
  Real IDA_Solver::solution(int k) const { 
    return NV_Ith_S(nv(y_),k);
  }
  
  const Real* IDA_Solver::solutionVector() const { return NV_DATA_S(nv(y_));}
  
  Real IDA_Solver::derivative(int k) const { 
    return NV_Ith_S(nv(ydot_),k);
  }
  
  const Real* IDA_Solver::derivativeVector() const { return NV_DATA_S(nv(ydot_));}
  
  
  void IDA_Solver::setTolerances(double reltol, double* abstol) {
    itol_ = IDA_SV;
    if (abstol_) N_VDestroy_Serial(nv(abstol_));
    abstol_ = reinterpret_cast<void*>(N_VNew_Serial(neq));
    for (int i=0; i < neq; i++) {
      NV_Ith_S(nv(abstol_), i) = abstol[i];
    }
    reltol_ = reltol_; 
  }
  
  void IDA_Solver::setTolerances(double reltol, double abstol) {
    itol_ = IDA_SS;
    reltol_ = reltol;
    abstols_ = abstol;
  }
  
  void IDA_Solver::setLinearSolverType(int solverType) {
    type_ = solverType;
  }
  
  void IDA_Solver::init(double t0) 
  {
    t0_  = t0;
    
    if (y_) N_VDestroy_Serial(nv(y_));
    if (ydot_) N_VDestroy_Serial(nv(ydot_));
    if (id_) N_VDestroy_Serial(nv(id_));
    if (constraints_) N_VDestroy_Serial(nv(constraints_));
    
    y_ = reinterpret_cast<void*>(N_VNew_Serial(neq_));
    ydot_ = reinterpret_cast<void*>(N_VNew_Serial(neq_));
    constraints_ = reinterpret_cast<void*>(N_VNew_Serial(neq_));
    
    for (int i=0; i<neq_; i++) {
      NV_Ith_S(nv(y_), i) = 0.0;
      NV_Ith_S(nv(ydot_), i) = 0.0;
      NV_Ith_S(nv(constraints_), i) = 0.0;
    }
    
    // get the initial conditions
    resid.getInitialConditions(t0, NV_DATA_S(nv(y_dot)), 
			       NV_DATA_S(nv(y_)));
    
    if (ida_mem_) IDAFree(&ida_mem);
    ida_mem = IDACreate();
    
        int flag = 0;
        if (itol_ == IDA_SV) {
	  // vector atol
	  flag = IDAMalloc(ida_mem_, ida_resid_, t0_, nv(y_), nv(ydot_), 
			   itol_, reltol_, nv(abstol_));
        }
        else {
	  // scalar atol
	  flag = IDAMalloc(ida_mem_, ida_resid_, t0, nv(y_), nv(ydot_), 
			   itol_, reltol_, &abstols_);
        }
        if (flag != IDA_SUCCESS) {
	  if (flag == IDA_MEFAIL) {
	    throw IDA_Err("Memory allocation failed.");  }
	  else if (flag == IDA_ILL_INPUT) {
	    throw IDA_Err("Illegal value for IDAMalloc input argument.");
	  }
	  else 
	    throw IDA_Err("IDAMalloc failed.");
        }
	
        //-----------------------------------
        // set the linear solver type
        //-----------------------------------
	
        if (type == 1) {
	  long int N = neq_;
	  IDADense(ida_mem_, N);
        }
        else if (type == 2) {
            long int N = neq_;
            long int nu = mupper_;
            long int nl = mlower_;
            IDABand(ida_mem_, N, nu_, nl_);
        }
        else {
            throw IDA_Err("unsupported linear solver type");
        }


        // pass a pointer to func in data 
        fdata = new FuncData(&func, func.nparams());

          flag = IDASetRdata(ida_mem_, (void*)fdata_);
        if (flag != IDA_SUCCESS) 
            throw IDA_Err("IDASetRdata failed.");

        // set options
        //if (maxord > 0)
        //    flag = CVodeSetMaxOrd(cvode_mem, maxord);
        //if (maxsteps > 0)
        //    flag = CVodeSetMaxNumSteps(cvode_mem, maxsteps);
        //if (hmax > 0)
        //    flag = CVodeSetMaxStep(cvode_mem, hmax);
    }

    void IDA_Solver::solve(double tout)
    {
        double t;
        int flag;
	flag = IDASolve(ida_mem_, tout, &t, nv(y_), nv(y_dot), IDA_NORMAL);
	if (flag != IDA_SUCCESS) 
	  throw IDA_Err(" IDA error encountered.");
    }

    double IDA_Solver::step(double tout)
    {
        double t;
        int flag;
	flag = IDASolve(ida_mem, tout, &t, nv(y_), nv(y_dot), IDA_ONE_STEP);
	if (flag != IDA_SUCCESS) 
	  throw IDA_Err(" IDA error encountered.");
        return t;
      }

    Real IDA_Solver::getOutputParameter(int flag) {
        switch (flag) {
        case REAL_WORKSPACE_SIZE:
            flag = IDAGetWorkSpace(ida_mem, &lenrw, &leniw);
            return Real<(lenrw);
        }
        
}


