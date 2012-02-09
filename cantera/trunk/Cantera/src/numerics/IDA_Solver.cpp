
/**
 *  @file IDA_Solver.cpp
 *
 */

// Copyright 2006  California Institute of Technology

#include "IDA_Solver.h"
#include "stringUtils.h"

#include <iostream>
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
    
    /**
     * A simple class to hold an array of parameter values and a pointer to 
     * an instance of a subclass of ResidEval.
     */
    class ResidData {




    public:
        ResidData(ResidEval* f, int npar = 0) {
            m_func = f;
        }
        virtual ~ResidData() {}
        ResidEval* m_func;
    };
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
        Cantera::ResidEval* f = d->m_func;
        f->eval(t, ydata, ydotdata, rdata);
        return 0;
    }

}

namespace Cantera {


    /**
     *  Constructor. Default settings: dense jacobian, no user-supplied
     *  Jacobian function, Newton iteration.
     */
    IDA_Solver::IDA_Solver(ResidEval& f) : DAE_Solver(f), 
                                           m_neq(0), 
                                           m_ida_mem(0), 
                                           m_t0(0.0), 
                                           m_y(0),
                                           m_ydot(0),
                                           m_abstol(0), 
                                           m_type(0), 
                                           m_itol(IDA_SS), 
                                           m_iter(0), 
                                           m_maxord(0),
                                           m_reltol(1.e-9), 
                                           m_abstols(1.e-15), 
                                           m_nabs(0), 
                                           m_hmax(0.0),
                                           m_maxsteps(20000), 
                                           m_mupper(0), 
                                           m_mlower(0) {}


    /// Destructor.
    IDA_Solver::~IDA_Solver()
    {   
        if (m_ida_mem) {
            IDAFree(&m_ida_mem);
        }
      if (m_y) N_VDestroy_Serial(nv(m_y));
      if (m_ydot) N_VDestroy_Serial(nv(m_ydot));
      if (m_abstol) N_VDestroy_Serial(nv(m_abstol));
      delete m_fdata;
    }
    
    doublereal IDA_Solver::solution(int k) const { 
        return NV_Ith_S(nv(m_y),k);
    }
 
    const doublereal* IDA_Solver::solutionVector() const { return NV_DATA_S(nv(m_y));}

    doublereal IDA_Solver::derivative(int k) const { 
        return NV_Ith_S(nv(m_ydot),k);
    }
 
    const doublereal* IDA_Solver::derivativeVector() const { return NV_DATA_S(nv(m_ydot));}


    void IDA_Solver::setTolerances(double reltol, double* abstol) {
        m_itol = IDA_SV;
        if (m_abstol) N_VDestroy_Serial(nv(m_abstol));
        m_abstol = reinterpret_cast<void*>(N_VNew_Serial(m_neq));
        for (int i=0; i < m_neq; i++) {
            NV_Ith_S(nv(m_abstol), i) = abstol[i];
        }
        m_reltol = reltol; 
    }

    void IDA_Solver::setTolerances(double reltol, double abstol) {
        m_itol = IDA_SS;
        m_reltol = reltol;
        m_abstols = abstol;
    }

    void IDA_Solver::setLinearSolverType(int solverType) {
        m_type = solverType;
    }

    void IDA_Solver::init(double t0) 
    {
        m_t0  = t0;

        if (m_y) N_VDestroy_Serial(nv(m_y));
        if (m_ydot) N_VDestroy_Serial(nv(m_ydot));
        if (m_id) N_VDestroy_Serial(nv(m_id));
        if (m_constraints) N_VDestroy_Serial(nv(m_constraints));

        m_y = reinterpret_cast<void*>(N_VNew_Serial(m_neq));
        m_ydot = reinterpret_cast<void*>(N_VNew_Serial(m_neq));
        m_constraints = reinterpret_cast<void*>(N_VNew_Serial(m_neq));

        for (int i=0; i<m_neq; i++) {
            NV_Ith_S(nv(m_y), i) = 0.0;
            NV_Ith_S(nv(m_ydot), i) = 0.0;
            NV_Ith_S(nv(m_constraints), i) = 0.0;
        }

        // get the initial conditions
        m_resid.getInitialConditions(m_t0, NV_DATA_S(nv(m_ydot)), 
            NV_DATA_S(nv(m_y)));
        
        if (m_ida_mem) IDAFree(&m_ida_mem);
        m_ida_mem = IDACreate();

        int flag = 0;
        if (m_itol == IDA_SV) {
            // vector atol
            flag = IDAMalloc(m_ida_mem, ida_resid, m_t0, nv(m_y), nv(m_ydot), 
                m_itol, m_reltol, nv(m_abstol));
        }
        else {
            // scalar atol
            flag = IDAMalloc(m_ida_mem, ida_resid, m_t0, nv(m_y), nv(m_ydot), 
                m_itol, m_reltol, &m_abstols);
        }
        if (flag != IDA_SUCCESS) {
            if (flag == IDA_MEM_FAIL) {
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

        if (m_type == 1) {
            long int N = m_neq;
            IDADense(m_ida_mem, N);
        }
        else if (m_type == 2) {
            long int N = m_neq;
            long int nu = m_mupper;
            long int nl = m_mlower;
            IDABand(m_ida_mem, N, nu, nl);
        }
        else {
            throw IDA_Err("unsupported linear solver type");
        }


        // pass a pointer to func in m_data 
        m_fdata = new ResidData(&func, func.nparams());

        flag = IDASetRdata(m_ida_mem, (void*)m_fdata);
        if (flag != IDA_SUCCESS) 
            throw IDA_Err("IDASetRdata failed.");

        // set options
        //if (m_maxord > 0)
        //    flag = CVodeSetMaxOrd(m_cvode_mem, m_maxord);
        //if (m_maxsteps > 0)
        //    flag = CVodeSetMaxNumSteps(m_cvode_mem, m_maxsteps);
        //if (m_hmax > 0)
        //    flag = CVodeSetMaxStep(m_cvode_mem, m_hmax);
    }

    void IDA_Solver::solve(double tout)
    {
        double t;
        int flag;
	flag = IDASolve(m_ida_mem, tout, &t, nv(m_y), nv(m_ydot), IDA_NORMAL);
	if (flag != IDA_SUCCESS) 
	  throw IDA_Err(" IDA error encountered.");
    }

    double IDA_Solver::step(double tout)
    {
        double t;
        int flag;
	flag = IDASolve(m_ida_mem, tout, &t, nv(m_y), nv(m_ydot), IDA_ONE_STEP);
	if (flag != IDA_SUCCESS) 
	  throw IDA_Err(" IDA error encountered.");
        return t;
      }

    doublereal IDA_Solver::getOutputParameter(int flag) {
        switch (flag) {
        case REAL_WORKSPACE_SIZE:
            flag = IDAGetWorkSpace(m_ida_mem, &lenrw, &leniw);
            return doublereal(lenrw);
        }
        
}


