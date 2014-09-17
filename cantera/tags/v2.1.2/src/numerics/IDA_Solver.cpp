/**
 *  @file IDA_Solver.cpp
 */

// Copyright 2006  California Institute of Technology

#include "cantera/numerics/IDA_Solver.h"
#include "cantera/base/stringUtils.h"

#if SUNDIALS_VERSION >= 24
#include "sundials/sundials_types.h"
#include "sundials/sundials_math.h"
#include "ida/ida.h"
#include "ida/ida_dense.h"
#include "ida/ida_spgmr.h"
#include "ida/ida_band.h"
#include "nvector/nvector_serial.h"

using namespace std;

#if SUNDIALS_VERSION < 25
typedef int sd_size_t;
#else
typedef long int sd_size_t;
#endif

namespace Cantera
{

/**
 * A simple class to hold an array of parameter values and a pointer to
 * an instance of a subclass of ResidEval.
 */
class ResidData
{

public:

    ResidData(ResidJacEval* f, IDA_Solver* s, int npar = 0) {
        m_func = f;
        m_solver = s;
    }

    virtual ~ResidData() {
    }

    ResidJacEval* m_func;
    IDA_Solver* m_solver;
};
}

extern "C" {
    //!  Function called by IDA to evaluate the residual, given y and ydot.
    /*!
     *  IDA allows passing in a void* pointer to access external data. Instead of requiring the user to provide a
     *  residual function directly to IDA (which would require using
     *  the sundials data types N_Vector, etc.), we define this function as the single function that IDA always calls. The
     *  real evaluation of the residual is done by an instance of a subclass of ResidEval, passed in to this
     *  function as a pointer in the parameters.
     *
     * FROM IDA WRITEUP -> What the IDA solver expects as a return flag from its residual routines ------
     *   A IDAResFn res should return a value of 0 if successful, a positive
     *   value if a recoverable error occured (e.g. yy has an illegal value),
     *   or a negative value if a nonrecoverable error occured. In the latter
     *   case, the program halts. If a recoverable error occured, the integrator
     *   will attempt to correct and retry.
     */
    static int ida_resid(realtype t, N_Vector y, N_Vector ydot, N_Vector r, void* f_data)
    {
        double* ydata = NV_DATA_S(y);
        double* ydotdata = NV_DATA_S(ydot);
        double* rdata = NV_DATA_S(r);
        Cantera::ResidData* d = (Cantera::ResidData*) f_data;
        Cantera::ResidJacEval* f = d->m_func;
        Cantera::IDA_Solver* s = d->m_solver;
        double delta_t = s->getCurrentStepFromIDA();
        // TODO evaluate  evalType. Assumed to be Base_ResidEval
        int retn = 0;
        int flag = f->evalResidNJ(t, delta_t, ydata, ydotdata, rdata);
        if (flag < 0) {
            // This signals to IDA that a nonrecoverable error has occurred.
            retn = flag;
        }
        return retn;
    }

    //! Function called by by IDA to evaluate the Jacobian, given y and ydot.
    /*!
     * typedef int (*IDADlsDenseJacFn)(sd_size_t N, realtype t, realtype c_j,
     *                             N_Vector y, N_Vector yp, N_Vector r,
     *                             DlsMat Jac, void *user_data,
     *                             N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
     *
     * A IDADlsDenseJacFn should return
     *     0 if successful,
     *     a positive int if a recoverable error occurred, or
     *     a negative int if a nonrecoverable error occurred.
     * In the case of a recoverable error return, the integrator will
     * attempt to recover by reducing the stepsize (which changes cj).
     */
    static int ida_jacobian(sd_size_t nrows, realtype t, realtype c_j,  N_Vector y, N_Vector ydot, N_Vector r,
                            DlsMat Jac,  void* f_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
    {
        doublereal* ydata    = NV_DATA_S(y);
        doublereal* ydotdata = NV_DATA_S(ydot);
        doublereal* rdata    = NV_DATA_S(r);
        Cantera::ResidData* d = (Cantera::ResidData*) f_data;
        Cantera::ResidJacEval* f = d->m_func;
        doublereal* const* colPts = Jac->cols;
        Cantera::IDA_Solver* s = d->m_solver;
        double delta_t = s->getCurrentStepFromIDA();
        // printf(" delta_t = %g 1/cj = %g\n", delta_t, 1.0/c_j);
        f->evalJacobianDP(t, delta_t, c_j,  ydata, ydotdata, colPts, rdata);
        return 0;
    }
}

namespace Cantera
{

IDA_Solver::IDA_Solver(ResidJacEval& f) :
    DAE_Solver(f),
    m_ida_mem(0),
    m_t0(0.0),
    m_y(0),
    m_ydot(0),
    m_id(0),
    m_constraints(0),
    m_abstol(0),
    m_type(0),
    m_itol(IDA_SS),
    m_iter(0),
    m_reltol(1.e-9),
    m_abstols(1.e-15),
    m_nabs(0),
    m_hmax(0.0),
    m_hmin(0.0),
    m_h0(0.0),
    m_maxsteps(20000),
    m_maxord(0),
    m_formJac(0),
    m_tstop(0.0),
    m_told_old(0.0),
    m_told(0.0),
    m_tcurrent(0.0),
    m_deltat(0.0),
    m_maxErrTestFails(-1),
    m_maxNonlinIters(0),
    m_maxNonlinConvFails(-1),
    m_setSuppressAlg(0),
    m_fdata(0),
    m_mupper(0),
    m_mlower(0)
{
}

IDA_Solver::~IDA_Solver()
{
    if (m_ida_mem) {
        IDAFree(&m_ida_mem);
    }
    if (m_y) {
        N_VDestroy_Serial(m_y);
    }
    if (m_ydot) {
        N_VDestroy_Serial(m_ydot);
    }
    if (m_abstol) {
        N_VDestroy_Serial(m_abstol);
    }
    if (m_constraints) {
        N_VDestroy_Serial(m_constraints);
    }
    delete m_fdata;
}

doublereal IDA_Solver::solution(int k) const
{
    return NV_Ith_S(m_y,k);
}

const doublereal* IDA_Solver::solutionVector() const
{
    return NV_DATA_S(m_y);
}

doublereal IDA_Solver::derivative(int k) const
{
    return NV_Ith_S(m_ydot,k);
}

const doublereal* IDA_Solver::derivativeVector() const
{
    return NV_DATA_S(m_ydot);
}

void IDA_Solver::setTolerances(double reltol, double* abstol)
{
    m_itol = IDA_SV;
    if (!m_abstol) {
        m_abstol = N_VNew_Serial(m_neq);
    }
    for (int i = 0; i < m_neq; i++) {
        NV_Ith_S(m_abstol, i) = abstol[i];
    }
    m_reltol = reltol;
    if (m_ida_mem) {
        int flag = IDASVtolerances(m_ida_mem, m_reltol, m_abstol);
        if (flag != IDA_SUCCESS) {
            throw IDA_Err("Memory allocation failed.");
        }
    }
}

void IDA_Solver::setTolerances(doublereal reltol, doublereal abstol)
{
    m_itol = IDA_SS;
    m_reltol = reltol;
    m_abstols = abstol;
    if (m_ida_mem) {
        int flag = IDASStolerances(m_ida_mem, m_reltol, m_abstols);
        if (flag != IDA_SUCCESS) {
            throw IDA_Err("Memory allocation failed.");
        }
    }
}

void IDA_Solver::setLinearSolverType(int solverType)
{
    m_type = solverType;
}

void IDA_Solver::setDenseLinearSolver()
{
    setLinearSolverType(0);
}

void IDA_Solver::setBandedLinearSolver(int m_upper, int m_lower)
{
    m_type = 2;
    m_upper = m_mupper;
    m_mlower = m_lower;
}

void IDA_Solver::setMaxOrder(int n)
{
    m_maxord = n;
}

void IDA_Solver::setMaxNumSteps(int n)
{
    m_maxsteps = n;
}

void IDA_Solver::setInitialStepSize(doublereal h0)
{
    m_h0 = h0;
}

void IDA_Solver::setStopTime(doublereal tstop)
{
    m_tstop = tstop;
}

doublereal IDA_Solver::getCurrentStepFromIDA()
{
    doublereal hcur;
    IDAGetCurrentStep(m_ida_mem, &hcur);
    return hcur;
}

void IDA_Solver::setJacobianType(int formJac)
{
    m_formJac = formJac;
    if (m_ida_mem) {
        if (m_formJac == 1) {
            int flag = IDADlsSetDenseJacFn(m_ida_mem, ida_jacobian);
            if (flag != IDA_SUCCESS) {
                throw IDA_Err("IDADlsSetDenseJacFn failed.");
            }
        }
    }
}

void IDA_Solver::setMaxErrTestFailures(int maxErrTestFails)
{
    m_maxErrTestFails = maxErrTestFails;
}

void IDA_Solver::setMaxNonlinIterations(int n)
{
    m_maxNonlinIters = n;
}

void IDA_Solver::setMaxNonlinConvFailures(int n)
{
    m_maxNonlinConvFails = n;
}

void IDA_Solver::inclAlgebraicInErrorTest(bool yesno)
{
    if (yesno) {
        m_setSuppressAlg = 0;
    } else {
        m_setSuppressAlg = 1;
    }
}

void IDA_Solver::init(doublereal t0)
{
    m_t0 = t0;
    m_told = t0;
    m_told_old = t0;
    m_tcurrent = t0;
    if (m_y) {
        N_VDestroy_Serial(m_y);
    }
    if (m_ydot) {
        N_VDestroy_Serial(m_ydot);
    }
    if (m_id) {
        N_VDestroy_Serial(m_id);
    }
    if (m_constraints) {
        N_VDestroy_Serial(m_constraints);
    }

    m_y = N_VNew_Serial(m_neq);
    m_ydot = N_VNew_Serial(m_neq);
    m_constraints = N_VNew_Serial(m_neq);

    for (int i=0; i<m_neq; i++) {
        NV_Ith_S(m_y, i) = 0.0;
        NV_Ith_S(m_ydot, i) = 0.0;
        NV_Ith_S(m_constraints, i) = 0.0;
    }

    // get the initial conditions
    m_resid.getInitialConditions(m_t0, NV_DATA_S(m_y), NV_DATA_S(m_ydot));

    if (m_ida_mem) {
        IDAFree(&m_ida_mem);
    }

    /* Call IDACreate */
    m_ida_mem = IDACreate();

    int flag = 0;



    if (m_itol == IDA_SV) {
#if SUNDIALS_VERSION <= 23
        // vector atol
        flag = IDAMalloc(m_ida_mem, ida_resid, m_t0, m_y, m_ydot,
                         m_itol, m_reltol, m_abstol);
        if (flag != IDA_SUCCESS) {
            if (flag == IDA_MEM_FAIL) {
                throw IDA_Err("Memory allocation failed.");
            } else if (flag == IDA_ILL_INPUT) {
                throw IDA_Err("Illegal value for IDAMalloc input argument.");
            }  else {
                throw IDA_Err("IDAMalloc failed.");
            }
        }

#elif SUNDIALS_VERSION >= 24
        flag = IDAInit(m_ida_mem, ida_resid, m_t0, m_y, m_ydot);
        if (flag != IDA_SUCCESS) {
            if (flag == IDA_MEM_FAIL) {
                throw IDA_Err("Memory allocation failed.");
            } else if (flag == IDA_ILL_INPUT) {
                throw IDA_Err("Illegal value for IDAMalloc input argument.");
            } else {
                throw IDA_Err("IDAMalloc failed.");
            }
        }
        flag = IDASVtolerances(m_ida_mem, m_reltol, m_abstol);
        if (flag != IDA_SUCCESS) {
            throw IDA_Err("Memory allocation failed.");
        }
#endif
    } else {
#if SUNDIALS_VERSION <= 23
        // scalar atol
        flag = IDAMalloc(m_ida_mem, ida_resid, m_t0, m_y, m_ydot,
                         m_itol, m_reltol, &m_abstols);
        if (flag != IDA_SUCCESS) {
            if (flag == IDA_MEM_FAIL) {
                throw IDA_Err("Memory allocation failed.");
            } else if (flag == IDA_ILL_INPUT) {
                throw IDA_Err("Illegal value for IDAMalloc input argument.");
            } else {
                throw IDA_Err("IDAMalloc failed.");
            }
        }

#elif SUNDIALS_VERSION >= 24
        flag = IDAInit(m_ida_mem, ida_resid, m_t0, m_y, m_ydot);
        if (flag != IDA_SUCCESS) {
            if (flag == IDA_MEM_FAIL) {
                throw IDA_Err("Memory allocation failed.");
            } else if (flag == IDA_ILL_INPUT) {
                throw IDA_Err("Illegal value for IDAMalloc input argument.");
            } else {
                throw IDA_Err("IDAMalloc failed.");
            }
        }
        flag = IDASStolerances(m_ida_mem, m_reltol, m_abstols);
        if (flag != IDA_SUCCESS) {
            throw IDA_Err("Memory allocation failed.");
        }
#endif
    }

    //-----------------------------------
    // set the linear solver type
    //-----------------------------------

    if (m_type == 1 || m_type == 0) {
        long int N = m_neq;
        flag = IDADense(m_ida_mem, N);
        if (flag) {
            throw IDA_Err("IDADense failed");
        }
    } else if (m_type == 2) {
        long int N = m_neq;
        long int nu = m_mupper;
        long int nl = m_mlower;
        IDABand(m_ida_mem, N, nu, nl);
    } else {
        throw IDA_Err("unsupported linear solver type");
    }

    if (m_formJac == 1) {
        flag = IDADlsSetDenseJacFn(m_ida_mem, ida_jacobian);
        if (flag != IDA_SUCCESS) {
            throw IDA_Err("IDADlsSetDenseJacFn failed.");
        }
    }

    // pass a pointer to func in m_data
    m_fdata = new ResidData(&m_resid, this, m_resid.nparams());
#if SUNDIALS_VERSION <= 23
    flag = IDASetRdata(m_ida_mem, (void*)m_fdata);
    if (flag != IDA_SUCCESS) {
        throw IDA_Err("IDASetRdata failed.");
    }
#elif SUNDIALS_VERSION >= 24
    flag = IDASetUserData(m_ida_mem, (void*)m_fdata);
    if (flag != IDA_SUCCESS) {
        throw IDA_Err("IDASetUserData failed.");
    }
#endif

    // set options
    if (m_maxord > 0) {
        flag = IDASetMaxOrd(m_ida_mem, m_maxord);
        if (flag != IDA_SUCCESS) {
            throw IDA_Err("IDASetMaxOrd failed.");
        }
    }
    if (m_maxsteps > 0) {
        flag = IDASetMaxNumSteps(m_ida_mem, m_maxsteps);
        if (flag != IDA_SUCCESS) {
            throw IDA_Err("IDASetMaxNumSteps failed.");
        }
    }
    if (m_h0 > 0.0) {
        flag = IDASetInitStep(m_ida_mem, m_h0);
        if (flag != IDA_SUCCESS) {
            throw IDA_Err("IDASetInitStep failed.");
        }
    }
    if (m_tstop > 0.0) {
        flag = IDASetStopTime(m_ida_mem, m_tstop);
        if (flag != IDA_SUCCESS) {
            throw IDA_Err("IDASetStopTime failed.");
        }
    }
    if (m_maxErrTestFails >= 0) {
        flag = IDASetMaxErrTestFails(m_ida_mem, m_maxErrTestFails);
        if (flag != IDA_SUCCESS) {
            throw IDA_Err("IDASetMaxErrTestFails failed.");
        }
    }
    if (m_maxNonlinIters > 0) {
        flag = IDASetMaxNonlinIters(m_ida_mem, m_maxNonlinIters);
        if (flag != IDA_SUCCESS) {
            throw IDA_Err("IDASetmaxNonlinIters failed.");
        }
    }
    if (m_maxNonlinConvFails >= 0) {
        flag = IDASetMaxConvFails(m_ida_mem, m_maxNonlinConvFails);
        if (flag != IDA_SUCCESS) {
            throw IDA_Err("IDASetMaxConvFails failed.");
        }
    }
    if (m_setSuppressAlg != 0) {
        flag = IDASetSuppressAlg(m_ida_mem, m_setSuppressAlg);
        if (flag != IDA_SUCCESS) {
            throw IDA_Err("IDASetSuppressAlg failed.");
        }
    }
}

void IDA_Solver::correctInitial_Y_given_Yp(doublereal* y, doublereal* yp,  doublereal tout)
{
    int icopt = IDA_Y_INIT;
    doublereal tout1 = tout;
    if (tout == 0.0) {
        double h0 = 1.0E-5;
        if (m_h0 > 0.0) {
            h0 = m_h0;
        }
        tout1 = m_t0 + h0;
    }

    int flag = IDACalcIC(m_ida_mem, icopt, tout1);
    if (flag != IDA_SUCCESS) {
        throw IDA_Err("IDACalcIC failed: error = " + int2str(flag));
    }


    flag = IDAGetConsistentIC(m_ida_mem, m_y, m_ydot);
    if (flag != IDA_SUCCESS) {
        throw IDA_Err("IDAGetSolution failed: error = " + int2str(flag));
    }
    doublereal* yy = NV_DATA_S(m_y);
    doublereal* yyp = NV_DATA_S(m_ydot);

    for (int i = 0; i < m_neq; i++) {
        y[i]  = yy[i];
        yp[i] = yyp[i];
    }
}

void IDA_Solver::correctInitial_YaYp_given_Yd(doublereal* y, doublereal* yp, doublereal tout)
{
    int icopt = IDA_YA_YDP_INIT;
    doublereal tout1 = tout;
    if (tout == 0.0) {
        double h0 = 1.0E-5;
        if (m_h0 > 0.0) {
            h0 = m_h0;
        }
        tout1 = m_t0 + h0;
    }

    int flag = IDACalcIC(m_ida_mem, icopt, tout1);
    if (flag != IDA_SUCCESS) {
        throw IDA_Err("IDACalcIC failed: error = " + int2str(flag));
    }


    flag = IDAGetConsistentIC(m_ida_mem, m_y, m_ydot);
    if (flag != IDA_SUCCESS) {
        throw IDA_Err("IDAGetSolution failed: error = " + int2str(flag));
    }
    doublereal* yy = NV_DATA_S(m_y);
    doublereal* yyp = NV_DATA_S(m_ydot);

    for (int i = 0; i < m_neq; i++) {
        y[i]  = yy[i];
        yp[i] = yyp[i];
    }
}

int IDA_Solver::solve(double tout)
{
    double tretn;
    int flag;
    flag = IDASetStopTime(m_ida_mem, tout);
    if (flag != IDA_SUCCESS) {
        throw IDA_Err(" IDA error encountered.");
    }
    do {
        if (tout <= m_tcurrent) {
            throw IDA_Err(" tout <= tcurrent");
        }
        m_told_old = m_told;
        m_told = m_tcurrent;
        flag = IDASolve(m_ida_mem, tout, &tretn, m_y, m_ydot, IDA_ONE_STEP);
        if (flag < 0) {
            throw IDA_Err(" IDA error encountered.");
        } else   if (flag == IDA_TSTOP_RETURN) {
            // we've reached our goal, and have actually integrated past it
        } else if (flag == IDA_ROOT_RETURN) {
            // not sure what to do with this yet
        } else if (flag == IDA_WARNING) {
            throw IDA_Err(" IDA Warning encountered.");
        }
        m_tcurrent = tretn;
        m_deltat = m_tcurrent - m_told;
    } while (tretn < tout);

    if (flag != IDA_SUCCESS && flag != IDA_TSTOP_RETURN) {
        throw IDA_Err(" IDA error encountered.");
    }
    return flag;
}

double IDA_Solver::step(double tout)
{
    double t;
    int flag;
    if (tout <= m_tcurrent) {
        throw IDA_Err(" tout <= tcurrent");
    }
    m_told_old = m_told;
    m_told = m_tcurrent;
    flag = IDASolve(m_ida_mem, tout, &t, m_y, m_ydot, IDA_ONE_STEP);
    if (flag < 0) {
        throw IDA_Err(" IDA error encountered.");
    } else   if (flag == IDA_TSTOP_RETURN) {
        // we've reached our goal, and have actually integrated past it
    } else if (flag == IDA_ROOT_RETURN) {
        // not sure what to do with this yet
    } else if (flag == IDA_WARNING) {
        throw IDA_Err(" IDA Warning encountered.");
    }
    m_tcurrent = t;
    m_deltat = m_tcurrent - m_told;
    return t;
}

doublereal IDA_Solver::getOutputParameter(int flag) const
{
    long int lenrw, leniw;
    switch (flag) {
    case REAL_WORKSPACE_SIZE:
        flag = IDAGetWorkSpace(m_ida_mem, &lenrw, &leniw);
        return doublereal(lenrw);
        break;
    }
    return 0.0;
}

}
#endif
