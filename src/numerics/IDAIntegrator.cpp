//! @file IDAIntegrator.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/IDAIntegrator.h"
#include "cantera/base/stringUtils.h"

#include "cantera/numerics/sundials_headers.h"

using namespace std;
namespace Cantera
{

extern "C" {
    //! Function called by IDA to evaluate the residual, given y and ydot.
    /*!
     * IDA allows passing in a void* pointer to access external data. Instead of
     * requiring the user to provide a residual function directly to IDA (which
     * would require using the sundials data types N_Vector, etc.), we define
     * this function as the single function that IDA always calls. The real
     * evaluation of the residual is done by an instance of a subclass of
     * ResidEval, passed in to this function as a pointer in the parameters.
     *
     * FROM IDA WRITEUP -> What the IDA solver expects as a return flag from its
     * residual routines:
     *
     * A IDAResFn res should return a value of 0 if successful, a positive value
     * if a recoverable error occured (e.g. yy has an illegal value), or a
     * negative value if a nonrecoverable error occured. In the latter case, the
     * program halts. If a recoverable error occured, the integrator will
     * attempt to correct and retry.
     */
    static int ida_rhs(realtype t, N_Vector y, N_Vector ydot, N_Vector r, void* f_data)
    {
        FuncEval* f = (FuncEval*) f_data;
        return f->eval_nothrow(t, NV_DATA_S(y), NV_DATA_S(ydot), NV_DATA_S(r));
    }

    //! Function called by IDA when an error is encountered instead of
    //! writing to stdout. Here, save the error message provided by IDA so
    //! that it can be included in the subsequently raised CanteraError.
    static void ida_err(int error_code, const char* module,
                        const char* function, char* msg, void* eh_data)
    {
        IDAIntegrator* integrator = (IDAIntegrator*) eh_data;
        integrator->m_error_message = msg;
        integrator->m_error_message += "\n";
    }
}

IDAIntegrator::IDAIntegrator() :
    m_ida_mem(0),
    m_linsol(0),
    m_linsol_matrix(0),
    m_func(0),
    m_t0(0.0),
    m_y(0),
    m_ydot(0),
    m_abstol(0),
    m_type(DENSE+NOJAC),
    m_itol(IDA_SS),
    m_maxord(0),
    m_reltol(1.e-9),
    m_abstols(1.e-15),
    m_nabs(0),
    m_hmax(0.0),
    m_hmin(0.0),
    m_maxsteps(20000),
    m_maxErrTestFails(-1),
    m_mupper(0),
    m_mlower(0),
    m_constraints(0),
    m_maxNonlinIters(0),
    m_maxNonlinConvFails(-1),
    m_setSuppressAlg(0),
    m_init_step(1.e-14)
{
}

IDAIntegrator::~IDAIntegrator()
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
}

double& IDAIntegrator::solution(size_t k)
{
    return NV_Ith_S(m_y,k);
}

double* IDAIntegrator::solution()
{
    return NV_DATA_S(m_y);
}

void IDAIntegrator::setTolerances(double reltol, size_t n, double* abstol)
{
    m_itol = IDA_SV;
    m_nabs = n;
    if (n != m_neq) {
        if (m_abstol) {
            N_VDestroy_Serial(m_abstol);
        }
        m_abstol = N_VNew_Serial(static_cast<sd_size_t>(n));
    }
    for (size_t i=0; i<n; i++) {
        NV_Ith_S(m_abstol, i) = abstol[i];
    }
    m_reltol = reltol;
}

void IDAIntegrator::setTolerances(double reltol, double abstol)
{
    m_itol = IDA_SS;
    m_reltol = reltol;
    m_abstols = abstol;
}

void IDAIntegrator::setSensitivityTolerances(double reltol, double abstol)
{
    m_reltolsens = reltol;
    m_abstolsens = abstol;
}


void IDAIntegrator::setProblemType(int probtype)
{
    m_type = probtype;
}

void IDAIntegrator::setMaxStepSize(double hmax)
{
    m_hmax = hmax;
    if (m_ida_mem) {
        IDASetMaxStep(m_ida_mem, hmax);
    }
}

void IDAIntegrator::setMaxSteps(int nmax)
{
    m_maxsteps = nmax;
    if (m_ida_mem) {
        IDASetMaxNumSteps(m_ida_mem, m_maxsteps);
    }
}

int IDAIntegrator::maxSteps()
{
    return m_maxsteps;
}

void IDAIntegrator::setMaxErrTestFails(int n)
{
    m_maxErrTestFails = n;
    if (m_ida_mem) {
        IDASetMaxErrTestFails(m_ida_mem, n);
    }
}

void IDAIntegrator::setMaxNonlinIterations(int n)
{
    m_maxNonlinIters = n;
    if (m_ida_mem)
    {
        int flag = IDASetMaxNonlinIters(m_ida_mem, m_maxNonlinIters);
        if (flag != IDA_SUCCESS) {
            throw CanteraError("IDAIntegrator::init",
                               "IDASetmaxNonlinIters failed.");
        }
    }

    if (m_setSuppressAlg != 0) {
        int flag = IDASetSuppressAlg(m_ida_mem, m_setSuppressAlg);
        if (flag != IDA_SUCCESS) {
            throw CanteraError("IDAIntegrator::init", "IDASetSuppressAlg failed.");
        }
    }
}

void IDAIntegrator::setMaxNonlinConvFailures(int n)
{
    m_maxNonlinIters = n;
    if (m_ida_mem)
    {
        int flag = IDASetMaxConvFails(m_ida_mem, m_maxNonlinConvFails);
        if (flag != IDA_SUCCESS) {
            throw CanteraError("IDAIntegrator::init",
                               "IDASetMaxConvFails failed.");
        }
    }
}

void IDAIntegrator::inclAlgebraicInErrorTest(bool yesno)
{
    if (yesno) {
        m_setSuppressAlg = 0;
    } else {
        m_setSuppressAlg = 1;
    }

    if (m_ida_mem)
    {
        int flag = IDASetSuppressAlg(m_ida_mem, m_setSuppressAlg);
        if (flag != IDA_SUCCESS) {
            throw CanteraError("IDAIntegrator::init", "IDASetSuppressAlg failed.");
        }
    }
}

void IDAIntegrator::initialize(double t0, FuncEval& func)
{
    m_neq = func.neq();
    m_t0 = t0;
    m_time = t0;
    m_func = &func;
    func.clearErrors();

    if (m_y) {
        N_VDestroy_Serial(m_y); // free solution vector if already allocated
    }
    m_y = N_VNew_Serial(static_cast<sd_size_t>(m_neq)); // allocate solution vector
    N_VConst(0.0, m_y);

    if (m_ydot)
    {
        N_VDestroy_Serial(m_ydot); // free derivative vector if already allocated
    }
    m_ydot = N_VNew_Serial(m_neq);
    N_VConst(0.0, m_ydot);

    // check abs tolerance array size
    if (m_itol == IDA_SV && m_nabs < m_neq) {
        throw CanteraError("IDAIntegrator::initialize",
                           "not enough absolute tolerance values specified.");
    }

    if (m_constraints) {
        N_VDestroy_Serial(m_constraints);
    }
    m_constraints = N_VNew_Serial(static_cast<sd_size_t>(m_neq));
    // set the constraints
    func.getConstraints(NV_DATA_S(m_constraints));


    // get the initial conditions
    func.getState(NV_DATA_S(m_y), NV_DATA_S(m_ydot));

    if (m_ida_mem) {
        IDAFree(&m_ida_mem);
    }

    //! Create the IDA solver
    m_ida_mem = IDACreate();
    if (!m_ida_mem) {
        throw CanteraError("IDAIntegrator::initialize",
                           "IDACreate failed.");
    }

    int flag = IDAInit(m_ida_mem, ida_rhs, m_t0, m_y, m_ydot);
    if (flag != IDA_SUCCESS) {
        if (flag == IDA_MEM_FAIL) {
            throw CanteraError("IDAIntegrator::initialize",
                               "Memory allocation failed.");
        } else if (flag == IDA_ILL_INPUT) {
            throw CanteraError("IDAIntegrator::initialize",
                               "Illegal value for IDAInit input argument.");
        } else {
            throw CanteraError("IDAIntegrator::initialize",
                               "IDAInit failed.");
        }
    }

    // set constraints
    flag = IDASetId(m_ida_mem, m_constraints);
    if (flag != IDA_SUCCESS) {
        if (flag == IDA_MEM_NULL) {
            throw CanteraError("IDAIntegrator::initialize",
                               "Memory allocation failed.");
        } else {
            throw CanteraError("IDAIntegrator::initialize",
                               "IDASetId failed.");
        }
    }

    flag = IDASetErrHandlerFn(m_ida_mem, &ida_err, this);
    if (flag != IDA_SUCCESS) {
        throw CanteraError("IDAIntegrator::initialize",
                            "IDASetErrHandlerFn failed.");
    }

    if (m_itol == IDA_SV) {
        flag = IDASVtolerances(m_ida_mem, m_reltol, m_abstol);
    } else {
        flag = IDASStolerances(m_ida_mem, m_reltol, m_abstols);
    }
    if (flag != IDA_SUCCESS) {
        if (flag == IDA_MEM_FAIL) {
            throw CanteraError("IDAIntegrator::initialize",
                               "Memory allocation failed.");
        } else if (flag == IDA_ILL_INPUT) {
            throw CanteraError("IDAIntegrator::initialize",
                               "Illegal value for IDAInit input argument.");
        } else {
            throw CanteraError("IDAIntegrator::initialize",
                               "IDAInit failed.");
        }
    }

    flag = IDASetUserData(m_ida_mem, &func);
    if (flag != IDA_SUCCESS) {
        throw CanteraError("IDAIntegrator::initialize",
                           "IDASetUserData failed.");
    }
    applyOptions();
}

void IDAIntegrator::reinitialize(double t0, FuncEval& func)
{
    m_t0 = t0;
    m_time = t0;
    func.getState(NV_DATA_S(m_y));
    m_func = &func;
    func.clearErrors();

    int result = IDAReInit(m_ida_mem, m_t0, m_y, m_ydot);
    if (result != IDA_SUCCESS) {
        throw CanteraError("IDAIntegrator::reinitialize",
                           "IDAReInit failed. result = {}", result);
    }
    applyOptions();
}

void IDAIntegrator::applyOptions()
{
    if (m_type == DENSE + NOJAC) {
        sd_size_t N = static_cast<sd_size_t>(m_neq);
        #if CT_SUNDIALS_VERSION >= 30
            SUNLinSolFree((SUNLinearSolver) m_linsol);
            SUNMatDestroy((SUNMatrix) m_linsol_matrix);
            m_linsol_matrix = SUNDenseMatrix(N, N);
            #if CT_SUNDIALS_USE_LAPACK
                m_linsol = SUNLapackDense(m_y, (SUNMatrix) m_linsol_matrix);
            #else
                m_linsol = SUNDenseLinearSolver(m_y, (SUNMatrix) m_linsol_matrix);
            #endif
            IDADlsSetLinearSolver(m_ida_mem, (SUNLinearSolver) m_linsol,
                                  (SUNMatrix) m_linsol_matrix);
        #else
            #if CT_SUNDIALS_USE_LAPACK
                IDALapackDense(m_ida_mem, N);
            #else
                IDADense(m_ida_mem, N);
            #endif
        #endif
    } else if (m_type == DIAG) {
        throw CanteraError("IDAIntegrator::applyOptions",
                           "Cannot use a diagonal matrix with IDA.");
    } else if (m_type == GMRES) {
        #if CT_SUNDIALS_VERSION >= 30
            m_linsol = SUNSPGMR(m_y, PREC_NONE, 0);
            IDASpilsSetLinearSolver(m_ida_mem, (SUNLinearSolver) m_linsol);
        #else
            IDASpgmr(m_ida_mem, PREC_NONE, 0);
        #endif
    } else if (m_type == BAND + NOJAC) {
        sd_size_t N = static_cast<sd_size_t>(m_neq);
        long int nu = m_mupper;
        long int nl = m_mlower;
        #if CT_SUNDIALS_VERSION >= 30
            SUNLinSolFree((SUNLinearSolver) m_linsol);
            SUNMatDestroy((SUNMatrix) m_linsol_matrix);
            m_linsol_matrix = SUNBandMatrix(N, nu, nl, nu+nl);
            #if CT_SUNDIALS_USE_LAPACK
                m_linsol = SUNLapackBand(m_y, (SUNMatrix) m_linsol_matrix);
            #else
                m_linsol = SUNBandLinearSolver(m_y, (SUNMatrix) m_linsol_matrix);
            #endif
            IDADlsSetLinearSolver(m_ida_mem, (SUNLinearSolver) m_linsol,
                                 (SUNMatrix) m_linsol_matrix);
        #else
            #if CT_SUNDIALS_USE_LAPACK
                IDALapackBand(m_ida_mem, N, nu, nl);
            #else
                IDABand(m_ida_mem, N, nu, nl);
            #endif
        #endif
    } else {
        throw CanteraError("IDAIntegrator::applyOptions",
                           "unsupported option");
    }

    if (m_init_step > 0) {
        IDASetInitStep(m_ida_mem, m_init_step);
    }

    if (m_maxord > 0) {
        IDASetMaxOrd(m_ida_mem, m_maxord);
    }
    if (m_maxsteps > 0) {
        IDASetMaxNumSteps(m_ida_mem, m_maxsteps);
    }
    if (m_hmax > 0) {
        IDASetMaxStep(m_ida_mem, m_hmax);
    }
    if (m_maxNonlinIters > 0) {
        int flag = IDASetMaxNonlinIters(m_ida_mem, m_maxNonlinIters);
        if (flag != IDA_SUCCESS) {
            throw CanteraError("IDAIntegrator::init",
                               "IDASetmaxNonlinIters failed.");
        }
    }
    if (m_maxNonlinConvFails > 0) {
        int flag = IDASetMaxConvFails(m_ida_mem, m_maxNonlinConvFails);
        if (flag != IDA_SUCCESS) {
            throw CanteraError("IDAIntegrator::init",
                               "IDASetMaxConvFails failed.");
        }
    }
    if (m_setSuppressAlg != 0) {
        int flag = IDASetSuppressAlg(m_ida_mem, m_setSuppressAlg);
        if (flag != IDA_SUCCESS) {
            throw CanteraError("IDAIntegrator::init", "IDASetSuppressAlg failed.");
        }
    }
}

void IDAIntegrator::sensInit(double t0, FuncEval& func)
{
    m_np = func.nparams();
    m_sens_ok = false;

    N_Vector y = N_VNew_Serial(static_cast<sd_size_t>(func.neq()));
    m_yS = N_VCloneVectorArray_Serial(static_cast<sd_size_t>(m_np), y);
    for (size_t n = 0; n < m_np; n++) {
        N_VConst(0.0, m_yS[n]);
    }
    N_VDestroy_Serial(y);
    N_Vector ydot = N_VNew_Serial(static_cast<sd_size_t>(func.neq()));
    m_ySdot = N_VCloneVectorArray_Serial(static_cast<sd_size_t>(m_np), ydot);
    for (size_t n = 0; n < m_np; n++) {
        N_VConst(0.0, m_ySdot[n]);
    }

    int flag = IDASensInit(m_ida_mem, static_cast<sd_size_t>(m_np),
                           IDA_STAGGERED, IDASensResFn(0), m_yS, m_ySdot);

    if (flag != IDA_SUCCESS) {
        throw CanteraError("IDAIntegrator::sensInit", "Error in IDASensInit");
    }
    vector_fp atol(m_np);
    for (size_t n = 0; n < m_np; n++) {
        // This scaling factor is tuned so that reaction and species enthalpy
        // sensitivities can be computed simultaneously with the same abstol.
        atol[n] = m_abstolsens / func.m_paramScales[n];
    }
    flag = IDASensSStolerances(m_ida_mem, m_reltolsens, atol.data());
    if (flag != IDA_SUCCESS) {
        throw CanteraError("IDAIntegrator::sensInit", "Error in IDASensSStolerances");
    }
}

void IDAIntegrator::integrate(double tout)
{
    if (tout == m_time) {
        return;
    }
    int flag = IDASolve(m_ida_mem, tout, &m_time, m_y, m_ydot, IDA_NORMAL);
    if (flag != IDA_SUCCESS) {
        string f_errs = m_func->getErrors();
        if (!f_errs.empty()) {
            f_errs = "Exceptions caught during RHS evaluation:\n" + f_errs;
        }
        throw CanteraError("IDAIntegrator::integrate",
            "IDA error encountered. Error code: {}\n{}\n"
            "{}"
            "Components with largest weighted error estimates:\n{}",
            flag, m_error_message, f_errs, getErrorInfo(10));
    }
}

double IDAIntegrator::step(double tout)
{
    int flag = IDASolve(m_ida_mem, tout, &m_time, m_y, m_ydot, IDA_ONE_STEP);
    if (flag != IDA_SUCCESS) {
        string f_errs = m_func->getErrors();
        if (!f_errs.empty()) {
            f_errs = "Exceptions caught during RHS evaluation:\n" + f_errs;
        }
        throw CanteraError("IDAIntegrator::step",
            "IDA error encountered. Error code: {}\n{}\n"
            "{}"
            "Components with largest weighted error estimates:\n{}",
            flag, f_errs, m_error_message, getErrorInfo(10));

    }
    return m_time;
}

int IDAIntegrator::nEvals() const
{
    // neither IDAGetNumRhsEvals, IDADlsGetNumRhsEvals, IDASpilsGetNumRhsEvals,
    // or IDAGetNumLinResEvals seem to exist
    throw CanteraError("IDAIntegrator::nEvals",
                       "Not implemented.");
}

double IDAIntegrator::sensitivity(size_t k, size_t p)
{
    if (m_time == m_t0) {
        // calls to IDAGetSens are only allowed after a successful time step.
        return 0.0;
    }
    if (!m_sens_ok && m_np) {
        int flag = IDAGetSens(m_ida_mem, &m_time, m_yS);
        if (flag != IDA_SUCCESS) {
            throw CanteraError("IDAIntegrator::sensitivity",
                               "IDAGetSens failed. Error code: {}", flag);
        }
        m_sens_ok = true;
    }

    if (k >= m_neq) {
        throw CanteraError("IDAIntegrator::sensitivity",
                           "sensitivity: k out of range ({})", k);
    }
    if (p >= m_np) {
        throw CanteraError("IDAIntegrator::sensitivity",
                           "sensitivity: p out of range ({})", p);
    }
    return NV_Ith_S(m_yS[p],k);
}

string IDAIntegrator::getErrorInfo(int N)
{
    N_Vector errs = N_VNew_Serial(static_cast<sd_size_t>(m_neq));
    N_Vector errw = N_VNew_Serial(static_cast<sd_size_t>(m_neq));
    IDAGetErrWeights(m_ida_mem, errw);
    IDAGetEstLocalErrors(m_ida_mem, errs);

    vector<tuple<double, double, size_t> > weightedErrors;
    for (size_t i=0; i<m_neq; i++) {
        double err = NV_Ith_S(errs, i) * NV_Ith_S(errw, i);
        weightedErrors.emplace_back(-abs(err), err, i);
    }
    N_VDestroy(errs);
    N_VDestroy(errw);

    N = std::min(N, static_cast<int>(m_neq));
    sort(weightedErrors.begin(), weightedErrors.end());
    fmt::memory_buffer s;
    for (int i=0; i<N; i++) {
        format_to(s, "{}: {}\n",
                get<2>(weightedErrors[i]), get<1>(weightedErrors[i]));
    }
    return to_string(s);
}

void IDAIntegrator::setMethod(MethodType t)
{
    if (t != BDF_Method) {
        // IDA only has the BDF method
        throw CanteraError("IDAIntegrator::setMethod", "unknown method");
    }
}

void IDAIntegrator::setIterator(IterType t)
{
    if (t != Newton_Iter) {
        // IDA only has the newton iterator
        throw CanteraError("IDAIntegrator::setIterator", "unknown iterator");
    }
}

}
