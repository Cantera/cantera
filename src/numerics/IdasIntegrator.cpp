//! @file IdasIntegrator.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/IdasIntegrator.h"
#include "cantera/base/stringUtils.h"

#include "cantera/numerics/sundials_headers.h"

using namespace std;

namespace {

N_Vector newNVector(size_t N, Cantera::SundialsContext& context)
{
#if CT_SUNDIALS_VERSION >= 60
    return N_VNew_Serial(static_cast<sd_size_t>(N), context.get());
#else
    return N_VNew_Serial(static_cast<sd_size_t>(N));
#endif
}

} // end anonymous namespace

namespace Cantera
{

extern "C" {

//! Function called by IDA to evaluate the residual, given y and ydot.
/*!
* IDA allows passing in a void* pointer to access external data. Instead of requiring
* the user to provide a residual function directly to IDA (which would require using the
* Sundials data types N_Vector, etc.), we define this function as the single function
* that IDA always calls. The real evaluation of the residual is done by an instance of a
* subclass of ResidEval, passed in to this function as a pointer in the parameters.
*
* FROM IDA WRITEUP -> What the IDA solver expects as a return flag from its
* residual routines:
*
* A IDAResFn res should return a value of 0 if successful, a positive value if a
* recoverable error occurred (e.g. yy has an illegal value), or a negative value if a
* nonrecoverable error occurred. In the latter case, the program halts. If a recoverable
* error occurred, the integrator will attempt to correct and retry.
*/
static int ida_rhs(realtype t, N_Vector y, N_Vector ydot, N_Vector r, void* f_data)
{
    FuncEval* f = (FuncEval*) f_data;
    return f->evalDaeNoThrow(t, NV_DATA_S(y), NV_DATA_S(ydot), NV_DATA_S(r));
}

//! Function called by IDA when an error is encountered instead of writing to stdout.
//! Here, save the error message provided by IDA so that it can be included in the
//! subsequently raised CanteraError.
static void ida_err(int error_code, const char* module,
                    const char* function, char* msg, void* eh_data)
{
    IdasIntegrator* integrator = (IdasIntegrator*) eh_data;
    integrator->m_error_message = msg;
    integrator->m_error_message += "\n";
}

}

IdasIntegrator::IdasIntegrator()
    : m_itol(IDA_SS)
{
}

IdasIntegrator::~IdasIntegrator()
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

double& IdasIntegrator::solution(size_t k)
{
    return NV_Ith_S(m_y, k);
}

double* IdasIntegrator::solution()
{
    return NV_DATA_S(m_y);
}

void IdasIntegrator::setTolerances(double reltol, size_t n, double* abstol)
{
    m_itol = IDA_SV;
    if (n != m_nabs) {
        m_nabs = n;
        if (m_abstol) {
            N_VDestroy_Serial(m_abstol);
        }
        m_abstol = newNVector(static_cast<sd_size_t>(n), m_sundials_ctx);
    }
    for (size_t i=0; i<n; i++) {
        NV_Ith_S(m_abstol, i) = abstol[i];
    }
    m_reltol = reltol;
}

void IdasIntegrator::setTolerances(double reltol, double abstol)
{
    m_itol = IDA_SS;
    m_reltol = reltol;
    m_abstols = abstol;
}

void IdasIntegrator::setSensitivityTolerances(double reltol, double abstol)
{
    m_reltolsens = reltol;
    m_abstolsens = abstol;
}


void IdasIntegrator::setLinearSolverType(const string& linearSolverType)
{
    m_type = linearSolverType;
}

void IdasIntegrator::setMaxOrder(int n)
{
    if (m_ida_mem) {
        int flag = IDASetMaxOrd(m_ida_mem, n);
        checkError(flag, "setMaxOrder", "IDASetMaxOrd");
    }
    m_maxord = n;
}

void IdasIntegrator::setMaxStepSize(double hmax)
{
    m_hmax = hmax;
    if (m_ida_mem) {
        int flag = IDASetMaxStep(m_ida_mem, hmax);
        checkError(flag, "setMaxStepSize", "IDASetMaxStep");
    }
}

void IdasIntegrator::setMaxSteps(int nmax)
{
    m_maxsteps = nmax;
    if (m_ida_mem) {
        IDASetMaxNumSteps(m_ida_mem, m_maxsteps);
    }
}

int IdasIntegrator::maxSteps()
{
    return m_maxsteps;
}

void IdasIntegrator::setMaxErrTestFails(int n)
{
    m_maxErrTestFails = n;
    if (m_ida_mem) {
        IDASetMaxErrTestFails(m_ida_mem, n);
    }
}

AnyMap IdasIntegrator::solverStats() const
{
    AnyMap stats;
    long int val;
    int lastOrder;

    int flag = IDAGetNumSteps(m_ida_mem, &val);
    checkError(flag, "solverStats", "IDAGetNumSteps");
    stats["steps"] = val;
    IDAGetNumResEvals(m_ida_mem, &val);
    stats["res_evals"] = val;
    IDAGetNumLinSolvSetups(m_ida_mem, &val);
    stats["lin_solve_setups"] = val;
    IDAGetNumErrTestFails(m_ida_mem, &val);
    stats["err_tests_fails"] = val;
    IDAGetLastOrder(m_ida_mem, &lastOrder);
    stats["last_order"] = lastOrder;
    IDAGetNumNonlinSolvIters(m_ida_mem, &val);
    stats["nonlinear_iters"] = val;
    IDAGetNumNonlinSolvConvFails(m_ida_mem, &val);
    stats["nonlinear_conv_fails"] = val;
    return stats;
}

void IdasIntegrator::setMaxNonlinIterations(int n)
{
    m_maxNonlinIters = n;
    if (m_ida_mem) {
        int flag = IDASetMaxNonlinIters(m_ida_mem, m_maxNonlinIters);
        checkError(flag, "setMaxNonlinIterations", "IDASetMaxNonlinIters");
    }
}

void IdasIntegrator::setMaxNonlinConvFailures(int n)
{
    m_maxNonlinConvFails = n;
    if (m_ida_mem) {
        int flag = IDASetMaxConvFails(m_ida_mem, m_maxNonlinConvFails);
        checkError(flag, "setMaxNonlinConvFailures", "IDASetMaxConvFails");
    }
}

void IdasIntegrator::includeAlgebraicInErrorTest(bool yesno)
{
    m_setSuppressAlg = !yesno;
    if (m_ida_mem) {
        int flag = IDASetSuppressAlg(m_ida_mem, m_setSuppressAlg);
        checkError(flag, "inclAlgebraicInErrorTest", "IDASetSuppressAlg");
    }
}

void IdasIntegrator::initialize(double t0, FuncEval& func)
{
    m_neq = func.neq();
    m_t0 = t0;
    m_time = t0;
    m_func = &func;
    func.clearErrors();

    if (m_y) {
        N_VDestroy_Serial(m_y); // free solution vector if already allocated
    }
    m_y = newNVector(static_cast<sd_size_t>(m_neq), m_sundials_ctx);
    N_VConst(0.0, m_y);

    if (m_ydot) {
        N_VDestroy_Serial(m_ydot); // free derivative vector if already allocated
    }
    m_ydot = newNVector(m_neq, m_sundials_ctx);
    N_VConst(0.0, m_ydot);

    // check abs tolerance array size
    if (m_itol == IDA_SV && m_nabs < m_neq) {
        throw CanteraError("IdasIntegrator::initialize",
                           "not enough absolute tolerance values specified.");
    }

    if (m_constraints) {
        N_VDestroy_Serial(m_constraints);
    }
    m_constraints = newNVector(static_cast<sd_size_t>(m_neq), m_sundials_ctx);
    // set the constraints
    func.getConstraints(NV_DATA_S(m_constraints));

    // get the initial conditions
    func.getStateDae(NV_DATA_S(m_y), NV_DATA_S(m_ydot));

    if (m_ida_mem) {
        IDAFree(&m_ida_mem);
    }

    //! Create the IDA solver
    #if CT_SUNDIALS_VERSION >= 60
        m_ida_mem = IDACreate(m_sundials_ctx.get());
    #else
        m_ida_mem = IDACreate();
    #endif
    if (!m_ida_mem) {
        throw CanteraError("IdasIntegrator::initialize", "IDACreate failed.");
    }

    int flag = IDAInit(m_ida_mem, ida_rhs, m_t0, m_y, m_ydot);
    if (flag != IDA_SUCCESS) {
        if (flag == IDA_MEM_FAIL) {
            throw CanteraError("IdasIntegrator::initialize",
                               "Memory allocation failed.");
        } else if (flag == IDA_ILL_INPUT) {
            throw CanteraError("IdasIntegrator::initialize",
                               "Illegal value for IDAInit input argument.");
        } else {
            throw CanteraError("IdasIntegrator::initialize", "IDAInit failed.");
        }
    }

    flag = IDASetErrHandlerFn(m_ida_mem, &ida_err, this);

    // set constraints
    flag = IDASetId(m_ida_mem, m_constraints);
    checkError(flag, "initialize", "IDASetId");

    if (m_itol == IDA_SV) {
        flag = IDASVtolerances(m_ida_mem, m_reltol, m_abstol);
        checkError(flag, "initialize", "IDASVtolerances");
    } else {
        flag = IDASStolerances(m_ida_mem, m_reltol, m_abstols);
        checkError(flag, "initialize", "IDASStolerances");
    }

    flag = IDASetUserData(m_ida_mem, &func);
    checkError(flag, "initialize", "IDASetUserData");

    if (func.nparams() > 0) {
        throw CanteraError("IdasIntegrator::initialize", "Sensitivity analysis "
                           "for DAE systems is not fully implemented");
        sensInit(t0, func);
        flag = IDASetSensParams(m_ida_mem, func.m_sens_params.data(),
                                func.m_paramScales.data(), NULL);
        checkError(flag, "initialize", "IDASetSensParams");
    }
    applyOptions();
}

void IdasIntegrator::reinitialize(double t0, FuncEval& func)
{
    m_t0 = t0;
    m_time = t0;
    func.getStateDae(NV_DATA_S(m_y), NV_DATA_S(m_ydot));
    m_func = &func;
    func.clearErrors();

    int result = IDAReInit(m_ida_mem, m_t0, m_y, m_ydot);
    checkError(result, "reinitialize", "IDAReInit");
    applyOptions();
}

void IdasIntegrator::applyOptions()
{
    if (m_type == "DENSE") {
        sd_size_t N = static_cast<sd_size_t>(m_neq);
        SUNLinSolFree((SUNLinearSolver) m_linsol);
        SUNMatDestroy((SUNMatrix) m_linsol_matrix);
        #if CT_SUNDIALS_VERSION >= 60
            m_linsol_matrix = SUNDenseMatrix(N, N, m_sundials_ctx.get());
        #else
            m_linsol_matrix = SUNDenseMatrix(N, N);
        #endif
        #if CT_SUNDIALS_VERSION >= 60
            m_linsol_matrix = SUNDenseMatrix(N, N, m_sundials_ctx.get());
        #else
            m_linsol_matrix = SUNDenseMatrix(N, N);
        #endif
        #if CT_SUNDIALS_USE_LAPACK
            #if CT_SUNDIALS_VERSION >= 60
                m_linsol = SUNLinSol_LapackDense(m_y, (SUNMatrix) m_linsol_matrix,
                                                 m_sundials_ctx.get());
            #else
                m_linsol = SUNLapackDense(m_y, (SUNMatrix) m_linsol_matrix);
            #endif
        #else
            #if CT_SUNDIALS_VERSION >= 60
                m_linsol = SUNLinSol_Dense(m_y, (SUNMatrix) m_linsol_matrix,
                                           m_sundials_ctx.get());
            #else
                m_linsol = SUNLinSol_Dense(m_y, (SUNMatrix) m_linsol_matrix);
            #endif
        #endif
        #if CT_SUNDIALS_VERSION >= 60
            IDASetLinearSolver(m_ida_mem, (SUNLinearSolver) m_linsol,
                               (SUNMatrix) m_linsol_matrix);
        #else
            IDADlsSetLinearSolver(m_ida_mem, (SUNLinearSolver) m_linsol,
                                  (SUNMatrix) m_linsol_matrix);
        #endif
    } else if (m_type == "GMRES") {
        #if CT_SUNDIALS_VERSION >= 60
            m_linsol = SUNLinSol_SPGMR(m_y, PREC_NONE, 0, m_sundials_ctx.get());
            IDASetLinearSolver(m_ida_mem, (SUNLinearSolver) m_linsol, nullptr);
        #elif CT_SUNDIALS_VERSION >= 40
            m_linsol = SUNLinSol_SPGMR(m_y, PREC_NONE, 0);
            IDASpilsSetLinearSolver(m_ida_mem, (SUNLinearSolver) m_linsol);
        #else
            m_linsol = SUNSPGMR(m_y, PREC_NONE, 0);
            IDASpilsSetLinearSolver(m_ida_mem, (SUNLinearSolver) m_linsol);
        #endif
    } else {
        throw CanteraError("IdasIntegrator::applyOptions",
                           "unsupported linear solver flag '{}'", m_type);
    }

    if (m_init_step > 0) {
        IDASetInitStep(m_ida_mem, m_init_step);
    }

    if (m_maxord > 0) {
        int flag = IDASetMaxOrd(m_ida_mem, m_maxord);
        checkError(flag, "applyOptions", "IDASetMaxOrd");
    }
    if (m_maxsteps > 0) {
        IDASetMaxNumSteps(m_ida_mem, m_maxsteps);
    }
    if (m_hmax > 0) {
        IDASetMaxStep(m_ida_mem, m_hmax);
    }
    if (m_maxNonlinIters > 0) {
        int flag = IDASetMaxNonlinIters(m_ida_mem, m_maxNonlinIters);
        checkError(flag, "applyOptions", "IDASetMaxNonlinIters");
    }
    if (m_maxNonlinConvFails > 0) {
        int flag = IDASetMaxConvFails(m_ida_mem, m_maxNonlinConvFails);
        checkError(flag, "applyOptions", "IDASetMaxConvFails");
    }
    if (m_setSuppressAlg != 0) {
        int flag = IDASetSuppressAlg(m_ida_mem, m_setSuppressAlg);
        checkError(flag, "applyOptions", "IDASetSuppressAlg");
    }
}

void IdasIntegrator::sensInit(double t0, FuncEval& func)
{
    m_np = func.nparams();
    m_sens_ok = false;

    N_Vector y = newNVector(static_cast<sd_size_t>(func.neq()), m_sundials_ctx);
    #if CT_SUNDIALS_VERSION >= 60
        m_yS = N_VCloneVectorArray(static_cast<int>(m_np), y);
    #else
        m_yS = N_VCloneVectorArray_Serial(static_cast<int>(m_np), y);
    #endif
    for (size_t n = 0; n < m_np; n++) {
        N_VConst(0.0, m_yS[n]);
    }
    N_VDestroy_Serial(y);
    N_Vector ydot = newNVector(static_cast<sd_size_t>(func.neq()), m_sundials_ctx);
    #if CT_SUNDIALS_VERSION >= 60
        m_ySdot = N_VCloneVectorArray(static_cast<int>(m_np), ydot);
    #else
        m_ySdot = N_VCloneVectorArray_Serial(static_cast<int>(m_np), ydot);
    #endif
    for (size_t n = 0; n < m_np; n++) {
        N_VConst(0.0, m_ySdot[n]);
    }

    int flag = IDASensInit(m_ida_mem, static_cast<sd_size_t>(m_np),
                           IDA_STAGGERED, IDASensResFn(0), m_yS, m_ySdot);
    checkError(flag, "sensInit", "IDASensInit");

    vector_fp atol(m_np);
    for (size_t n = 0; n < m_np; n++) {
        // This scaling factor is tuned so that reaction and species enthalpy
        // sensitivities can be computed simultaneously with the same abstol.
        atol[n] = m_abstolsens / func.m_paramScales[n];
    }
    flag = IDASensSStolerances(m_ida_mem, m_reltolsens, atol.data());
    checkError(flag, "sensInit", "IDASensSStolerances");
}

void IdasIntegrator::integrate(double tout)
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
        throw CanteraError("IdasIntegrator::integrate",
            "IDA error encountered. Error code: {}\n{}\n"
            "{}"
            "Components with largest weighted error estimates:\n{}",
            flag, m_error_message, f_errs, getErrorInfo(10));
    }
}

double IdasIntegrator::step(double tout)
{
    int flag = IDASolve(m_ida_mem, tout, &m_time, m_y, m_ydot, IDA_ONE_STEP);
    if (flag != IDA_SUCCESS) {
        string f_errs = m_func->getErrors();
        if (!f_errs.empty()) {
            f_errs = "Exceptions caught during RHS evaluation:\n" + f_errs;
        }
        throw CanteraError("IdasIntegrator::step",
            "IDA error encountered. Error code: {}\n{}\n"
            "{}"
            "Components with largest weighted error estimates:\n{}",
            flag, f_errs, m_error_message, getErrorInfo(10));

    }
    return m_time;
}

double IdasIntegrator::sensitivity(size_t k, size_t p)
{
    if (m_time == m_t0) {
        // calls to IDAGetSens are only allowed after a successful time step.
        return 0.0;
    }
    if (!m_sens_ok && m_np) {
        int flag = IDAGetSens(m_ida_mem, &m_time, m_yS);
        checkError(flag, "sensitivity", "IDAGetSens");
        m_sens_ok = true;
    }

    if (k >= m_neq) {
        throw CanteraError("IdasIntegrator::sensitivity",
                           "sensitivity: k out of range ({})", k);
    }
    if (p >= m_np) {
        throw CanteraError("IdasIntegrator::sensitivity",
                           "sensitivity: p out of range ({})", p);
    }
    return NV_Ith_S(m_yS[p],k);
}

string IdasIntegrator::getErrorInfo(int N)
{
    N_Vector errs = newNVector(static_cast<sd_size_t>(m_neq), m_sundials_ctx);
    N_Vector errw = newNVector(static_cast<sd_size_t>(m_neq), m_sundials_ctx);
    IDAGetErrWeights(m_ida_mem, errw);
    IDAGetEstLocalErrors(m_ida_mem, errs);

    vector<tuple<double, double, size_t>> weightedErrors;
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
        fmt_append(s, "{}: {}\n", get<2>(weightedErrors[i]), get<1>(weightedErrors[i]));
    }
    return to_string(s);
}

void IdasIntegrator::checkError(long flag, const string& ctMethod,
                                const string& idaMethod) const
{
    if (flag == IDA_SUCCESS) {
        return;
    } else if (flag == IDA_MEM_NULL) {
        throw CanteraError("IdasIntegrator::" + ctMethod,
                           "IDAS integrator is not initialized");
    } else {
        const char* flagname = IDAGetReturnFlagName(flag);
        throw CanteraError("IdasIntegrator::" + ctMethod,
            "{} returned error code {} ({}):\n{}",
            idaMethod, flag, flagname, m_error_message);
    }
}

void IdasIntegrator::setMethod(MethodType t)
{
    if (t != BDF_Method) {
        // IDA only has the BDF method
        throw CanteraError("IdasIntegrator::setMethod", "unknown method");
    }
}

}
