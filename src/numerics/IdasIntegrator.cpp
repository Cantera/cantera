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
    return N_VNew_Serial(static_cast<sd_size_t>(N), context.get());
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
* that IDA always calls. The real evaluation of the residual is done by the
* FuncEval::evalDae() method of an instance of a subclass of FuncEval that is passed
* into this function as the `f_data` parameter.
*/
static int ida_rhs(sunrealtype t, N_Vector y, N_Vector ydot, N_Vector r, void* f_data)
{
    FuncEval* f = (FuncEval*) f_data;
    return f->evalDaeNoThrow(t, asSpan(y), asSpan(ydot), asSpan(r));
}

#if SUNDIALS_VERSION_MAJOR >= 7
    //! Function called by CVodes when an error is encountered instead of
    //! writing to stdout. Here, save the error message provided by CVodes so
    //! that it can be included in the subsequently raised CanteraError. Used by
    //! SUNDIALS 7.0 and newer.
    static void sundials_err(int line, const char *func, const char *file,
                            const char *msg, SUNErrCode err_code,
                            void *err_user_data, SUNContext sunctx)
    {
        IdasIntegrator* integrator = (IdasIntegrator*) err_user_data;
        integrator->m_error_message = fmt::format("{}: {}\n", func, msg);
    }
#else
    //! Function called by IDA when an error is encountered instead of writing to
    //! stdout. Here, save the error message provided by IDA so that it can be included
    //! in the subsequently raised CanteraError.
    static void ida_err(int error_code, const char* module,
                        const char* function, char* msg, void* eh_data)
    {
        IdasIntegrator* integrator = (IdasIntegrator*) eh_data;
        integrator->m_error_message = msg;
        integrator->m_error_message += "\n";
    }
#endif

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

span<double> IdasIntegrator::solution()
{
    return asSpan(m_y);
}

void IdasIntegrator::setTolerances(double reltol, span<const double> abstol)
{
    size_t n = abstol.size();
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
    // AnyMap to return stats
    AnyMap stats;

    // long int linear solver stats provided by IDAS
    long int steps = 0, stepSolveFails = 0, resEvals = 0, errTestFails = 0,
             jacEvals = 0, linSetup = 0, linResEvals = 0, linIters = 0,
             linConvFails = 0, precEvals = 0, precSolves = 0, jtSetupEvals = 0,
             jTimesEvals = 0, nonlinIters = 0, nonlinConvFails = 0;
    int lastOrder = 0;

    int flag = IDAGetNumSteps(m_ida_mem, &steps);
    checkError(flag, "solverStats", "IDAGetNumSteps");

    // Remaining stats are best-effort; leave corresponding values at zero if the
    // selected linear solver does not report a given counter.
    IDAGetNumStepSolveFails(m_ida_mem, &stepSolveFails);
    IDAGetNumResEvals(m_ida_mem, &resEvals);
    IDAGetNumNonlinSolvIters(m_ida_mem, &nonlinIters);
    IDAGetNumNonlinSolvConvFails(m_ida_mem, &nonlinConvFails);
    IDAGetNumErrTestFails(m_ida_mem, &errTestFails);
    IDAGetLastOrder(m_ida_mem, &lastOrder);

    IDAGetNumJacEvals(m_ida_mem, &jacEvals);
    IDAGetNumLinResEvals(m_ida_mem, &linResEvals);
    IDAGetNumLinSolvSetups(m_ida_mem, &linSetup);
    IDAGetNumLinIters(m_ida_mem, &linIters);
    IDAGetNumLinConvFails(m_ida_mem, &linConvFails);
    IDAGetNumPrecEvals(m_ida_mem, &precEvals);
    IDAGetNumPrecSolves(m_ida_mem, &precSolves);
    IDAGetNumJTSetupEvals(m_ida_mem, &jtSetupEvals);
    IDAGetNumJtimesEvals(m_ida_mem, &jTimesEvals);

    stats["steps"] = steps;
    stats["step_solve_fails"] = stepSolveFails;
    stats["res_evals"] = resEvals;
    stats["rhs_evals"] = resEvals;
    stats["nonlinear_iters"] = nonlinIters;
    stats["nonlinear_conv_fails"] = nonlinConvFails;
    stats["err_test_fails"] = errTestFails;
    stats["last_order"] = lastOrder;

    stats["jac_evals"] = jacEvals;
    stats["lin_solve_setups"] = linSetup;
    stats["lin_rhs_evals"] = linResEvals;
    stats["lin_iters"] = linIters;
    stats["lin_conv_fails"] = linConvFails;
    stats["prec_evals"] = precEvals;
    stats["prec_solves"] = precSolves;
    stats["jt_vec_setup_evals"] = jtSetupEvals;
    stats["jt_vec_prod_evals"] = jTimesEvals;
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
    m_tInteg = t0;
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
    func.getConstraints(asSpan(m_constraints));

    // get the initial conditions
    func.getStateDae(asSpan(m_y), asSpan(m_ydot));

    if (m_ida_mem) {
        IDAFree(&m_ida_mem);
    }

    //! Create the IDA solver
    m_ida_mem = IDACreate(m_sundials_ctx.get());
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

    #if SUNDIALS_VERSION_MAJOR >= 7
        flag = SUNContext_PushErrHandler(m_sundials_ctx.get(), &sundials_err, this);
    #else
        flag = IDASetErrHandlerFn(m_ida_mem, &ida_err, this);
    #endif

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
    m_tInteg = t0;
    func.getStateDae(asSpan(m_y), asSpan(m_ydot));
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
        m_linsol_matrix = SUNDenseMatrix(N, N, m_sundials_ctx.get());
        #if CT_SUNDIALS_USE_LAPACK
            m_linsol = SUNLinSol_LapackDense(m_y, (SUNMatrix) m_linsol_matrix,
                                             m_sundials_ctx.get());
        #else
            m_linsol = SUNLinSol_Dense(m_y, (SUNMatrix) m_linsol_matrix,
                                        m_sundials_ctx.get());
        #endif
        IDASetLinearSolver(m_ida_mem, (SUNLinearSolver) m_linsol,
                           (SUNMatrix) m_linsol_matrix);
    } else if (m_type == "GMRES") {
        SUNLinSolFree((SUNLinearSolver) m_linsol);
        m_linsol = SUNLinSol_SPGMR(m_y, SUN_PREC_NONE, 0, m_sundials_ctx.get());
        IDASetLinearSolver(m_ida_mem, (SUNLinearSolver) m_linsol, nullptr);
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
    m_yS = N_VCloneVectorArray(static_cast<int>(m_np), y);
    for (size_t n = 0; n < m_np; n++) {
        N_VConst(0.0, m_yS[n]);
    }
    N_VDestroy_Serial(y);
    N_Vector ydot = newNVector(static_cast<sd_size_t>(func.neq()), m_sundials_ctx);
    m_ySdot = N_VCloneVectorArray(static_cast<int>(m_np), ydot);
    for (size_t n = 0; n < m_np; n++) {
        N_VConst(0.0, m_ySdot[n]);
    }

    int flag = IDASensInit(m_ida_mem, static_cast<sd_size_t>(m_np),
                           IDA_STAGGERED, IDASensResFn(0), m_yS, m_ySdot);
    checkError(flag, "sensInit", "IDASensInit");

    vector<double> atol(m_np);
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
    } else if (tout < m_time) {
        throw CanteraError("IdasIntegrator::integrate",
                           "Cannot integrate backwards in time.\n"
                           "Requested time {} < current time {}",
                           tout, m_time);
    }
    int nsteps = 0;
    while (m_tInteg < tout) {
        if (nsteps >= m_maxsteps) {
            throw CanteraError("IdasIntegrator::integrate",
                "Maximum number of timesteps ({}) taken without reaching output "
                "time ({}).\nCurrent integrator time: {}",
                nsteps, tout, m_time);
        }
        int flag = IDASolve(m_ida_mem, tout, &m_tInteg, m_y, m_ydot, IDA_ONE_STEP);
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
        nsteps++;
    }
    int flag = IDAGetDky(m_ida_mem, tout, 0, m_y);
    checkError(flag, "integrate", "IDAGetDky");
    m_time = tout;
}

double IdasIntegrator::step(double tout)
{
    int flag = IDASolve(m_ida_mem, tout, &m_tInteg, m_y, m_ydot, IDA_ONE_STEP);
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
    m_time = m_tInteg;
    return m_time;
}

double IdasIntegrator::sensitivity(size_t k, size_t p)
{
    if (m_time == m_t0) {
        // calls to IDAGetSens are only allowed after a successful time step.
        return 0.0;
    }
    if (!m_sens_ok && m_np) {
        int flag = IDAGetSensDky(m_ida_mem, m_time, 0, m_yS);
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
