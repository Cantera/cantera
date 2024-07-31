//! @file CVodesIntegrator.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/CVodesIntegrator.h"
#include "cantera/base/stringUtils.h"

#include <iostream>
using namespace std;

#include "cantera/numerics/sundials_headers.h"

namespace {

N_Vector newNVector(size_t N, Cantera::SundialsContext& context)
{
#if SUNDIALS_VERSION_MAJOR >= 6
    return N_VNew_Serial(static_cast<sd_size_t>(N), context.get());
#else
    return N_VNew_Serial(static_cast<sd_size_t>(N));
#endif
}

} // end anonymous namespace

namespace Cantera
{

extern "C" {
    /**
     * Function called by cvodes to evaluate ydot given y.  The CVODE integrator
     * allows passing in a void* pointer to access external data. This pointer
     * is cast to a pointer to a instance of class FuncEval. The equations to be
     * integrated should be specified by deriving a class from FuncEval that
     * evaluates the desired equations.
     * @ingroup odeGroup
     */
    static int cvodes_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* f_data)
    {
        FuncEval* f = (FuncEval*) f_data;
        return f->evalNoThrow(t, NV_DATA_S(y), NV_DATA_S(ydot));
    }

    //! Function called by CVodes when an error is encountered instead of
    //! writing to stdout. Here, save the error message provided by CVodes so
    //! that it can be included in the subsequently raised CanteraError. Used by
    //! SUNDIALS 6.x and older.
    static void cvodes_err(int error_code, const char* module,
                           const char* function, char* msg, void* eh_data)
    {
        CVodesIntegrator* integrator = (CVodesIntegrator*) eh_data;
        integrator->m_error_message = msg;
        integrator->m_error_message += "\n";
    }

    //! Function called by CVodes when an error is encountered instead of
    //! writing to stdout. Here, save the error message provided by CVodes so
    //! that it can be included in the subsequently raised CanteraError. Used by
    //! SUNDIALS 7.0 and newer.
    #if SUNDIALS_VERSION_MAJOR >= 7
        static void sundials_err(int line, const char *func, const char *file,
                                const char *msg, SUNErrCode err_code,
                                void *err_user_data, SUNContext sunctx)
        {
            CVodesIntegrator* integrator = (CVodesIntegrator*) err_user_data;
            integrator->m_error_message = fmt::format("{}: {}\n", func, msg);
        }
    #endif

    static int cvodes_prec_setup(sunrealtype t, N_Vector y, N_Vector ydot,
                                 sunbooleantype jok, sunbooleantype *jcurPtr,
                                 sunrealtype gamma, void *f_data)
    {
        FuncEval* f = (FuncEval*) f_data;
        if (!jok) {
            *jcurPtr = true; // jacobian data was recomputed
            return f->preconditioner_setup_nothrow(t, NV_DATA_S(y), gamma);
        } else {
            f->updatePreconditioner(gamma); // updates preconditioner with new gamma
            *jcurPtr = false; // indicates that Jacobian data was not recomputed
            return 0; // no error because not recomputed
        }
    }

    static int cvodes_prec_solve(sunrealtype t, N_Vector y, N_Vector ydot, N_Vector r,
                                 N_Vector z, sunrealtype gamma, sunrealtype delta,
                                 int lr, void* f_data)
    {
        FuncEval* f = (FuncEval*) f_data;
        return f->preconditioner_solve_nothrow(NV_DATA_S(r),NV_DATA_S(z));
    }
}

CVodesIntegrator::CVodesIntegrator()
    : m_itol(CV_SS)
    , m_method(CV_BDF)
{
}

CVodesIntegrator::~CVodesIntegrator()
{
    if (m_cvode_mem) {
        if (m_np > 0) {
            CVodeSensFree(m_cvode_mem);
        }
        CVodeFree(&m_cvode_mem);
    }

    SUNLinSolFree((SUNLinearSolver) m_linsol);
    SUNMatDestroy((SUNMatrix) m_linsol_matrix);

    if (m_y) {
        N_VDestroy_Serial(m_y);
    }
    if (m_abstol) {
        N_VDestroy_Serial(m_abstol);
    }
    if (m_dky) {
        N_VDestroy_Serial(m_dky);
    }
    if (m_yS) {
        #if SUNDIALS_VERSION_MAJOR >= 6
            N_VDestroyVectorArray(m_yS, static_cast<int>(m_np));
        #else
            N_VDestroyVectorArray_Serial(m_yS, static_cast<int>(m_np));
        #endif
    }
}


double& CVodesIntegrator::solution(size_t k)
{
    return NV_Ith_S(m_y, k);
}

double* CVodesIntegrator::solution()
{
    return NV_DATA_S(m_y);
}

void CVodesIntegrator::setTolerances(double reltol, size_t n, double* abstol)
{
    m_itol = CV_SV;
    m_nabs = n;
    if (n != m_neq) {
        if (m_abstol) {
            N_VDestroy_Serial(m_abstol);
        }
        m_abstol = newNVector(n, m_sundials_ctx);
    }
    for (size_t i=0; i<n; i++) {
        NV_Ith_S(m_abstol, i) = abstol[i];
    }
    m_reltol = reltol;
}

void CVodesIntegrator::setTolerances(double reltol, double abstol)
{
    m_itol = CV_SS;
    m_reltol = reltol;
    m_abstols = abstol;
}

void CVodesIntegrator::setSensitivityTolerances(double reltol, double abstol)
{
    m_reltolsens = reltol;
    m_abstolsens = abstol;
}

void CVodesIntegrator::setMethod(MethodType t)
{
    if (t == BDF_Method) {
        m_method = CV_BDF;
    } else if (t == Adams_Method) {
        m_method = CV_ADAMS;
    } else {
        throw CanteraError("CVodesIntegrator::setMethod", "unknown method");
    }
}

void CVodesIntegrator::setMaxStepSize(double hmax)
{
    m_hmax = hmax;
    if (m_cvode_mem) {
        CVodeSetMaxStep(m_cvode_mem, hmax);
    }
}

void CVodesIntegrator::setMinStepSize(double hmin)
{
    m_hmin = hmin;
    if (m_cvode_mem) {
        CVodeSetMinStep(m_cvode_mem, hmin);
    }
}

void CVodesIntegrator::setMaxSteps(int nmax)
{
    m_maxsteps = nmax;
    if (m_cvode_mem) {
        CVodeSetMaxNumSteps(m_cvode_mem, m_maxsteps);
    }
}

int CVodesIntegrator::maxSteps()
{
    return m_maxsteps;
}

void CVodesIntegrator::setMaxErrTestFails(int n)
{
    m_maxErrTestFails = n;
    if (m_cvode_mem) {
        CVodeSetMaxErrTestFails(m_cvode_mem, n);
    }
}

void CVodesIntegrator::sensInit(double t0, FuncEval& func)
{
    m_np = func.nparams();
    m_sens_ok = false;

    N_Vector y = newNVector(func.neq(), m_sundials_ctx);
    #if SUNDIALS_VERSION_MAJOR >= 6
        m_yS = N_VCloneVectorArray(static_cast<int>(m_np), y);
    #else
        m_yS = N_VCloneVectorArray_Serial(static_cast<int>(m_np), y);
    #endif
    for (size_t n = 0; n < m_np; n++) {
        N_VConst(0.0, m_yS[n]);
    }
    N_VDestroy_Serial(y);

    int flag = CVodeSensInit(m_cvode_mem, static_cast<int>(m_np),
                             CV_STAGGERED, CVSensRhsFn(0), m_yS);
    checkError(flag, "sensInit", "CVodeSensInit");

    vector<double> atol(m_np);
    for (size_t n = 0; n < m_np; n++) {
        // This scaling factor is tuned so that reaction and species enthalpy
        // sensitivities can be computed simultaneously with the same abstol.
        atol[n] = m_abstolsens / func.m_paramScales[n];
    }
    flag = CVodeSensSStolerances(m_cvode_mem, m_reltolsens, atol.data());
    checkError(flag, "sensInit", "CVodeSensSStolerances");
}

void CVodesIntegrator::initialize(double t0, FuncEval& func)
{
    m_neq = func.neq();
    m_t0 = t0;
    m_time = t0;
    m_tInteg = t0;
    m_func = &func;
    func.clearErrors();
    // Initialize preconditioner if applied
    if (m_prec_side != PreconditionerSide::NO_PRECONDITION) {
        m_preconditioner->initialize(m_neq);
    }
    if (m_y) {
        N_VDestroy_Serial(m_y); // free solution vector if already allocated
    }
    m_y = newNVector(m_neq, m_sundials_ctx); // allocate solution vector
    N_VConst(0.0, m_y);
    if (m_dky) {
        N_VDestroy_Serial(m_dky); // free derivative vector if already allocated
    }
    m_dky = newNVector(m_neq, m_sundials_ctx); // allocate derivative vector
    N_VConst(0.0, m_dky);
    // check abs tolerance array size
    if (m_itol == CV_SV && m_nabs < m_neq) {
        throw CanteraError("CVodesIntegrator::initialize",
                           "not enough absolute tolerance values specified.");
    }

    func.getState(NV_DATA_S(m_y));

    if (m_cvode_mem) {
        CVodeFree(&m_cvode_mem);
    }

    //! Specify the method and the iteration type. Cantera Defaults:
    //!        CV_BDF  - Use BDF methods
    //!        CV_NEWTON - use Newton's method
    #if SUNDIALS_VERSION_MAJOR < 4
        m_cvode_mem = CVodeCreate(m_method, CV_NEWTON);
    #elif SUNDIALS_VERSION_MAJOR < 6
        m_cvode_mem = CVodeCreate(m_method);
    #else
        m_cvode_mem = CVodeCreate(m_method, m_sundials_ctx.get());
    #endif
    if (!m_cvode_mem) {
        throw CanteraError("CVodesIntegrator::initialize",
                           "CVodeCreate failed.");
    }

    int flag = CVodeInit(m_cvode_mem, cvodes_rhs, m_t0, m_y);
    if (flag != CV_SUCCESS) {
        if (flag == CV_MEM_FAIL) {
            throw CanteraError("CVodesIntegrator::initialize",
                               "Memory allocation failed.");
        } else if (flag == CV_ILL_INPUT) {
            throw CanteraError("CVodesIntegrator::initialize",
                               "Illegal value for CVodeInit input argument.");
        } else {
            throw CanteraError("CVodesIntegrator::initialize",
                               "CVodeInit failed.");
        }
    }
    #if SUNDIALS_VERSION_MAJOR >= 7
        flag = SUNContext_PushErrHandler(m_sundials_ctx.get(), &sundials_err, this);
    #else
        flag = CVodeSetErrHandlerFn(m_cvode_mem, &cvodes_err, this);
    #endif

    if (m_itol == CV_SV) {
        flag = CVodeSVtolerances(m_cvode_mem, m_reltol, m_abstol);
        checkError(flag, "initialize", "CVodeSVtolerances");
    } else {
        flag = CVodeSStolerances(m_cvode_mem, m_reltol, m_abstols);
        checkError(flag, "initialize", "CVodeSStolerances");
    }

    flag = CVodeSetUserData(m_cvode_mem, &func);
    checkError(flag, "initialize", "CVodeSetUserData");

    if (func.nparams() > 0) {
        sensInit(t0, func);
        flag = CVodeSetSensParams(m_cvode_mem, func.m_sens_params.data(),
                                  func.m_paramScales.data(), NULL);
        checkError(flag, "initialize", "CVodeSetSensParams");
    }
    applyOptions();
}

void CVodesIntegrator::reinitialize(double t0, FuncEval& func)
{
    m_t0 = t0;
    m_time = t0;
    m_tInteg = t0;
    func.getState(NV_DATA_S(m_y));
    m_func = &func;
    func.clearErrors();
    // reinitialize preconditioner if applied
    if (m_prec_side != PreconditionerSide::NO_PRECONDITION) {
        m_preconditioner->initialize(m_neq);
    }
    int result = CVodeReInit(m_cvode_mem, m_t0, m_y);
    checkError(result, "reinitialize", "CVodeReInit");
    applyOptions();
}

void CVodesIntegrator::applyOptions()
{
    if (m_type == "DENSE") {
        sd_size_t N = static_cast<sd_size_t>(m_neq);
        SUNLinSolFree((SUNLinearSolver) m_linsol);
        SUNMatDestroy((SUNMatrix) m_linsol_matrix);
        #if SUNDIALS_VERSION_MAJOR >= 6
            m_linsol_matrix = SUNDenseMatrix(N, N, m_sundials_ctx.get());
        #else
            m_linsol_matrix = SUNDenseMatrix(N, N);
        #endif
        if (m_linsol_matrix == nullptr) {
            throw CanteraError("CVodesIntegrator::applyOptions",
                "Unable to create SUNDenseMatrix of size {0} x {0}", N);
        }
        int flag;
        #if SUNDIALS_VERSION_MAJOR >= 6
            #if CT_SUNDIALS_USE_LAPACK
                m_linsol = SUNLinSol_LapackDense(m_y, (SUNMatrix) m_linsol_matrix,
                                                    m_sundials_ctx.get());
            #else
                m_linsol = SUNLinSol_Dense(m_y, (SUNMatrix) m_linsol_matrix,
                                            m_sundials_ctx.get());
            #endif
            flag = CVodeSetLinearSolver(m_cvode_mem, (SUNLinearSolver) m_linsol,
                                        (SUNMatrix) m_linsol_matrix);
        #else
            #if CT_SUNDIALS_USE_LAPACK
                m_linsol = SUNLapackDense(m_y, (SUNMatrix) m_linsol_matrix);
            #else
                m_linsol = SUNDenseLinearSolver(m_y, (SUNMatrix) m_linsol_matrix);
            #endif
            flag = CVDlsSetLinearSolver(m_cvode_mem, (SUNLinearSolver) m_linsol,
                                        (SUNMatrix) m_linsol_matrix);
        #endif
        if (m_linsol == nullptr) {
            throw CanteraError("CVodesIntegrator::applyOptions",
                "Error creating Sundials dense linear solver object");
        } else if (flag != CV_SUCCESS) {
            throw CanteraError("CVodesIntegrator::applyOptions",
                "Error connecting linear solver to CVODES. "
                "Sundials error code: {}", flag);
        }

        // throw preconditioner error for DENSE + NOJAC
        if (m_prec_side != PreconditionerSide::NO_PRECONDITION) {
            throw CanteraError("CVodesIntegrator::applyOptions",
                "Preconditioning is not available with the specified problem type.");
        }
    } else if (m_type == "DIAG") {
        CVDiag(m_cvode_mem);
        // throw preconditioner error for DIAG
        if (m_prec_side != PreconditionerSide::NO_PRECONDITION) {
            throw CanteraError("CVodesIntegrator::applyOptions",
                "Preconditioning is not available with the specified problem type.");
        }
    } else if (m_type == "GMRES") {
        SUNLinSolFree((SUNLinearSolver) m_linsol);
        #if SUNDIALS_VERSION_MAJOR >= 6
            m_linsol = SUNLinSol_SPGMR(m_y, SUN_PREC_NONE, 0, m_sundials_ctx.get());
            CVodeSetLinearSolver(m_cvode_mem, (SUNLinearSolver) m_linsol, nullptr);
        #elif SUNDIALS_VERSION_MAJOR >= 4
            m_linsol = SUNLinSol_SPGMR(m_y, PREC_NONE, 0);
            CVSpilsSetLinearSolver(m_cvode_mem, (SUNLinearSolver) m_linsol);
        # else
            m_linsol = SUNSPGMR(m_y, PREC_NONE, 0);
            CVSpilsSetLinearSolver(m_cvode_mem, (SUNLinearSolver) m_linsol);
        #endif
        // set preconditioner if used
        #if SUNDIALS_VERSION_MAJOR >= 4
            if (m_prec_side != PreconditionerSide::NO_PRECONDITION) {
                SUNLinSol_SPGMRSetPrecType((SUNLinearSolver) m_linsol,
                    static_cast<int>(m_prec_side));
                CVodeSetPreconditioner(m_cvode_mem, cvodes_prec_setup,
                    cvodes_prec_solve);
            }
        #else
            if (m_prec_side != PreconditionerSide::NO_PRECONDITION) {
                SUNSPGMRSetPrecType((SUNLinearSolver) m_linsol,
                    static_cast<int>(m_prec_side));
                CVSpilsSetPreconditioner(m_cvode_mem, cvodes_prec_setup,
                    cvodes_prec_solve);
            }
        #endif
    } else if (m_type == "BAND") {
        sd_size_t N = static_cast<sd_size_t>(m_neq);
        sd_size_t nu = m_mupper;
        sd_size_t nl = m_mlower;
        SUNLinSolFree((SUNLinearSolver) m_linsol);
        SUNMatDestroy((SUNMatrix) m_linsol_matrix);
        #if SUNDIALS_VERSION_MAJOR >= 6
            m_linsol_matrix = SUNBandMatrix(N, nu, nl, m_sundials_ctx.get());
        #elif SUNDIALS_VERSION_MAJOR >= 4
            m_linsol_matrix = SUNBandMatrix(N, nu, nl);
        #else
            m_linsol_matrix = SUNBandMatrix(N, nu, nl, nu+nl);
        #endif
        if (m_linsol_matrix == nullptr) {
            throw CanteraError("CVodesIntegrator::applyOptions",
                "Unable to create SUNBandMatrix of size {} with bandwidths "
                "{} and {}", N, nu, nl);
        }
        #if SUNDIALS_VERSION_MAJOR >= 6
            #if CT_SUNDIALS_USE_LAPACK
                m_linsol = SUNLinSol_LapackBand(m_y, (SUNMatrix) m_linsol_matrix,
                                                m_sundials_ctx.get());
            #else
                m_linsol = SUNLinSol_Band(m_y, (SUNMatrix) m_linsol_matrix,
                                            m_sundials_ctx.get());
            #endif
                CVodeSetLinearSolver(m_cvode_mem, (SUNLinearSolver) m_linsol,
                                    (SUNMatrix) m_linsol_matrix);
        #else
            #if CT_SUNDIALS_USE_LAPACK
                m_linsol = SUNLapackBand(m_y, (SUNMatrix) m_linsol_matrix);
            #else
                m_linsol = SUNBandLinearSolver(m_y, (SUNMatrix) m_linsol_matrix);
            #endif
            CVDlsSetLinearSolver(m_cvode_mem, (SUNLinearSolver) m_linsol,
                                (SUNMatrix) m_linsol_matrix);
        #endif
    } else {
        throw CanteraError("CVodesIntegrator::applyOptions",
                           "unsupported linear solver flag '{}'", m_type);
    }

    if (m_maxord > 0) {
        int flag = CVodeSetMaxOrd(m_cvode_mem, m_maxord);
        checkError(flag, "applyOptions", "CVodeSetMaxOrd");
    }
    if (m_maxsteps > 0) {
        CVodeSetMaxNumSteps(m_cvode_mem, m_maxsteps);
    }
    if (m_hmax > 0) {
        CVodeSetMaxStep(m_cvode_mem, m_hmax);
    }
    if (m_hmin > 0) {
        CVodeSetMinStep(m_cvode_mem, m_hmin);
    }
    if (m_maxErrTestFails > 0) {
        CVodeSetMaxErrTestFails(m_cvode_mem, m_maxErrTestFails);
    }
}

void CVodesIntegrator::integrate(double tout)
{
    if (tout == m_time) {
        return;
    } else if (tout < m_time) {
        throw CanteraError("CVodesIntegrator::integrate",
                           "Cannot integrate backwards in time.\n"
                           "Requested time {} < current time {}",
                           tout, m_time);
    }
    int nsteps = 0;
    while (m_tInteg < tout) {
        if (nsteps >= m_maxsteps) {
            string f_errs = m_func->getErrors();
            if (!f_errs.empty()) {
                f_errs = "\nExceptions caught during RHS evaluation:\n" + f_errs;
            }
            throw CanteraError("CVodesIntegrator::integrate",
                "Maximum number of timesteps ({}) taken without reaching output "
                "time ({}).\nCurrent integrator time: {}{}",
                nsteps, tout, m_tInteg, f_errs);
        }
        int flag = CVode(m_cvode_mem, tout, m_y, &m_tInteg, CV_ONE_STEP);
        if (flag != CV_SUCCESS) {
            string f_errs = m_func->getErrors();
            if (!f_errs.empty()) {
                f_errs = "Exceptions caught during RHS evaluation:\n" + f_errs;
            }
            throw CanteraError("CVodesIntegrator::integrate",
                "CVodes error encountered. Error code: {}\n{}\n"
                "{}"
                "Components with largest weighted error estimates:\n{}",
                flag, m_error_message, f_errs, getErrorInfo(10));
        }
        nsteps++;
    }
    int flag = CVodeGetDky(m_cvode_mem, tout, 0, m_y);
    checkError(flag, "integrate", "CVodeGetDky");
    m_time = tout;
    m_sens_ok = false;
}

double CVodesIntegrator::step(double tout)
{
    int flag = CVode(m_cvode_mem, tout, m_y, &m_tInteg, CV_ONE_STEP);
    if (flag != CV_SUCCESS) {
        string f_errs = m_func->getErrors();
        if (!f_errs.empty()) {
            f_errs = "Exceptions caught during RHS evaluation:\n" + f_errs;
        }
        throw CanteraError("CVodesIntegrator::step",
            "CVodes error encountered. Error code: {}\n{}\n"
            "{}"
            "Components with largest weighted error estimates:\n{}",
            flag, f_errs, m_error_message, getErrorInfo(10));

    }
    m_sens_ok = false;
    m_time = m_tInteg;
    return m_time;
}

double* CVodesIntegrator::derivative(double tout, int n)
{
    int flag = CVodeGetDky(m_cvode_mem, tout, n, m_dky);
    checkError(flag, "derivative", "CVodeGetDky");
    return NV_DATA_S(m_dky);
}

int CVodesIntegrator::lastOrder() const
{
    int ord;
    CVodeGetLastOrder(m_cvode_mem, &ord);
    return ord;
}

int CVodesIntegrator::nEvals() const
{
    long int ne;
    CVodeGetNumRhsEvals(m_cvode_mem, &ne);
    return ne;
}

AnyMap CVodesIntegrator::solverStats() const
{
    // AnyMap to return stats
    AnyMap stats;

    // long int linear solver stats provided by CVodes
    long int steps = 0, rhsEvals = 0, errTestFails = 0, jacEvals = 0, linSetup = 0,
             linRhsEvals = 0, linIters = 0, linConvFails = 0, precEvals = 0,
             precSolves = 0, jtSetupEvals = 0, jTimesEvals = 0, nonlinIters = 0,
             nonlinConvFails = 0, orderReductions = 0;
    int lastOrder = 0;
;
    #if SUNDIALS_VERSION_MAJOR >= 4
        CVodeGetNumSteps(m_cvode_mem, &steps);
        CVodeGetNumRhsEvals(m_cvode_mem, &rhsEvals);
        CVodeGetNonlinSolvStats(m_cvode_mem, &nonlinIters, &nonlinConvFails);
        CVodeGetNumErrTestFails(m_cvode_mem, &errTestFails);
        CVodeGetLastOrder(m_cvode_mem, &lastOrder);
        CVodeGetNumStabLimOrderReds(m_cvode_mem, &orderReductions);
        CVodeGetNumJacEvals(m_cvode_mem, &jacEvals);
        CVodeGetNumLinRhsEvals(m_cvode_mem, &linRhsEvals);
        CVodeGetNumLinSolvSetups(m_cvode_mem, &linSetup);
        CVodeGetNumLinIters(m_cvode_mem, &linIters);
        CVodeGetNumLinConvFails(m_cvode_mem, &linConvFails);
        CVodeGetNumPrecEvals(m_cvode_mem, &precEvals);
        CVodeGetNumPrecSolves(m_cvode_mem, &precSolves);
        CVodeGetNumJTSetupEvals(m_cvode_mem, &jtSetupEvals);
        CVodeGetNumJtimesEvals(m_cvode_mem, &jTimesEvals);
    #else
        warn_user("CVodesIntegrator::solverStats", "Function not"
                  "supported with sundials versions less than 4.");
    #endif

    #if SUNDIALS_VERSION_MAJOR >= 6
        long int stepSolveFails = 0;
        CVodeGetNumStepSolveFails(m_cvode_mem, &stepSolveFails);
        stats["step_solve_fails"] = stepSolveFails;
    #endif

    stats["steps"] = steps;
    stats["rhs_evals"] = rhsEvals;
    stats["nonlinear_iters"] = nonlinIters;
    stats["nonlinear_conv_fails"] = nonlinConvFails;
    stats["err_test_fails"] = errTestFails;
    stats["last_order"] = lastOrder;
    stats["stab_order_reductions"] = orderReductions;

    stats["jac_evals"] = jacEvals;
    stats["lin_solve_setups"] = linSetup;
    stats["lin_rhs_evals"] = linRhsEvals;
    stats["lin_iters"] = linIters;
    stats["lin_conv_fails"] = linConvFails;
    stats["prec_evals"] = precEvals;
    stats["prec_solves"] = precSolves;
    stats["jt_vec_setup_evals"] = jtSetupEvals;
    stats["jt_vec_prod_evals"] = jTimesEvals;
    return stats;
}

double CVodesIntegrator::sensitivity(size_t k, size_t p)
{
    if (m_time == m_t0) {
        // calls to CVodeGetSens are only allowed after a successful time step.
        return 0.0;
    }
    if (!m_sens_ok && m_np) {
        int flag = CVodeGetSensDky(m_cvode_mem, m_time, 0, m_yS);
        checkError(flag, "sensitivity", "CVodeGetSens");
        m_sens_ok = true;
    }

    if (k >= m_neq) {
        throw CanteraError("CVodesIntegrator::sensitivity",
                           "sensitivity: k out of range ({})", k);
    }
    if (p >= m_np) {
        throw CanteraError("CVodesIntegrator::sensitivity",
                           "sensitivity: p out of range ({})", p);
    }
    return NV_Ith_S(m_yS[p],k);
}

string CVodesIntegrator::getErrorInfo(int N)
{
    N_Vector errs = newNVector(m_neq, m_sundials_ctx);
    N_Vector errw = newNVector(m_neq, m_sundials_ctx);
    CVodeGetErrWeights(m_cvode_mem, errw);
    CVodeGetEstLocalErrors(m_cvode_mem, errs);

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
        fmt_append(s, "{}: {}\n",
                   get<2>(weightedErrors[i]), get<1>(weightedErrors[i]));
    }
    return to_string(s);
}

void CVodesIntegrator::checkError(long flag, const string& ctMethod,
                                  const string& cvodesMethod) const
{
    if (flag == CV_SUCCESS) {
        return;
    } else if (flag == CV_MEM_NULL) {
        throw CanteraError("CVodesIntegrator::" + ctMethod,
                           "CVODES integrator is not initialized");
    } else {
        const char* flagname = CVodeGetReturnFlagName(flag);
        throw CanteraError("CVodesIntegrator::" + ctMethod,
            "{} returned error code {} ({}):\n{}",
            cvodesMethod, flag, flagname, m_error_message);
    }
}

}
