//! @file CVodesIntegrator.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/CVodesIntegrator.h"
#include "cantera/base/stringUtils.h"

#include <iostream>
using namespace std;

#include "sundials/sundials_types.h"
#include "sundials/sundials_math.h"
#include "sundials/sundials_nvector.h"
#include "nvector/nvector_serial.h"
#include "cvodes/cvodes.h"
#if CT_SUNDIALS_VERSION >= 30
    #if CT_SUNDIALS_USE_LAPACK
        #include "sunlinsol/sunlinsol_lapackdense.h"
        #include "sunlinsol/sunlinsol_lapackband.h"
    #else
        #include "sunlinsol/sunlinsol_dense.h"
        #include "sunlinsol/sunlinsol_band.h"
    #endif
    #include "sunlinsol/sunlinsol_spgmr.h"
    #include "cvodes/cvodes_direct.h"
    #include "cvodes/cvodes_diag.h"
    #include "cvodes/cvodes_spils.h"
#else
    #if CT_SUNDIALS_USE_LAPACK
        #include "cvodes/cvodes_lapack.h"
    #else
        #include "cvodes/cvodes_dense.h"
        #include "cvodes/cvodes_band.h"
    #endif
    #include "cvodes/cvodes_diag.h"
    #include "cvodes/cvodes_spgmr.h"
#endif

#define CV_SS 1
#define CV_SV 2

#if CT_SUNDIALS_VERSION < 25
typedef int sd_size_t;
#else
typedef long int sd_size_t;
#endif

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
    /**
     * Function called by cvodes to evaluate ydot given y.  The CVODE integrator
     * allows passing in a void* pointer to access external data. This pointer
     * is cast to a pointer to a instance of class FuncEval. The equations to be
     * integrated should be specified by deriving a class from FuncEval that
     * evaluates the desired equations.
     * @ingroup odeGroup
     */
    static int cvodes_rhs(realtype t, N_Vector y, N_Vector ydot, void* f_data)
    {
        FuncEval* f = (FuncEval*) f_data;
        return f->eval_nothrow(t, NV_DATA_S(y), NV_DATA_S(ydot));
    }

    //! Function called by CVodes when an error is encountered instead of
    //! writing to stdout. Here, save the error message provided by CVodes so
    //! that it can be included in the subsequently raised CanteraError.
    static void cvodes_err(int error_code, const char* module,
                           const char* function, char* msg, void* eh_data)
    {
        CVodesIntegrator* integrator = (CVodesIntegrator*) eh_data;
        integrator->m_error_message = msg;
        integrator->m_error_message += "\n";
    }

    #if CT_SUNDIALS_VERSION >= 30
        static int cvodes_prec_setup(realtype t, N_Vector y, N_Vector ydot, booleantype
        jok, booleantype *jcurPtr, realtype gamma, void *f_data)
    #else
        static int cvodes_prec_setup(realtype t, N_Vector y, N_Vector ydot, booleantype
        jok, booleantype *jcurPtr, realtype gamma, void *f_data, N_Vector tmp1,
        N_Vector tmp2, N_Vector tmp3)
    #endif
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

    #if CT_SUNDIALS_VERSION >= 30
        static int cvodes_prec_solve(realtype t, N_Vector y, N_Vector ydot, N_Vector r,
        N_Vector z, realtype gamma, realtype delta, int lr, void *f_data)

    #else
        static int cvodes_prec_solve(realtype t, N_Vector y, N_Vector ydot, N_Vector r,
        N_Vector z, realtype gamma, realtype delta, int lr, void *f_data, N_Vector tmp)
    #endif
        {
            FuncEval* f = (FuncEval*) f_data;
            return f->preconditioner_solve_nothrow(NV_DATA_S(r),NV_DATA_S(z));
        }
}

CVodesIntegrator::CVodesIntegrator() :
    m_neq(0),
    m_cvode_mem(0),
    m_linsol(0),
    m_linsol_matrix(0),
    m_func(0),
    m_t0(0.0),
    m_y(0),
    m_abstol(0),
    m_dky(0),
    m_type("DENSE"),
    m_itol(CV_SS),
    m_method(CV_BDF),
    m_maxord(0),
    m_reltol(1.e-9),
    m_abstols(1.e-15),
    m_reltolsens(1.0e-5),
    m_abstolsens(1.0e-4),
    m_nabs(0),
    m_hmax(0.0),
    m_hmin(0.0),
    m_maxsteps(20000),
    m_maxErrTestFails(0),
    m_yS(nullptr),
    m_np(0),
    m_mupper(0), m_mlower(0),
    m_sens_ok(false)
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

    #if CT_SUNDIALS_VERSION >= 30
        SUNLinSolFree((SUNLinearSolver) m_linsol);
        SUNMatDestroy((SUNMatrix) m_linsol_matrix);
    #endif

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
        #if CT_SUNDIALS_VERSION >= 60
            N_VDestroyVectorArray(m_yS, static_cast<sd_size_t>(m_np));
        #else
            N_VDestroyVectorArray_Serial(m_yS, static_cast<sd_size_t>(m_np));
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

void CVodesIntegrator::setProblemType(int probtype)
{
    warn_deprecated("CVodesIntegrator::setProblemType()",
        "To be removed. Set linear solver type with setLinearSolverType");
    if (probtype == DIAG)
    {
        setLinearSolverType("DIAG");
    } else if (probtype == DENSE + NOJAC) {
        setLinearSolverType("DENSE");
    } else if (probtype == BAND + NOJAC) {
        setLinearSolverType("BAND");
    } else if (probtype == GMRES) {
        setLinearSolverType("GMRES");
    } else {
        setLinearSolverType("Invalid Option");
    }
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

void CVodesIntegrator::setMaxStepSize(doublereal hmax)
{
    m_hmax = hmax;
    if (m_cvode_mem) {
        CVodeSetMaxStep(m_cvode_mem, hmax);
    }
}

void CVodesIntegrator::setMinStepSize(doublereal hmin)
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
    #if CT_SUNDIALS_VERSION >= 60
        m_yS = N_VCloneVectorArray(static_cast<sd_size_t>(m_np), y);
    #else
        m_yS = N_VCloneVectorArray_Serial(static_cast<sd_size_t>(m_np), y);
    #endif
    for (size_t n = 0; n < m_np; n++) {
        N_VConst(0.0, m_yS[n]);
    }
    N_VDestroy_Serial(y);

    int flag = CVodeSensInit(m_cvode_mem, static_cast<sd_size_t>(m_np),
                             CV_STAGGERED, CVSensRhsFn(0), m_yS);

    if (flag != CV_SUCCESS) {
        throw CanteraError("CVodesIntegrator::sensInit", "Error in CVodeSensInit");
    }
    vector_fp atol(m_np);
    for (size_t n = 0; n < m_np; n++) {
        // This scaling factor is tuned so that reaction and species enthalpy
        // sensitivities can be computed simultaneously with the same abstol.
        atol[n] = m_abstolsens / func.m_paramScales[n];
    }
    flag = CVodeSensSStolerances(m_cvode_mem, m_reltolsens, atol.data());
}

void CVodesIntegrator::initialize(double t0, FuncEval& func)
{
    m_neq = func.neq();
    m_t0 = t0;
    m_time = t0;
    m_func = &func;
    func.clearErrors();
    // Initialize preconditioner if applied
    if (m_prec_type != PreconditionerType::NO_PRECONDITION) {
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
    #if CT_SUNDIALS_VERSION < 40
        m_cvode_mem = CVodeCreate(m_method, CV_NEWTON);
    #elif CT_SUNDIALS_VERSION < 60
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
    CVodeSetErrHandlerFn(m_cvode_mem, &cvodes_err, this);

    if (m_itol == CV_SV) {
        flag = CVodeSVtolerances(m_cvode_mem, m_reltol, m_abstol);
    } else {
        flag = CVodeSStolerances(m_cvode_mem, m_reltol, m_abstols);
    }
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

    flag = CVodeSetUserData(m_cvode_mem, &func);
    if (flag != CV_SUCCESS) {
        throw CanteraError("CVodesIntegrator::initialize",
                           "CVodeSetUserData failed.");
    }
    if (func.nparams() > 0) {
        sensInit(t0, func);
        flag = CVodeSetSensParams(m_cvode_mem, func.m_sens_params.data(),
                                  func.m_paramScales.data(), NULL);
        if (flag != CV_SUCCESS) {
            throw CanteraError("CVodesIntegrator::initialize",
                               "CVodeSetSensParams failed.");
        }
    }
    applyOptions();
}

void CVodesIntegrator::reinitialize(double t0, FuncEval& func)
{
    m_t0 = t0;
    m_time = t0;
    func.getState(NV_DATA_S(m_y));
    m_func = &func;
    func.clearErrors();
    // reinitialize preconditioner if applied
    if (m_prec_type != PreconditionerType::NO_PRECONDITION) {
        m_preconditioner->initialize(m_neq);
    }
    int result = CVodeReInit(m_cvode_mem, m_t0, m_y);
    if (result != CV_SUCCESS) {
        throw CanteraError("CVodesIntegrator::reinitialize",
                           "CVodeReInit failed. result = {}", result);
    }
    applyOptions();
}

void CVodesIntegrator::applyOptions()
{
    if (m_type == "DENSE") {
        sd_size_t N = static_cast<sd_size_t>(m_neq);
        #if CT_SUNDIALS_VERSION >= 30
            SUNLinSolFree((SUNLinearSolver) m_linsol);
            SUNMatDestroy((SUNMatrix) m_linsol_matrix);
            #if CT_SUNDIALS_VERSION >= 60
                m_linsol_matrix = SUNDenseMatrix(N, N, m_sundials_ctx.get());
            #else
                m_linsol_matrix = SUNDenseMatrix(N, N);
            #endif
            if (m_linsol_matrix == nullptr) {
                throw CanteraError("CVodesIntegrator::applyOptions",
                    "Unable to create SUNDenseMatrix of size {0} x {0}", N);
            }
            int flag;
            #if CT_SUNDIALS_VERSION >= 60
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
        #else
            #if CT_SUNDIALS_USE_LAPACK
                CVLapackDense(m_cvode_mem, N);
            #else
                CVDense(m_cvode_mem, N);
            #endif

        #endif
        // throw preconditioner error for DENSE + NOJAC
        if (m_prec_type != PreconditionerType::NO_PRECONDITION) {
            throw CanteraError("CVodesIntegrator::applyOptions",
                "Preconditioning is not available with the specified problem type.");
        }
    } else if (m_type == "DIAG") {
        CVDiag(m_cvode_mem);
        // throw preconditioner error for DIAG
        if (m_prec_type != PreconditionerType::NO_PRECONDITION) {
            throw CanteraError("CVodesIntegrator::applyOptions",
                "Preconditioning is not available with the specified problem type.");
        }
    } else if (m_type == "GMRES") {
        #if CT_SUNDIALS_VERSION >= 30
            #if CT_SUNDIALS_VERSION >= 60
                m_linsol = SUNLinSol_SPGMR(m_y, PREC_NONE, 0, m_sundials_ctx.get());
                CVodeSetLinearSolver(m_cvode_mem, (SUNLinearSolver) m_linsol, nullptr);
            #elif CT_SUNDIALS_VERSION >= 40
                m_linsol = SUNLinSol_SPGMR(m_y, PREC_NONE, 0);
                CVSpilsSetLinearSolver(m_cvode_mem, (SUNLinearSolver) m_linsol);
            # else
                m_linsol = SUNSPGMR(m_y, PREC_NONE, 0);
                CVSpilsSetLinearSolver(m_cvode_mem, (SUNLinearSolver) m_linsol);
            #endif
            // set preconditioner if used
            #if CT_SUNDIALS_VERSION >= 40
                if (m_prec_type != PreconditionerType::NO_PRECONDITION) {
                    SUNLinSol_SPGMRSetPrecType((SUNLinearSolver) m_linsol,
                        static_cast<int>(m_prec_type));
                    CVodeSetPreconditioner(m_cvode_mem, cvodes_prec_setup,
                        cvodes_prec_solve);
                }
            #else
                if (m_prec_type != PreconditionerType::NO_PRECONDITION) {
                    SUNSPGMRSetPrecType((SUNLinearSolver) m_linsol,
                        static_cast<int>(m_prec_type));
                    CVSpilsSetPreconditioner(m_cvode_mem, cvodes_prec_setup,
                        cvodes_prec_solve);
                }
            #endif
        #else
            if (m_prec_type != PreconditionerType::NO_PRECONDITION) {
                CVSpgmr(m_cvode_mem, static_cast<int>(m_prec_type), 0);
                CVSpilsSetPreconditioner(m_cvode_mem, cvodes_prec_setup,
                    cvodes_prec_solve);
            } else {
                CVSpgmr(m_cvode_mem, PREC_NONE, 0);
            }
        #endif
    } else if (m_type == "BAND") {
        sd_size_t N = static_cast<sd_size_t>(m_neq);
        long int nu = m_mupper;
        long int nl = m_mlower;
        #if CT_SUNDIALS_VERSION >= 30
            SUNLinSolFree((SUNLinearSolver) m_linsol);
            SUNMatDestroy((SUNMatrix) m_linsol_matrix);
            #if CT_SUNDIALS_VERSION >= 60
                m_linsol_matrix = SUNBandMatrix(N, nu, nl, m_sundials_ctx.get());
            #elif CT_SUNDIALS_VERSION >= 40
                m_linsol_matrix = SUNBandMatrix(N, nu, nl);
            #else
                m_linsol_matrix = SUNBandMatrix(N, nu, nl, nu+nl);
            #endif
            if (m_linsol_matrix == nullptr) {
                throw CanteraError("CVodesIntegrator::applyOptions",
                    "Unable to create SUNBandMatrix of size {} with bandwidths "
                    "{} and {}", N, nu, nl);
            }
            #if CT_SUNDIALS_VERSION >= 60
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
        #else
            #if CT_SUNDIALS_USE_LAPACK
                CVLapackBand(m_cvode_mem, N, nu, nl);
            #else
                CVBand(m_cvode_mem, N, nu, nl);
            #endif
        #endif
    } else {
        throw CanteraError("CVodesIntegrator::applyOptions",
                           "unsupported option");
    }

    if (m_maxord > 0) {
        CVodeSetMaxOrd(m_cvode_mem, m_maxord);
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
    }
    int flag = CVode(m_cvode_mem, tout, m_y, &m_time, CV_NORMAL);
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
    m_sens_ok = false;
}

double CVodesIntegrator::step(double tout)
{
    int flag = CVode(m_cvode_mem, tout, m_y, &m_time, CV_ONE_STEP);
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
    return m_time;
}

double* CVodesIntegrator::derivative(double tout, int n)
{
    int flag = CVodeGetDky(m_cvode_mem, tout, n, m_dky);
    if (flag != CV_SUCCESS) {
        string f_errs = m_func->getErrors();
        if (!f_errs.empty()) {
            f_errs = "Exceptions caught evaluating derivative:\n" + f_errs;
        }
        throw CanteraError("CVodesIntegrator::derivative",
            "CVodes error encountered. Error code: {}\n{}\n"
            "{}",
            flag, m_error_message, f_errs);
    }
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

AnyMap CVodesIntegrator::nonlinearSolverStats() const
{
    long int iters = 0, conv_fails = 0;
    // AnyMap to return stats
    AnyMap stats;
    #if CT_SUNDIALS_VERSION >= 40
        CVodeGetNonlinSolvStats(m_cvode_mem, &iters, &conv_fails);
    #else
        warn_user("CVodesIntegrator::nonlinearSolverStats", "Function not"
                  "supported with sundials versions less than 4.");
    #endif
    // Add values to AnyMap
    stats["nonlinear_iters"] = iters;
    stats["nonlinear_conv_fails"] = conv_fails;
    return stats;
}

AnyMap CVodesIntegrator::linearSolverStats() const
{
    // long int linear solver stats provided by CVodes
    long int jacEvals = 0, linRhsEvals = 0, linIters = 0, linConvFails = 0,
             precEvals = 0, precSolves = 0, jtSetupEvals = 0, jTimesEvals = 0;
    // AnyMap to return stats
    AnyMap stats;
    #if CT_SUNDIALS_VERSION >= 60
        CVodeGetLinSolveStats(m_cvode_mem, &jacEvals, &linRhsEvals, &linIters,
            &linConvFails, &precEvals, &precSolves, &jtSetupEvals, &jTimesEvals);
    #elif CT_SUNDIALS_VERSION >= 40
        CVodeGetNumJacEvals(m_cvode_mem, &jacEvals);
        CVodeGetNumLinRhsEvals(m_cvode_mem, &linRhsEvals);
        CVodeGetNumLinIters(m_cvode_mem, &linIters);
        CVodeGetNumLinConvFails(m_cvode_mem, &linConvFails);
        CVodeGetNumPrecEvals(m_cvode_mem, &precEvals);
        CVodeGetNumPrecSolves(m_cvode_mem, &precSolves);
        CVodeGetNumJTSetupEvals(m_cvode_mem, &jtSetupEvals);
        CVodeGetNumJtimesEvals(m_cvode_mem, &jTimesEvals);
    #else
        warn_user("CVodesIntegrator::linearSolverStats", "Function not"
                  "supported with sundials versions less than 4.");
    #endif
    // Add data to AnyMap
    stats["jac_evals"] = jacEvals;
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
        int flag = CVodeGetSens(m_cvode_mem, &m_time, m_yS);
        if (flag != CV_SUCCESS) {
            throw CanteraError("CVodesIntegrator::sensitivity",
                               "CVodeGetSens failed. Error code: {}", flag);
        }
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
        fmt_append(s, "{}: {}\n",
                   get<2>(weightedErrors[i]), get<1>(weightedErrors[i]));
    }
    return to_string(s);
}

}
