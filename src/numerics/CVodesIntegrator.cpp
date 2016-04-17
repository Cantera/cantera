//! @file CVodesIntegrator.cpp

// Copyright 2001  California Institute of Technology
#include "cantera/numerics/CVodesIntegrator.h"
#include "cantera/base/stringUtils.h"

#include <iostream>
using namespace std;

#include "sundials/sundials_types.h"
#include "sundials/sundials_math.h"
#include "sundials/sundials_nvector.h"
#include "nvector/nvector_serial.h"
#include "cvodes/cvodes.h"
#if SUNDIALS_USE_LAPACK
    #include "cvodes/cvodes_lapack.h"
#else
    #include "cvodes/cvodes_dense.h"
    #include "cvodes/cvodes_band.h"
#endif
#include "cvodes/cvodes_diag.h"
#include "cvodes/cvodes_spgmr.h"

#define CV_SS 1
#define CV_SV 2

#if SUNDIALS_VERSION < 25
typedef int sd_size_t;
#else
typedef long int sd_size_t;
#endif

namespace Cantera
{

class FuncData
{
public:
    FuncData(FuncEval* f, size_t npar = 0) {
        m_pars.resize(npar, 1.0);
        m_func = f;
    }
    virtual ~FuncData() {}
    vector_fp m_pars;
    FuncEval* m_func;
};

extern "C" {
    /**
     * Function called by cvodes to evaluate ydot given y.  The CVODE integrator
     * allows passing in a void* pointer to access external data. This pointer
     * is cast to a pointer to a instance of class FuncEval. The equations to be
     * integrated should be specified by deriving a class from FuncEval that
     * evaluates the desired equations.
     * @ingroup odeGroup
     */
    static int cvodes_rhs(realtype t, N_Vector y, N_Vector ydot,
                          void* f_data)
    {
        try {
            FuncData* d = (FuncData*)f_data;
            FuncEval* f = d->m_func;
            f->eval(t, NV_DATA_S(y), NV_DATA_S(ydot), d->m_pars.data());
        } catch (CanteraError& err) {
            std::cerr << err.what() << std::endl;
            return 1; // possibly recoverable error
        } catch (...) {
            std::cerr << "cvodes_rhs: unhandled exception" << std::endl;
            return -1; // unrecoverable error
        }
        return 0; // successful evaluation
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
}

CVodesIntegrator::CVodesIntegrator() :
    m_neq(0),
    m_cvode_mem(0),
    m_t0(0.0),
    m_y(0),
    m_abstol(0),
    m_type(DENSE+NOJAC),
    m_itol(CV_SS),
    m_method(CV_BDF),
    m_iter(CV_NEWTON),
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
    if (m_y) {
        N_VDestroy_Serial(m_y);
    }
    if (m_abstol) {
        N_VDestroy_Serial(m_abstol);
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
        m_abstol = N_VNew_Serial(static_cast<sd_size_t>(n));
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
    m_type = probtype;
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

void CVodesIntegrator::setMaxErrTestFails(int n)
{
    m_maxErrTestFails = n;
    if (m_cvode_mem) {
        CVodeSetMaxErrTestFails(m_cvode_mem, n);
    }
}

void CVodesIntegrator::setIterator(IterType t)
{
    if (t == Newton_Iter) {
        m_iter = CV_NEWTON;
    } else if (t == Functional_Iter) {
        m_iter = CV_FUNCTIONAL;
    } else {
        throw CanteraError("CVodesIntegrator::setIterator", "unknown iterator");
    }
}

void CVodesIntegrator::sensInit(double t0, FuncEval& func)
{
    m_np = func.nparams();
    m_sens_ok = false;

    N_Vector y = N_VNew_Serial(static_cast<sd_size_t>(func.neq()));
    m_yS = N_VCloneVectorArray_Serial(static_cast<sd_size_t>(m_np), y);
    for (size_t n = 0; n < m_np; n++) {
        N_VConst(0.0, m_yS[n]);
    }

    int flag = CVodeSensInit(m_cvode_mem, static_cast<sd_size_t>(m_np),
                             CV_STAGGERED, CVSensRhsFn(0), m_yS);

    if (flag != CV_SUCCESS) {
        throw CanteraError("CVodesIntegrator::sensInit", "Error in CVodeSensInit");
    }
    vector_fp atol(m_np, m_abstolsens);
    flag = CVodeSensSStolerances(m_cvode_mem, m_reltolsens, atol.data());
}

void CVodesIntegrator::initialize(double t0, FuncEval& func)
{
    m_neq = func.neq();
    m_t0 = t0;
    m_time = t0;

    if (m_y) {
        N_VDestroy_Serial(m_y); // free solution vector if already allocated
    }
    m_y = N_VNew_Serial(static_cast<sd_size_t>(m_neq)); // allocate solution vector
    N_VConst(0.0, m_y);
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
    m_cvode_mem = CVodeCreate(m_method, m_iter);
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

    // pass a pointer to func in m_data
    m_fdata.reset(new FuncData(&func, func.nparams()));

    flag = CVodeSetUserData(m_cvode_mem, m_fdata.get());
    if (flag != CV_SUCCESS) {
        throw CanteraError("CVodesIntegrator::initialize",
                           "CVodeSetUserData failed.");
    }
    if (func.nparams() > 0) {
        sensInit(t0, func);
        flag = CVodeSetSensParams(m_cvode_mem, m_fdata->m_pars.data(),
                                  NULL, NULL);
    }
    applyOptions();
}

void CVodesIntegrator::reinitialize(double t0, FuncEval& func)
{
    m_t0 = t0;
    m_time = t0;
    func.getState(NV_DATA_S(m_y));

    int result = CVodeReInit(m_cvode_mem, m_t0, m_y);
    if (result != CV_SUCCESS) {
        throw CanteraError("CVodesIntegrator::reinitialize",
                           "CVodeReInit failed. result = {}", result);
    }
    applyOptions();
}

void CVodesIntegrator::applyOptions()
{
    if (m_type == DENSE + NOJAC) {
        sd_size_t N = static_cast<sd_size_t>(m_neq);
        #if SUNDIALS_USE_LAPACK
            CVLapackDense(m_cvode_mem, N);
        #else
            CVDense(m_cvode_mem, N);
        #endif
    } else if (m_type == DIAG) {
        CVDiag(m_cvode_mem);
    } else if (m_type == GMRES) {
        CVSpgmr(m_cvode_mem, PREC_NONE, 0);
    } else if (m_type == BAND + NOJAC) {
        sd_size_t N = static_cast<sd_size_t>(m_neq);
        long int nu = m_mupper;
        long int nl = m_mlower;
        #if SUNDIALS_USE_LAPACK
            CVLapackBand(m_cvode_mem, N, nu, nl);
        #else
            CVBand(m_cvode_mem, N, nu, nl);
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
    int flag = CVode(m_cvode_mem, tout, m_y, &m_time, CV_NORMAL);
    if (flag != CV_SUCCESS) {
        throw CanteraError("CVodesIntegrator::integrate",
            "CVodes error encountered. Error code: {}\n{}\n"
            "Components with largest weighted error estimates:\n{}",
            flag, m_error_message, getErrorInfo(10));
    }
    m_sens_ok = false;
}

double CVodesIntegrator::step(double tout)
{
    int flag = CVode(m_cvode_mem, tout, m_y, &m_time, CV_ONE_STEP);
    if (flag != CV_SUCCESS) {
        throw CanteraError("CVodesIntegrator::step",
            "CVodes error encountered. Error code: {}\n{}\n"
            "Components with largest weighted error estimates:\n{}",
            flag, m_error_message, getErrorInfo(10));

    }
    m_sens_ok = false;
    return m_time;
}

int CVodesIntegrator::nEvals() const
{
    long int ne;
    CVodeGetNumRhsEvals(m_cvode_mem, &ne);
    return ne;
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
    N_Vector errs = N_VNew_Serial(static_cast<sd_size_t>(m_neq));
    N_Vector errw = N_VNew_Serial(static_cast<sd_size_t>(m_neq));
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
    fmt::MemoryWriter s;
    for (int i=0; i<N; i++) {
        s.write("{}: {}\n",
                get<2>(weightedErrors[i]), get<1>(weightedErrors[i]));
    }
    return s.str();
}

}
