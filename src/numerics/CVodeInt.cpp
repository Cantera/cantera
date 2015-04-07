/**
 *  @file CVodeInt.cpp
 */

// Copyright 2001  California Institute of Technology

#include "CVodeInt.h"
using namespace std;

// cvode includes
#include "../../ext/cvode/include/llnltyps.h"
#include "../../ext/cvode/include/llnlmath.h"
#include "../../ext/cvode/include/cvode.h"
#include "../../ext/cvode/include/cvdense.h"
#include "../../ext/cvode/include/cvdiag.h"
#include "../../ext/cvode/include/cvspgmr.h"
#include "../../ext/cvode/include/cvode.h"

#include "cantera/base/stringUtils.h"

extern "C" {

    /**
     *  Function called by cvode to evaluate ydot given y.  The cvode
     *  integrator allows passing in a void* pointer to access
     *  external data. This pointer is cast to a pointer to a instance
     *  of class FuncEval. The equations to be integrated should be
     *  specified by deriving a class from FuncEval that evaluates the
     *  desired equations.
     *  @ingroup odeGroup
     */
    static void cvode_rhs(integer N, real t, N_Vector y, N_Vector ydot,
                          void* f_data)
    {
        double* ydata = N_VDATA(y);
        double* ydotdata = N_VDATA(ydot);
        Cantera::FuncEval* f = (Cantera::FuncEval*)f_data;
        f->eval(t, ydata, ydotdata, NULL);
    }

    /**
     *  Function called by cvode to evaluate the Jacobian matrix.
     *  (temporary)
     *  @ingroup odeGroup
     */
    static void cvode_jac(integer N, DenseMat J, RhsFn f, void* f_data,
                          real t, N_Vector y, N_Vector fy, N_Vector ewt, real h, real uround,
                          void* jac_data, long int* nfePtr, N_Vector vtemp1, N_Vector vtemp2,
                          N_Vector vtemp3)
    {
        // get pointers to start of data
        double* ydata = N_VDATA(y);
        double* fydata = N_VDATA(fy);
        double* ewtdata = N_VDATA(ewt);
        double* ydot = N_VDATA(vtemp1);

        Cantera::FuncEval* func = (Cantera::FuncEval*)f_data;

        int i,j;
        double* col_j;
        double ysave, dy;
        for (j=0; j < N; j++) {
            col_j = (J->data)[j];
            ysave = ydata[j];
            dy = 1.0/ewtdata[j];
            ydata[j] = ysave + dy;
            dy = ydata[j] - ysave;
            func->eval(t, ydata, ydot, NULL);
            for (i=0; i < N; i++) {
                col_j[i] = (ydot[i] - fydata[i])/dy;
            }
            ydata[j] = ysave;
        }
    }
}

namespace Cantera
{
CVodeInt::CVodeInt() : m_neq(0),
    m_cvode_mem(0),
    m_t0(0.0),
    m_y(0),
    m_abstol(0),
    m_type(DENSE+NOJAC),
    m_itol(0),
    m_method(BDF),
    m_iter(NEWTON),
    m_maxord(0),
    m_reltol(1.e-9),
    m_abstols(1.e-15),
    m_nabs(0),
    m_hmax(0.0),
    m_maxsteps(20000)
{
    m_ropt.resize(OPT_SIZE,0.0);
    m_iopt = new long[OPT_SIZE];
    fill(m_iopt, m_iopt+OPT_SIZE,0);
}

CVodeInt::~CVodeInt()
{
    if (m_cvode_mem) {
        CVodeFree(m_cvode_mem);
    }
    if (m_y) {
        N_VFree(m_y);
    }
    if (m_abstol) {
        N_VFree(m_abstol);
    }
    delete[] m_iopt;
}

double& CVodeInt::solution(size_t k)
{
    return N_VIth(m_y, int(k));
}
double* CVodeInt::solution()
{
    return N_VDATA(m_y);
}

void CVodeInt::setTolerances(double reltol, size_t n, double* abstol)
{
    m_itol = 1;
    m_nabs = int(n);
    if (m_nabs != m_neq) {
        if (m_abstol) {
            N_VFree(m_abstol);
        }
        m_abstol = N_VNew(m_nabs, 0);
    }
    for (int i=0; i<m_nabs; i++) {
        N_VIth(m_abstol, i) = abstol[i];
    }
    m_reltol = reltol;
}

void CVodeInt::setTolerances(double reltol, double abstol)
{
    m_itol = 0;
    m_reltol = reltol;
    m_abstols = abstol;
}

void CVodeInt::setProblemType(int probtype)
{
    m_type = probtype;
}

void CVodeInt::setMethod(MethodType t)
{
    if (t == BDF_Method) {
        m_method = BDF;
    } else if (t == Adams_Method) {
        m_method = ADAMS;
    } else {
        throw CVodeErr("unknown method");
    }
}

void CVodeInt::setMaxStepSize(doublereal hmax)
{
    m_hmax = hmax;
    m_ropt[HMAX] = hmax;
}

void CVodeInt::setMinStepSize(doublereal hmin)
{
    m_hmin = hmin;
    m_ropt[HMIN] = hmin;
}

void CVodeInt::setMaxSteps(int nmax)
{
    m_maxsteps = nmax;
    m_iopt[MXSTEP] = m_maxsteps;
}

void CVodeInt::setIterator(IterType t)
{
    if (t == Newton_Iter) {
        m_iter = NEWTON;
    } else if (t == Functional_Iter) {
        m_iter = FUNCTIONAL;
    } else {
        throw CVodeErr("unknown iterator");
    }
}

void CVodeInt::initialize(double t0, FuncEval& func)
{
    m_neq = int(func.neq());
    m_t0  = t0;

    if (m_y) {
        N_VFree(m_y);    // free solution vector if already allocated
    }
    m_y = N_VNew(m_neq, 0);   // allocate solution vector
    // check abs tolerance array size
    if (m_itol == 1 && m_nabs < m_neq) {
        throw CVodeErr("not enough absolute tolerance values specified.");
    }
    func.getInitialConditions(m_t0, m_neq, N_VDATA(m_y));

    // set options
    m_iopt[MXSTEP] = m_maxsteps;
    m_iopt[MAXORD] = m_maxord;
    m_ropt[HMAX]   = m_hmax;

    if (m_cvode_mem) {
        CVodeFree(m_cvode_mem);
    }

    // pass a pointer to func in m_data
    m_data = (void*)&func;

    if (m_itol) {
        m_cvode_mem = CVodeMalloc(m_neq, cvode_rhs, m_t0, m_y, m_method,
                                  m_iter, m_itol, &m_reltol,
                                  m_abstol, m_data, NULL, 1, m_iopt,
                                  DATA_PTR(m_ropt), NULL);
    } else {
        m_cvode_mem = CVodeMalloc(m_neq, cvode_rhs, m_t0, m_y, m_method,
                                  m_iter, m_itol, &m_reltol,
                                  &m_abstols, m_data, NULL, 1, m_iopt,
                                  DATA_PTR(m_ropt), NULL);
    }

    if (!m_cvode_mem) {
        throw CVodeErr("CVodeMalloc failed.");
    }

    if (m_type == DENSE + NOJAC) {
        CVDense(m_cvode_mem, NULL, NULL);
    } else if (m_type == DENSE + JAC) {
        CVDense(m_cvode_mem, cvode_jac, NULL);
    } else if (m_type == DIAG) {
        CVDiag(m_cvode_mem);
    } else if (m_type == GMRES) {
        CVSpgmr(m_cvode_mem, NONE, MODIFIED_GS, 0, 0.0,
                NULL, NULL, NULL);
    } else {
        throw CVodeErr("unsupported option");
    }
}

void CVodeInt::reinitialize(double t0, FuncEval& func)
{
    m_t0  = t0;
    func.getInitialConditions(m_t0, m_neq, N_VDATA(m_y));

    // set options
    m_iopt[MXSTEP] = m_maxsteps;
    m_iopt[MAXORD] = m_maxord;
    m_ropt[HMAX]   = m_hmax;

    //if (m_cvode_mem) CVodeFree(m_cvode_mem);

    // pass a pointer to func in m_data
    m_data = (void*)&func;
    int result;
    if (m_itol) {
        result = CVReInit(m_cvode_mem, cvode_rhs, m_t0, m_y, m_method,
                          m_iter, m_itol, &m_reltol,
                          m_abstol, m_data, NULL, 1, m_iopt,
                          DATA_PTR(m_ropt), NULL);
    } else {
        result = CVReInit(m_cvode_mem, cvode_rhs, m_t0, m_y, m_method,
                          m_iter, m_itol, &m_reltol,
                          &m_abstols, m_data, NULL, 1, m_iopt,
                          DATA_PTR(m_ropt), NULL);
    }

    if (result != 0) {
        throw CVodeErr("CVReInit failed.");
    }

    if (m_type == DENSE + NOJAC) {
        CVDense(m_cvode_mem, NULL, NULL);
    } else if (m_type == DENSE + JAC) {
        CVDense(m_cvode_mem, cvode_jac, NULL);
    } else if (m_type == DIAG) {
        CVDiag(m_cvode_mem);
    } else if (m_type == GMRES) {
        CVSpgmr(m_cvode_mem, NONE, MODIFIED_GS, 0, 0.0,
                NULL, NULL, NULL);
    } else {
        throw CVodeErr("unsupported option");
    }
}

void CVodeInt::integrate(double tout)
{
    double t;
    int flag;
    flag = CVode(m_cvode_mem, tout, m_y, &t, NORMAL);
    if (flag != SUCCESS) {
        throw CVodeErr(" CVode error encountered. Error code: " + int2str(flag));
    }
}

double CVodeInt::step(double tout)
{
    double t;
    int flag;
    flag = CVode(m_cvode_mem, tout, m_y, &t, ONE_STEP);
    if (flag != SUCCESS) {
        throw CVodeErr(" CVode error encountered. Error code: " + int2str(flag));
    }
    return t;
}

int CVodeInt::nEvals() const
{
    return m_iopt[NFE];
}
}
