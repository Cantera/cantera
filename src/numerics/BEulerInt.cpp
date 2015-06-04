/**
 *  @file BEulerInt.cpp
 */

/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "cantera/numerics/BEulerInt.h"
#include "cantera/numerics/SquareMatrix.h"
#include "cantera/base/global.h"

using namespace std;

namespace Cantera
{

BEulerErr::BEulerErr(const std::string& msg) :
    CanteraError("BEulerInt", msg)
{
}

BEulerInt::BEulerInt() :
    m_iter(Newton_Iter),
    m_method(BEulerVarStep),
    m_jacFormMethod(BEULER_JAC_NUM),
    m_rowScaling(true),
    m_colScaling(false),
    m_matrixConditioning(false),
    m_itol(0),
    m_reltol(1.e-4),
    m_abstols(1.e-10),
    m_hmax(0.0),
    m_maxord(0),
    m_time_step_num(0),
    m_time_step_attempts(0),
    m_max_time_step_attempts(11000000),
    m_numInitialConstantDeltaTSteps(0),
    m_failure_counter(0),
    m_min_newt_its(0),
    m_printSolnStepInterval(1),
    m_printSolnNumberToTout(1),
    m_printSolnFirstSteps(0),
    m_dumpJacobians(false),
    m_neq(0),
    m_t0(0.0),
    m_time_final(0.0),
    time_n(0.0),
    time_nm1(0.0),
    time_nm2(0.0),
    delta_t_n(0.0),
    delta_t_nm1(0.0),
    delta_t_nm2(0.0),
    delta_t_np1(1.0E-8),
    delta_t_max(1.0E300),
    m_func(0),
    tdjac_ptr(0),
    m_print_flag(3),
    m_nfe(0),
    m_nJacEval(0),
    m_numTotalNewtIts(0),
    m_numTotalLinearSolves(0),
    m_numTotalConvFails(0),
    m_numTotalTruncFails(0),
    num_failures(0)
{
    warn_deprecated("class BEulerInt", "To be removed after Cantera 2.2.");
}

BEulerInt::~BEulerInt()
{
    delete tdjac_ptr;
}

void BEulerInt::setTolerances(double reltol, size_t n, double* abstol)
{
    m_itol = 1;
    m_abstol.resize(m_neq);
    if (static_cast<int>(n) != m_neq) {
        printf("ERROR n is wrong\n");
        exit(-1);
    }
    for (int i = 0; i < m_neq; i++) {
        m_abstol[i] = abstol[i];
    }
    m_reltol = reltol;
}

void BEulerInt::setTolerances(double reltol, double abstol)
{
    m_itol = 0;
    m_reltol = reltol;
    m_abstols = abstol;
}

void BEulerInt::setProblemType(int jacFormMethod)
{
    m_jacFormMethod = jacFormMethod;
}

void BEulerInt::setMethodBEMT(BEulerMethodType t)
{
    m_method = t;
}

void BEulerInt::setMaxStep(doublereal hmax)
{
    m_hmax = hmax;
}

void BEulerInt::setMaxNumTimeSteps(int maxNumTimeSteps)
{
    m_max_time_step_attempts = maxNumTimeSteps;
}

void BEulerInt::setNumInitialConstantDeltaTSteps(int num)
{
    m_numInitialConstantDeltaTSteps = num;
}

void BEulerInt::setPrintSolnOptions(int printSolnStepInterval,
                                    int printSolnNumberToTout,
                                    int printSolnFirstSteps,
                                    bool dumpJacobians)
{
    m_printSolnStepInterval = printSolnStepInterval;
    m_printSolnNumberToTout = printSolnNumberToTout;
    m_printSolnFirstSteps   = printSolnFirstSteps;
    m_dumpJacobians         = dumpJacobians;
}

void BEulerInt::setIterator(IterType t)
{
    m_iter = t;
}

void BEulerInt::setNonLinOptions(int min_newt_its, bool matrixConditioning,
                                 bool colScaling, bool rowScaling)
{
    m_min_newt_its = min_newt_its;
    m_matrixConditioning = matrixConditioning;
    m_colScaling = colScaling;
    m_rowScaling = rowScaling;
    if (m_colScaling && m_colScales.empty()) {
        m_colScales.assign(m_neq, 1.0);
    }
    if (m_rowScaling && m_rowScales.empty()) {
        m_rowScales.assign(m_neq, 1.0);
    }
}

void BEulerInt::setInitialTimeStep(double deltaT)
{
    delta_t_np1 = deltaT;
}

void BEulerInt::setPrintFlag(int print_flag)
{
    m_print_flag = print_flag;
}

void BEulerInt::initializeRJE(double t0, ResidJacEval& func)
{
    m_neq = func.nEquations();
    m_t0  = t0;
    internalMalloc();

    /*
     * Get the initial conditions.
     */
    func.getInitialConditions(m_t0, &m_y_n[0], &m_ydot_n[0]);

    // Store a pointer to the residual routine in the object
    m_func = &func;

    /*
     * Initialize the various time counters in the object
     */
    time_n = t0;
    time_nm1 = time_n;
    time_nm2 = time_nm1;
    delta_t_n = 0.0;
    delta_t_nm1 = 0.0;
}

void BEulerInt::reinitializeRJE(double t0, ResidJacEval& func)
{
    m_neq = func.nEquations();
    m_t0  = t0;
    internalMalloc();
    /*
     * At the initial time, get the initial conditions and time and store
     * them into internal storage in the object, my[].
     */
    m_t0  = t0;
    func.getInitialConditions(m_t0, &m_y_n[0], &m_ydot_n[0]);
    /**
     * Set up the internal weights that are used for testing convergence
     */
    setSolnWeights();

    // Store a pointer to the function
    m_func = &func;

}

double BEulerInt::getPrintTime(double time_current)
{
    double tnext;
    if (m_printSolnNumberToTout > 0) {
        double dt = (m_time_final - m_t0) / m_printSolnNumberToTout;
        for (int i = 0; i <= m_printSolnNumberToTout; i++) {
            tnext = m_t0 + dt * i;
            if (tnext >= time_current) {
                return tnext;
            }
        }
    }
    return 1.0E300;
}

int BEulerInt::nEvals() const
{
    return m_nfe;
}

void BEulerInt::internalMalloc()
{
    m_ewt.assign(m_neq, 0.0);
    m_y_n.assign(m_neq, 0.0);
    m_y_nm1.assign(m_neq, 0.0);
    m_y_pred_n.assign(m_neq, 0.0);
    m_ydot_n.assign(m_neq, 0.0);
    m_ydot_nm1.assign(m_neq, 0.0);
    m_resid.assign(m_neq, 0.0);
    m_residWts.assign(m_neq, 0.0);
    m_wksp.assign(m_neq, 0.0);
    if (m_rowScaling) {
        m_rowScales.assign(m_neq, 1.0);
    }
    if (m_colScaling) {
        m_colScales.assign(m_neq, 1.0);
    }
    tdjac_ptr = new SquareMatrix(m_neq);
}

void BEulerInt::setSolnWeights()
{
    int i;
    if (m_itol == 1) {
        for (i = 0; i < m_neq; i++) {
            m_ewt[i] = m_abstol[i] + m_reltol * 0.5 *
                       (fabs(m_y_n[i]) + fabs(m_y_pred_n[i]));
        }
    } else {
        for (i = 0; i < m_neq; i++) {
            m_ewt[i] = m_abstols + m_reltol * 0.5 *
                       (fabs(m_y_n[i]) + fabs(m_y_pred_n[i]));
        }
    }
}

void BEulerInt::setColumnScales()
{
    m_func->calcSolnScales(time_n, &m_y_n[0], &m_y_nm1[0], &m_colScales[0]);
}

void BEulerInt::computeResidWts(GeneralMatrix& jac)
{
    /*
     * We compute residual weights here, which we define as the L_0 norm
     * of the Jacobian Matrix, weighted by the solution weights.
     * This is the proper way to guage the magnitude of residuals. However,
     * it does need the evaluation of the Jacobian, and the implementation
     * below is slow, but doesn't take up much memory.
     *
     * Here a small weighting indicates that the change in solution is
     * very sensitive to that equation.
     */
    int i, j;
    double* data = &(*(jac.begin()));
    double value;
    for (i = 0; i < m_neq; i++) {
        m_residWts[i] = fabs(data[i] * m_ewt[0]);
        for (j = 1; j < m_neq; j++) {
            value = fabs(data[j*m_neq + i] * m_ewt[j]);
            m_residWts[i] = std::max(m_residWts[i], value);
        }
    }
}

double BEulerInt::filterNewStep(double timeCurrent, double* y_current, double* ydot_current)
{
    return 0.0;
}

/*
 * Print out for relevant time step information
 */
static void print_time_step1(int order, int n_time_step, double time,
                             double delta_t_n, double delta_t_nm1,
                             bool step_failed, int num_failures)
{
    const char* string = 0;
    if (order == 0) {
        string = "Backward Euler";
    } else if (order == 1) {
        string = "Forward/Backward Euler";
    } else if (order == 2) {
        string = "Adams-Bashforth/TR";
    }
    writeline('=', 80, true, true);
    printf("\nStart of Time Step: %5d       Time_n = %9.5g Time_nm1 = %9.5g\n",
           n_time_step, time, time - delta_t_n);
    printf("\tIntegration method = %s\n", string);
    if (step_failed) {
        printf("\tPreviously attempted step was a failure\n");
    }
    if (delta_t_n > delta_t_nm1) {
        string = "(Increased from previous iteration)";
    } else if (delta_t_n < delta_t_nm1) {
        string = "(Decreased from previous iteration)";
    } else {
        string = "(same as previous iteration)";
    }
    printf("\tdelta_t_n        = %8.5e %s", delta_t_n, string);
    if (num_failures > 0) {
        printf("\t(Bad_History Failure Counter = %d)", num_failures);
    }
    printf("\n\tdelta_t_nm1      = %8.5e\n", delta_t_nm1);
}

/*
 * Print out for relevant time step information
 */
static void print_time_step2(int  time_step_num, int order,
                             double time, double time_error_factor,
                             double delta_t_n, double delta_t_np1)
{
    printf("\tTime Step Number %5d was a success: time = %10g\n", time_step_num,
           time);
    printf("\t\tEstimated Error\n");
    printf("\t\t--------------------   =   %8.5e\n", time_error_factor);
    printf("\t\tTolerated Error\n\n");
    printf("\t- Recommended next delta_t (not counting history) = %g\n",
           delta_t_np1);
    writeline('=', 80, true, true);
}

/*
 * Print Out descriptive information on why the current step failed
 */
static void print_time_fail(bool convFailure, int time_step_num,
                            double time, double delta_t_n,
                            double delta_t_np1, double  time_error_factor)
{
    writeline('=', 80, true, true);
    if (convFailure) {
        printf("\tTime Step Number %5d experienced a convergence "
               "failure\n", time_step_num);
        printf("\tin the non-linear or linear solver\n");
        printf("\t\tValue of time at failed step           = %g\n", time);
        printf("\t\tdelta_t of the   failed step           = %g\n",
               delta_t_n);
        printf("\t\tSuggested value of delta_t to try next = %g\n",
               delta_t_np1);
    } else {
        printf("\tTime Step Number %5d experienced a truncation error "
               "failure!\n", time_step_num);
        printf("\t\tValue of time at failed step           = %g\n", time);
        printf("\t\tdelta_t of the   failed step           = %g\n",
               delta_t_n);
        printf("\t\tSuggested value of delta_t to try next = %g\n",
               delta_t_np1);
        printf("\t\tCalculated truncation error factor  = %g\n",
               time_error_factor);
    }
    writeline('=', 80, true, true);
}

/*
 * Print out the final results and counters
 */
static void print_final(double time, int step_failed,
                        int time_step_num, int num_newt_its,
                        int total_linear_solves, int numConvFails,
                        int numTruncFails, int nfe, int nJacEval)
{
    writeline('=', 80, true, true);
    printf("TIME INTEGRATION ROUTINE HAS FINISHED: ");
    if (step_failed) {
        printf(" IT WAS A FAILURE\n");
    } else {
        printf(" IT WAS A SUCCESS\n");
    }
    printf("\tEnding time                   = %g\n", time);
    printf("\tNumber of time steps          = %d\n", time_step_num);
    printf("\tNumber of newt its            = %d\n", num_newt_its);
    printf("\tNumber of linear solves       = %d\n", total_linear_solves);
    printf("\tNumber of convergence failures= %d\n", numConvFails);
    printf("\tNumber of TimeTruncErr fails  = %d\n", numTruncFails);
    printf("\tNumber of Function evals      = %d\n", nfe);
    printf("\tNumber of Jacobian evals/solvs= %d\n", nJacEval);
    writeline('=', 80, true, true);
}

/*
 * Header info for one line comment about a time step
 */
static void print_lvl1_Header(int nTimes)
{
    printf("\n");
    if (nTimes) {
        writeline('-', 80);
    }
    printf("time       Time              Time                     Time  ");
    if (nTimes == 0) {
        printf("     START");
    } else {
        printf("    (continued)");
    }
    printf("\n");

    printf("step      (sec)              step  Newt   Aztc bktr  trunc  ");
    printf("\n");

    printf(" No.               Rslt      size    Its  Its  stps  error     |");
    printf("  comment");
    writeline('-', 80, true);
}

/*
 * One line entry about time step
 *   rslt -> 4 letter code
 */
static void print_lvl1_summary(
    int time_step_num, double time, const char* rslt,  double delta_t_n,
    int newt_its, int aztec_its, int bktr_stps, double  time_error_factor,
    const char* comment)
{
    printf("%6d %11.6g %4s %10.4g %4d %4d %4d %11.4g",
           time_step_num, time, rslt, delta_t_n, newt_its, aztec_its,
           bktr_stps, time_error_factor);
    if (comment) {
        printf(" | %s", comment);
    }
    printf("\n");
}

/*
 *   This routine subtracts 2 numbers. If the difference is less
 *   than 1.0E-14 times the magnitude of the smallest number,
 *   then diff returns an exact zero.
 *   It also returns an exact zero if the difference is less than
 *   1.0E-300.
 *
 *   returns:  a - b
 *
 *   This routine is used in numerical differencing schemes in order
 *   to avoid roundoff errors resulting in creating Jacobian terms.
 *   Note: This is a slow routine. However, Jacobian errors may cause
 *         loss of convergence. Therefore, in practice this routine
 *         has proved cost-effective.
 */
double subtractRD(double a, double b)
{
    double diff = a - b;
    double d = std::min(fabs(a), fabs(b));
    d *= 1.0E-14;
    double ad = fabs(diff);
    if (ad < 1.0E-300) {
        diff = 0.0;
    }
    if (ad < d) {
        diff = 0.0;
    }
    return diff;
}

void BEulerInt::beuler_jac(GeneralMatrix& J, double* const f,
                           double time_curr, double CJ,
                           double* const y,
                           double* const ydot,
                           int num_newt_its)
{
    int i, j;
    double* col_j;
    double ysave, ydotsave, dy;
    /**
     * Clear the factor flag
     */
    J.clearFactorFlag();


    if (m_jacFormMethod & BEULER_JAC_ANAL) {
        /********************************************************************
         * Call the function to get a Jacobian.
         */
        m_func->evalJacobian(time_curr, delta_t_n, CJ, y, ydot, J, f);
        m_nJacEval++;
        m_nfe++;
    }  else {
        /*******************************************************************
         * Generic algorithm to calculate a numerical Jacobian
         */
        /*
         * Calculate the current value of the RHS given the
         * current conditions.
         */

        m_func->evalResidNJ(time_curr, delta_t_n, y, ydot, f, JacBase_ResidEval);
        m_nfe++;
        m_nJacEval++;


        /*
         * Malloc a vector and call the function object to return a set of
         * deltaY's that are appropriate for calculating the numerical
         * derivative.
         */
         vector_fp dyVector(m_neq);
        m_func->calcDeltaSolnVariables(time_curr, y, &m_y_nm1[0], &dyVector[0],
                                       &m_ewt[0]);
#ifdef DEBUG_HKM
        bool print_NumJac = false;
        if (print_NumJac) {
            FILE* idy = fopen("NumJac.csv", "w");
            fprintf(idy, "Unk          m_ewt        y     "
                    "dyVector      ResN\n");
            for (int iii = 0; iii < m_neq; iii++) {
                fprintf(idy, " %4d       %16.8e   %16.8e   %16.8e  %16.8e \n",
                        iii,   m_ewt[iii],  y[iii], dyVector[iii], f[iii]);
            }
            fclose(idy);
        }
#endif
        /*
         * Loop over the variables, formulating a numerical derivative
         * of the dense matrix.
         * For the delta in the variable, we will use a variety of approaches
         * The original approach was to use the error tolerance amount.
         * This may not be the best approach, as it could be overly large in
         * some instances and overly small in others.
         * We will first protect from being overly small, by using the usual
         * sqrt of machine precision approach, i.e., 1.0E-7,
         * to bound the lower limit of the delta.
         */
        for (j = 0; j < m_neq; j++) {


            /*
             * Get a pointer into the column of the matrix
             */


            col_j = (double*) J.ptrColumn(j);
            ysave = y[j];
            dy = dyVector[j];

            y[j] = ysave + dy;
            dy = y[j] - ysave;
            ydotsave = ydot[j];
            ydot[j] += dy * CJ;
            /*
             * Call the functon
             */


            m_func->evalResidNJ(time_curr, delta_t_n, y, ydot, &m_wksp[0],
                                JacDelta_ResidEval, j, dy);
            m_nfe++;
            double diff;
            for (i = 0; i < m_neq; i++) {
                diff = subtractRD(m_wksp[i], f[i]);
                col_j[i] = diff / dy;
            }

            y[j] = ysave;
            ydot[j] = ydotsave;

        }
    }


}

void BEulerInt::calc_y_pred(int order)
{
    int i;
    double c1, c2;
    switch (order) {
    case 0:
    case 1:
        c1 = delta_t_n;
        for (i = 0; i < m_neq; i++) {
            m_y_pred_n[i] = m_y_n[i] + c1 * m_ydot_n[i];
        }
        break;
    case 2:
        c1 = delta_t_n * (2.0 + delta_t_n / delta_t_nm1) / 2.0;
        c2 = (delta_t_n * delta_t_n) / (delta_t_nm1 * 2.0);
        for (i = 0; i < m_neq; i++) {
            m_y_pred_n[i] = m_y_n[i] + c1 * m_ydot_n[i] - c2 * m_ydot_nm1[i];
        }
        break;
    }

    /*
     * Filter the predictions.
     */
    m_func->filterSolnPrediction(time_n, &m_y_pred_n[0]);

}

void BEulerInt::calc_ydot(int order, double* y_curr, double* ydot_curr)
{
    int    i;
    double c1;
    switch (order) {
    case 0:
    case 1:             /* First order forward Euler/backward Euler */
        c1 = 1.0 / delta_t_n;
        for (i = 0; i < m_neq; i++) {
            ydot_curr[i] = c1 * (y_curr[i] - m_y_nm1[i]);
        }
        return;
    case 2:             /* Second order Adams-Bashforth / Trapezoidal Rule */
        c1 = 2.0 / delta_t_n;
        for (i = 0; i < m_neq; i++) {
            ydot_curr[i] = c1 * (y_curr[i] - m_y_nm1[i])  - m_ydot_nm1[i];
        }
        return;
    }
}

double BEulerInt::time_error_norm()
{
    int    i;
    double rel_norm, error;
#ifdef DEBUG_HKM
#define NUM_ENTRIES 5
    if (m_print_flag > 2) {
        int imax[NUM_ENTRIES], j, jnum;
        double dmax;
        bool used;
        printf("\t\ttime step truncation error contributors:\n");
        printf("\t\t    I       entry   actual   predicted   "
               "    weight       ydot\n");
        printf("\t\t");
        writeline('-', 70);
        for (j = 0; j < NUM_ENTRIES; j++) {
            imax[j] = -1;
        }
        for (jnum = 0; jnum < NUM_ENTRIES; jnum++) {
            dmax = -1.0;
            for (i = 0; i < m_neq; i++) {
                used = false;
                for (j = 0; j < jnum; j++) {
                    if (imax[j] == i) {
                        used = true;
                    }
                }
                if (!used) {
                    error     = (m_y_n[i] - m_y_pred_n[i]) /  m_ewt[i];
                    rel_norm = sqrt(error * error);
                    if (rel_norm > dmax) {
                        imax[jnum] = i;
                        dmax = rel_norm;
                    }
                }
            }
            if (imax[jnum] >= 0) {
                i = imax[jnum];
                printf("\t\t%4d %12.4e %12.4e %12.4e %12.4e %12.4e\n",
                       i, dmax, m_y_n[i], m_y_pred_n[i], m_ewt[i], m_ydot_n[i]);
            }
        }
        printf("\t\t");
        writeline('-', 70);
    }
#endif
    rel_norm = 0.0;
    for (i = 0; i < m_neq; i++) {
        error     = (m_y_n[i] - m_y_pred_n[i]) /  m_ewt[i];
        rel_norm += (error * error);
    }
    return sqrt(rel_norm / m_neq);
}

double BEulerInt::time_step_control(int order, double time_error_factor)
{
    double factor = 0.0, power = 0.0, delta_t;
    const char*  yo = "time_step_control";

    /*
     * Special case time_error_factor so that zeroes don't cause a problem.
     */
    time_error_factor = std::max(1.0E-50, time_error_factor);

    /*
     * Calculate the factor for the change in magnitude of time step.
     */
    switch (order) {
    case 1:
        factor = 1.0/(2.0 *(time_error_factor));
        power  = 0.5;
        break;
    case 2:
        factor = 1.0/(3.0 * (1.0 + delta_t_nm1 / delta_t_n)
                      * (time_error_factor));
        power  = 0.3333333333333333;
    }
    factor = pow(factor, power);
    if (factor < 0.5) {
        if (m_print_flag > 1) {
            printf("\t%s: WARNING - Current time step will be chucked\n", yo);
            printf("\t\tdue to a time step truncation error failure.\n");
        }
        delta_t = - 0.5 * delta_t_n;
    } else {
        factor  = std::min(factor, 1.5);
        delta_t = factor * delta_t_n;
    }
    return delta_t;
}

double BEulerInt::integrateRJE(double tout, double time_init)
{
    double time_current;
    bool weAreNotFinished = true;
    m_time_final = tout;
    int flag = SUCCESS;
    /**
     * Initialize the time step number to zero. step will increment so that
     * the first time step is number 1
     */
    m_time_step_num = 0;


    /*
     * Do the integration a step at a time
     */
    int istep = 0;
    int printStep = 0;
    bool doPrintSoln = false;
    time_current = time_init;
    time_n = time_init;
    time_nm1 = time_init;
    time_nm2 = time_init;
    m_func->evalTimeTrackingEqns(time_current, 0.0, &m_y_n[0], &m_ydot_n[0]);
    double print_time = getPrintTime(time_current);
    if (print_time == time_current) {
        m_func->writeSolution(4, time_current, delta_t_n,
                              istep, &m_y_n[0], &m_ydot_n[0]);
    }
    /*
     * We print out column headers here for the case of
     */
    if (m_print_flag == 1) {
        print_lvl1_Header(0);
    }
    /*
     * Call a different user routine at the end of each step,
     * that will probably print to a file.
     */
    m_func->user_out2(0, time_current, 0.0, &m_y_n[0], &m_ydot_n[0]);

    do {

        print_time = getPrintTime(time_current);
        if (print_time >= tout) {
            print_time = tout;
        }

        /************************************************************
         * Step the solution
         */
        time_current = step(tout);
        istep++;
        printStep++;
        /***********************************************************/
        if (time_current < 0.0) {
            if (time_current == -1234.) {
                time_current = 0.0;
            } else {
                time_current = -time_current;
            }
            flag = FAILURE;
        }

        if (flag != FAILURE) {
            bool retn =
                m_func->evalStoppingCritera(time_current, delta_t_n,
                                            &m_y_n[0], &m_ydot_n[0]);
            if (retn) {
                weAreNotFinished = false;
                doPrintSoln = true;
            }
        }

        /*
         * determine conditional printing of soln
         */
        if (time_current >= print_time) {
            doPrintSoln = true;
        }
        if (m_printSolnStepInterval == printStep) {
            doPrintSoln = true;
        }
        if (m_printSolnFirstSteps > istep) {
            doPrintSoln = true;
        }

        /*
         * Evaluate time integrated quantities that are calculated at the
         * end of every successful time step.
         */
        if (flag != FAILURE) {
            m_func->evalTimeTrackingEqns(time_current, delta_t_n,
                                         &m_y_n[0], &m_ydot_n[0]);
        }

        /*
         * Call the printout routine.
         */
        if (doPrintSoln) {
            m_func->writeSolution(1, time_current, delta_t_n,
                                  istep, &m_y_n[0], &m_ydot_n[0]);
            printStep = 0;
            doPrintSoln = false;
            if (m_print_flag == 1) {
                print_lvl1_Header(1);
            }
        }
        /*
         * Call a different user routine at the end of each step,
         * that will probably print to a file.
         */
        if (flag == FAILURE) {
            m_func->user_out2(-1, time_current, delta_t_n, &m_y_n[0], &m_ydot_n[0]);
        } else {
            m_func->user_out2(1, time_current, delta_t_n, &m_y_n[0], &m_ydot_n[0]);
        }

    } while (time_current < tout &&
             m_time_step_attempts <  m_max_time_step_attempts &&
             flag == SUCCESS && weAreNotFinished);

    /*
     * Check current time against the max solution time.
     */
    if (time_current >= tout) {
        printf("Simulation completed time integration in %d time steps\n",
               m_time_step_num);
        printf("Final Time: %e\n\n", time_current);
    } else if (m_time_step_attempts >= m_max_time_step_attempts) {
        printf("Simulation ran into time step attempt limit in"
               "%d time steps\n",
               m_time_step_num);
        printf("Final Time: %e\n\n", time_current);
    } else if (flag == FAILURE) {
        printf("ERROR: time stepper failed at time = %g\n", time_current);
    }

    /*
     * Print out the final results and counters.
     */
    print_final(time_n, flag, m_time_step_num, m_numTotalNewtIts,
                m_numTotalLinearSolves, m_numTotalConvFails,
                m_numTotalTruncFails, m_nfe, m_nJacEval);

    /*
     * Call a different user routine at the end of each step,
     * that will probably print to a file.
     */
    m_func->user_out2(2, time_current, delta_t_n, &m_y_n[0], &m_ydot_n[0]);


    if (flag != SUCCESS) {
        throw BEulerErr(" BEuler error encountered.");
    }
    return time_current;
}

double BEulerInt::step(double t_max)
{
    double CJ;
    bool step_failed = false;
    bool giveUp = false;
    bool convFailure = false;
    const char* rslt;
    double time_error_factor = 0.0;
    double normFilter = 0.0;
    int numTSFailures = 0;
    int bktr_stps = 0;
    int nonlinearloglevel = m_print_flag;
    int num_newt_its = 0;
    int aztec_its = 0;
    string comment;
    /*
     * Increment the time counter - May have to be taken back,
     * if time step is found to be faulty.
     */
    m_time_step_num++;

    /**
     * Loop here until we achieve a successful step or we set the giveUp
     * flag indicating that repeated errors have occurred.
     */
    do {
        m_time_step_attempts++;
        comment.clear();

        /*
         * Possibly adjust the delta_t_n value for this time step from the
         * recommended delta_t_np1 value determined in the previous step
         *  due to maximum time step constraints or other occurences,
         * known to happen at a given time.
         */
        if ((time_n + delta_t_np1) >= t_max) {
            delta_t_np1 =t_max - time_n;
        }

        if (delta_t_np1 >= delta_t_max) {
            delta_t_np1 = delta_t_max;
        }

        /*
         * Increment the delta_t counters and the time for the current
         * time step.
         */

        delta_t_nm2 = delta_t_nm1;
        delta_t_nm1 = delta_t_n;
        delta_t_n   = delta_t_np1;
        time_n     += delta_t_n;

        /*
         * Determine the integration order of the current step.
         *
         * Special case for start-up of time integration procedure
         *           First time step = Do a predictor step as we
         *                             have recently added an initial
         *                             ydot input option. And, setting ydot=0
         *                             is equivalent to not doing a
         *                             predictor step.
         *           Second step     = If 2nd order method, do a first order
         *                             step for this time-step, only.
         *
         *           If 2nd order method with a constant time step, the
         *           first and second steps are 1/10 the specified step, and
         *           the third step is 8/10 the specified step.  This reduces
         *           the error asociated with using lower order
         *           integration on the first two steps. (RCS 11-6-97)
         *
         * If the previous time step failed for one reason or another,
         * do a linear step. It's more robust.
         */
        if (m_time_step_num == 1) {
            m_order = 1;                          /* Backward Euler          */
        } else if (m_time_step_num == 2) {
            m_order = 1;                          /* Forward/Backward Euler  */
        } else if (step_failed) {
            m_order = 1;                          /* Forward/Backward Euler  */
        } else if (m_time_step_num > 2) {
            m_order = 1;                          /* Specified
                         Predictor/Corrector
                         - not implemented */
        }

        /*
         * Print out an initial statement about the step.
         */
        if (m_print_flag > 1) {
            print_time_step1(m_order, m_time_step_num, time_n, delta_t_n,
                             delta_t_nm1, step_failed, m_failure_counter);
        }

        /*
         * Calculate the predicted solution, m_y_pred_n, for the current
         * time step.
         */
        calc_y_pred(m_order);

        /*
         * HKM - Commented this out. I may need it for particles later.
         * If Solution bounds checking is turned on, we need to crop the
         * predicted solution to make sure bounds are enforced
         *
         *
         * cropNorm = 0.0;
         * if (Cur_Realm->Realm_Nonlinear.Constraint_Backtracking_Flag ==
         * Constraint_Backtrack_Enable) {
         * cropNorm = cropPredictor(mesh, x_pred_n, abs_time_error,
         *                          m_reltol);
         */

        /*
         * Save the old solution, before overwriting with the new solution
         * - use
         */
        m_y_nm1 = m_y_n;

        /*
         * Use the predicted value as the initial guess for the corrector
         * loop, for
         * every step other than the first step.
         */
        if (m_order > 0) {
            m_y_n = m_y_pred_n;
        }

        /*
         * Save the old time derivative, if necessary, before it is
         * overwritten.
         * This overwrites ydot_nm1, losing information from the previous time
         * step.
         */
        m_ydot_nm1 = m_ydot_n;

        /*
         * Calculate the new time derivative, ydot_n, that is consistent
         * with the
         * initial guess for the corrected solution vector.
         *
         */
        calc_ydot(m_order, &m_y_n[0], &m_ydot_n[0]);

        /*
         * Calculate CJ, the coefficient for the Jacobian corresponding to the
         * derivative of the residual wrt to the acceleration vector.
         */
        if (m_order < 2) {
            CJ = 1.0 / delta_t_n;
        } else {
            CJ = 2.0 / delta_t_n;
        }

        /*
         * Calculate a new Solution Error Weighting vector
         */
        setSolnWeights();

        /*
         * Solve the system of equations at the current time step.
         * Note - x_corr_n and x_dot_n are considered to be updated,
         * on return from this solution.
         */
        int ierror = solve_nonlinear_problem(&m_y_n[0], &m_ydot_n[0],
                                             CJ, time_n, *tdjac_ptr, num_newt_its,
                                             aztec_its, bktr_stps,
                                             nonlinearloglevel);
        /*
         * Set the appropriate flags if a convergence failure is detected.
         */
        if (ierror < 0) {                    /* Step failed */
            convFailure = true;
            step_failed = true;
            rslt = "fail";
            m_numTotalConvFails++;
            m_failure_counter +=3;
            if (m_print_flag > 1) {
                printf("\tStep is Rejected, nonlinear problem didn't converge,"
                       "ierror = %d\n", ierror);
            }
        } else {                             /* Step succeeded */
            convFailure = false;
            step_failed = false;
            rslt = "done";

            /*
             *  Apply a filter to a new successful step
             */
            normFilter = filterNewStep(time_n, &m_y_n[0], &m_ydot_n[0]);
            if (normFilter > 1.0) {
                convFailure = true;
                step_failed = true;
                rslt = "filt";
                if (m_print_flag > 1) {
                    printf("\tStep is Rejected, too large filter adjustment = %g\n",
                           normFilter);
                }
            } else if (normFilter > 0.0) {
                if (normFilter > 0.3) {
                    if (m_print_flag > 1) {
                        printf("\tStep was filtered, norm = %g, next "
                               "time step adjusted\n",  normFilter);
                    }
                } else {
                    if (m_print_flag > 1) {
                        printf("\tStep was filtered, norm = %g\n", normFilter);
                    }
                }
            }
        }

        /*
         * Calculate the time step truncation error for the current step.
         */
        if (!step_failed) {
            time_error_factor = time_error_norm();
        } else {
            time_error_factor = 1000.;
        }

        /*
         * Dynamic time step control- delta_t_n, delta_t_nm1 are set here.
         */
        if (step_failed) {
            /*
             * For convergence failures, decrease the step-size by a factor of
             *  4 and try again.
             */
            delta_t_np1 = 0.25 * delta_t_n;
        } else if (m_method == BEulerVarStep) {

            /*
             * If we are doing a predictor/corrector method, and we are
             * past a certain number of time steps given by the input file
             * then either correct the DeltaT for the next time step or
             *
             */
            if ((m_order > 0) &&
                    (m_time_step_num > m_numInitialConstantDeltaTSteps)) {
                delta_t_np1 = time_step_control(m_order, time_error_factor);
                if (normFilter > 0.1) {
                    if (delta_t_np1 > delta_t_n) {
                        delta_t_np1 = delta_t_n;
                    }
                }

                /*
                 * Check for Current time step failing due to violation of
                 * time step
                 * truncation bounds.
                 */
                if (delta_t_np1 < 0.0) {
                    m_numTotalTruncFails++;
                    step_failed   = true;
                    delta_t_np1   = -delta_t_np1;
                    m_failure_counter += 2;
                    comment += "TIME TRUNC FAILURE";
                    rslt = "TRNC";
                }

                /*
                 * Prevent churning of the time step by not increasing the
                 * time step,
                 * if the recent "History" of the time step behavior is still bad
                 */
                else if (m_failure_counter > 0) {
                    delta_t_np1 = std::min(delta_t_np1, delta_t_n);
                }
            } else {
                delta_t_np1 = delta_t_n;
            }

            /* Decrease time step if a lot of Newton Iterations are
             * taken.
             * The idea being if more or less Newton iteration are taken
             * than the target number of iterations, then adjust the time
             * step downwards so that the target number of iterations or lower
             * is achieved. This
             * should prevent step failure by too many Newton iterations because
             * the time step becomes too large.  CCO
             * hkm -> put in num_new_its min of 3 because the time step
             *        was being altered even when num_newt_its == 1
             */
            int max_Newton_steps = 10000;
            int target_num_iter  = 5;
            if (num_newt_its > 3000 && !step_failed) {
                if (max_Newton_steps != target_num_iter) {
                    double iter_diff        = num_newt_its     - target_num_iter;
                    double iter_adjust_zone = max_Newton_steps - target_num_iter;
                    double target_time_step = delta_t_n
                                              *(1.0 - iter_diff*fabs(iter_diff)/
                                                ((2.0*iter_adjust_zone*iter_adjust_zone)));
                    target_time_step = std::max(0.5*delta_t_n, target_time_step);
                    if (target_time_step < delta_t_np1) {
                        printf("\tNext time step will be decreased from %g to %g"
                               " because of new its restraint\n",
                               delta_t_np1, target_time_step);
                        delta_t_np1 = target_time_step;
                    }
                }
            }


        }

        /*
         * The final loop in the time stepping algorithm depends on whether the
         * current step was a success or not.
         */
        if (step_failed) {
            /*
             * Increment the counter indicating the number of consecutive
             * failures
             */
            numTSFailures++;
            /*
             * Print out a statement about the failure of the time step.
             */
            if (m_print_flag > 1) {
                print_time_fail(convFailure, m_time_step_num, time_n, delta_t_n,
                                delta_t_np1, time_error_factor);
            } else if (m_print_flag == 1) {
                print_lvl1_summary(m_time_step_num, time_n, rslt, delta_t_n,
                                   num_newt_its, aztec_its, bktr_stps,
                                   time_error_factor,
                                   comment.c_str());
            }

            /*
             * Change time step counters back to the previous step before
             * the failed
             * time step occurred.
             */
            time_n     -= delta_t_n;
            delta_t_n   = delta_t_nm1;
            delta_t_nm1 = delta_t_nm2;

            /*
             * Replace old solution vector and time derivative solution vector.
             */
             m_y_n = m_y_nm1;
             m_ydot_n = m_ydot_nm1;
            /*
             * Decide whether to bail on the whole loop
             */
            if (numTSFailures > 35) {
                giveUp = true;
            }
        }

        /*
         * Do processing for a successful step.
         */
        else {

            /*
             * Decrement the number of consequative failure counter.
             */
            m_failure_counter = std::max(0, m_failure_counter-1);

            /*
             * Print out final results of a successfull time step.
             */
            if (m_print_flag > 1) {
                print_time_step2(m_time_step_num, m_order, time_n, time_error_factor,
                                 delta_t_n, delta_t_np1);
            } else if (m_print_flag == 1) {
                print_lvl1_summary(m_time_step_num, time_n, "    ", delta_t_n,
                                   num_newt_its, aztec_its, bktr_stps, time_error_factor,
                                   comment.c_str());
            }

            /*
             * Output information at the end of every successful time step, if
             * requested.
             *
             * fill in
             */


        }
    } while (step_failed && !giveUp);

    /*
     * Send back the overall result of the time step.
     */
    if (step_failed) {
        if (time_n == 0.0) {
            return -1234.0;
        }
        return -time_n;
    }
    return time_n;
}

//-----------------------------------------------------------
//                 Constants
//-----------------------------------------------------------

const double DampFactor = 4;
const int NDAMP = 10;

//-----------------------------------------------------------
//                 MultiNewton methods
//-----------------------------------------------------------

double BEulerInt::soln_error_norm(const double* const delta_y,
                                  bool printLargest)
{
    int    i;
    double sum_norm = 0.0, error;
    for (i = 0; i < m_neq; i++) {
        error     = delta_y[i] / m_ewt[i];
        sum_norm += (error * error);
    }
    sum_norm = sqrt(sum_norm / m_neq);
    if (printLargest) {
        const int num_entries = 8;
        double dmax1, normContrib;
        int j;
        vector_int imax(num_entries, -1);
        printf("\t\tPrintout of Largest Contributors to norm "
               "of value (%g)\n", sum_norm);
        printf("\t\t         I    ysoln  deltaY  weightY  "
               "Error_Norm**2\n");
        printf("\t\t   ");
        writeline('-', 80);
        for (int jnum = 0; jnum < num_entries; jnum++) {
            dmax1 = -1.0;
            for (i = 0; i < m_neq; i++) {
                bool used = false;
                for (j = 0; j < jnum; j++) {
                    if (imax[j] == i) {
                        used = true;
                    }
                }
                if (!used) {
                    error = delta_y[i] / m_ewt[i];
                    normContrib = sqrt(error * error);
                    if (normContrib > dmax1) {
                        imax[jnum] = i;
                        dmax1 = normContrib;
                    }
                }
            }
            i = imax[jnum];
            if (i >= 0) {
                printf("\t\t   %4d %12.4e %12.4e %12.4e %12.4e\n",
                       i, m_y_n[i], delta_y[i], m_ewt[i], dmax1);
            }
        }
        printf("\t\t   ");
        writeline('-', 80);
    }
    return sum_norm;
}
#ifdef DEBUG_HKM_JAC
SquareMatrix jacBack();
#endif

void BEulerInt::doNewtonSolve(double time_curr, double* y_curr,
                              double* ydot_curr, double* delta_y,
                              GeneralMatrix& jac, int loglevel)
{
    int irow, jcol;

    m_func->evalResidNJ(time_curr, delta_t_n, y_curr,
                        ydot_curr, delta_y, Base_ResidEval);
    m_nfe++;
    int sz = m_func->nEquations();
    for (int n = 0; n < sz; n++) {
        delta_y[n] = -delta_y[n];
    }


    /*
     * Column scaling -> We scale the columns of the Jacobian
     * by the nominal important change in the solution vector
     */
    if (m_colScaling) {
        if (!jac.factored()) {
            /*
             * Go get new scales
             */
            setColumnScales();

            /*
             * Scale the new Jacobian
             */
            double* jptr = &(*(jac.begin()));
            for (jcol = 0; jcol < m_neq; jcol++) {
                for (irow = 0; irow < m_neq; irow++) {
                    *jptr *= m_colScales[jcol];
                    jptr++;
                }
            }
        }
    }

    if (m_matrixConditioning) {
        if (jac.factored()) {
            m_func->matrixConditioning(0, sz, delta_y);
        } else {
            double* jptr = &(*(jac.begin()));
            m_func->matrixConditioning(jptr, sz, delta_y);
        }
    }

    /*
     * row sum scaling -> Note, this is an unequivocal success
     *      at keeping the small numbers well balanced and
     *      nonnegative.
     */
    if (m_rowScaling) {
        if (! jac.factored()) {
            /*
             * Ok, this is ugly. jac.begin() returns an vector<double> iterator
             * to the first data location.
             * Then &(*()) reverts it to a double *.
             */
            double* jptr = &(*(jac.begin()));
            for (irow = 0; irow < m_neq; irow++) {
                m_rowScales[irow] = 0.0;
            }
            for (jcol = 0; jcol < m_neq; jcol++) {
                for (irow = 0; irow < m_neq; irow++) {
                    m_rowScales[irow] += fabs(*jptr);
                    jptr++;
                }
            }

            jptr = &(*(jac.begin()));
            for (jcol = 0; jcol < m_neq; jcol++) {
                for (irow = 0; irow < m_neq; irow++) {
                    *jptr /= m_rowScales[irow];
                    jptr++;
                }
            }
        }
        for (irow = 0; irow < m_neq; irow++) {
            delta_y[irow] /= m_rowScales[irow];
        }
    }

#ifdef DEBUG_HKM_JAC
    bool  printJacContributions = false;
    if (m_time_step_num > 304) {
        printJacContributions = false;
    }
    int focusRow = 10;
    int numRows = 2;
    double RRow[2];
    bool freshJac = true;
    RRow[0] = delta_y[focusRow];
    RRow[1] = delta_y[focusRow+1];
    double Pcutoff = 1.0E-70;
    if (!jac.factored()) {
        jacBack = jac;
    } else {
        freshJac = false;
    }
#endif
    /*
     * Solve the system -> This also involves inverting the
     * matrix
     */
    (void) jac.solve(delta_y);


    /*
     * reverse the column scaling if there was any.
     */
    if (m_colScaling) {
        for (irow = 0; irow < m_neq; irow++) {
            delta_y[irow] *= m_colScales[irow];
        }
    }

#ifdef DEBUG_HKM_JAC
    if (printJacContributions) {
        for (int iNum = 0; iNum < numRows; iNum++) {
            if (iNum > 0) {
                focusRow++;
            }
            double dsum = 0.0;
            vector_fp& Jdata = jacBack.data();
            double dRow = Jdata[m_neq * focusRow + focusRow];
            printf("\n Details on delta_Y for row %d \n", focusRow);
            printf("  Value before = %15.5e, delta = %15.5e,"
                   "value after = %15.5e\n", y_curr[focusRow],
                   delta_y[focusRow],
                   y_curr[focusRow] +  delta_y[focusRow]);
            if (!freshJac) {
                printf("    Old Jacobian\n");
            }
            printf("     col          delta_y            aij     "
                   "contrib   \n");
            printf("--------------------------------------------------"
                   "---------------------------------------------\n");
            printf(" Res(%d) %15.5e  %15.5e  %15.5e  (Res = %g)\n",
                   focusRow, delta_y[focusRow],
                   dRow, RRow[iNum] / dRow, RRow[iNum]);
            dsum +=  RRow[iNum] / dRow;
            for (int ii = 0; ii < m_neq; ii++) {
                if (ii != focusRow) {
                    double aij =  Jdata[m_neq * ii + focusRow];
                    double contrib = aij * delta_y[ii] * (-1.0) / dRow;
                    dsum += contrib;
                    if (fabs(contrib) > Pcutoff) {
                        printf("%6d  %15.5e  %15.5e  %15.5e\n", ii,
                               delta_y[ii]  , aij, contrib);
                    }
                }
            }
            printf("--------------------------------------------------"
                   "---------------------------------------------\n");
            printf("        %15.5e                   %15.5e\n",
                   delta_y[focusRow], dsum);
        }
    }

#endif

    m_numTotalLinearSolves++;
}

double BEulerInt::boundStep(const double* const y,
                            const double* const step0, int loglevel)
{
    int i, i_lower = -1, ifbd = 0, i_fbd = 0;
    double fbound = 1.0, f_lowbounds = 1.0, f_delta_bounds = 1.0;
    double ff, y_new, ff_alt;
    for (i = 0; i < m_neq; i++) {
        y_new = y[i] + step0[i];
        if ((y_new < (-0.01 * m_ewt[i])) && y[i] >= 0.0) {
            ff = 0.9 * (y[i] / (y[i] - y_new));
            if (ff < f_lowbounds) {
                f_lowbounds = ff;
                i_lower = i;
            }
        }
        /*
         * Now do a delta bounds
         * Increase variables by a factor of 2 only
         * decrease variables by a factor of 5 only
         */
        ff = 1.0;
        if ((fabs(y_new) > 2.0 * fabs(y[i])) &&
                (fabs(y_new-y[i]) > m_ewt[i])) {
            ff = fabs(y[i]/(y_new - y[i]));
            ff_alt = fabs(m_ewt[i] / (y_new - y[i]));
            ff = std::max(ff, ff_alt);
            ifbd = 1;
        }
        if ((fabs(5.0 * y_new) < fabs(y[i])) &&
                (fabs(y_new - y[i]) > m_ewt[i])) {
            ff = y[i]/(y_new-y[i]) * (1.0 - 5.0)/5.0;
            ff_alt = fabs(m_ewt[i] / (y_new - y[i]));
            ff = std::max(ff, ff_alt);
            ifbd = 0;
        }
        if (ff < f_delta_bounds) {
            f_delta_bounds = ff;
            i_fbd = ifbd;
        }
        f_delta_bounds = std::min(f_delta_bounds, ff);
    }
    fbound = std::min(f_lowbounds, f_delta_bounds);
    /*
     * Report on any corrections
     */
    if (loglevel > 1) {
        if (fbound != 1.0) {
            if (f_lowbounds < f_delta_bounds) {
                printf("\t\tboundStep: Variable %d causing lower bounds "
                       "damping of %g\n",
                       i_lower, f_lowbounds);
            } else {
                if (ifbd) {
                    printf("\t\tboundStep: Decrease of Variable %d causing "
                           "delta damping of %g\n",
                           i_fbd, f_delta_bounds);
                } else {
                    printf("\t\tboundStep: Increase of variable %d causing"
                           "delta damping of %g\n",
                           i_fbd, f_delta_bounds);
                }
            }
        }
    }
    return fbound;
}

int BEulerInt::dampStep(double time_curr, const double* y0,
                        const double* ydot0, const double* step0,
                        double* y1, double* ydot1, double* step1,
                        double& s1, GeneralMatrix& jac,
                        int& loglevel, bool writetitle,
                        int& num_backtracks)
{
    // Compute the weighted norm of the undamped step size step0
    double s0 = soln_error_norm(step0);

    // Compute the multiplier to keep all components in bounds
    // A value of one indicates that there is no limitation
    // on the current step size in the nonlinear method due to
    // bounds constraints (either negative values of delta
    // bounds constraints.
    double fbound = boundStep(y0, step0, loglevel);

    // if fbound is very small, then y0 is already close to the
    // boundary and step0 points out of the allowed domain. In
    // this case, the Newton algorithm fails, so return an error
    // condition.
    if (fbound < 1.e-10) {
        if (loglevel > 1) {
            printf("\t\t\tdampStep: At limits.\n");
        }
        return -3;
    }

    //--------------------------------------------
    //           Attempt damped step
    //--------------------------------------------

    // damping coefficient starts at 1.0
    double damp = 1.0;
    int j, m;
    double ff;
    num_backtracks = 0;
    for (m = 0; m < NDAMP; m++) {

        ff = fbound*damp;

        // step the solution by the damped step size
        /*
         * Whenever we update the solution, we must also always
         * update the time derivative.
         */
        for (j = 0; j < m_neq; j++) {
            y1[j] = y0[j] + ff*step0[j];
        }
        calc_ydot(m_order, y1, ydot1);

        // compute the next undamped step, step1[], that would result
        // if y1[] were accepted.

        doNewtonSolve(time_curr, y1, ydot1, step1, jac, loglevel);

#ifdef DEBUG_HKM
        for (j = 0; j < m_neq; j++) {
            checkFinite(step1[j]);
            checkFinite(y1[j]);
        }
#endif
        // compute the weighted norm of step1
        s1 = soln_error_norm(step1);

        // write log information
        if (loglevel > 3) {
            print_solnDelta_norm_contrib((const double*) step0,
                                         "DeltaSolnTrial",
                                         (const double*) step1,
                                         "DeltaSolnTrialTest",
                                         "dampNewt: Important Entries for "
                                         "Weighted Soln Updates:",
                                         y0, y1, ff, 5);
        }
        if (loglevel > 1) {
            printf("\t\t\tdampNewt: s0 = %g, s1 = %g, fbound = %g,"
                   "damp = %g\n",  s0, s1, fbound, damp);
        }
#ifdef DEBUG_HKM
        if (loglevel > 2) {
            if (s1 > 1.00000001 * s0 && s1 > 1.0E-5) {
                printf("WARNING: Possible Jacobian Problem "
                       "-> turning on more debugging for this step!!!\n");
                print_solnDelta_norm_contrib((const double*) step0,
                                             "DeltaSolnTrial",
                                             (const double*) step1,
                                             "DeltaSolnTrialTest",
                                             "dampNewt: Important Entries for "
                                             "Weighted Soln Updates:",
                                             y0, y1, ff, 5);
                loglevel = 4;
            }
        }
#endif

        // if the norm of s1 is less than the norm of s0, then
        // accept this damping coefficient. Also accept it if this
        // step would result in a converged solution. Otherwise,
        // decrease the damping coefficient and try again.

        if (s1 < 1.0E-5 || s1 < s0) {
            if (loglevel > 2) {
                if (s1 > s0) {
                    if (s1 > 1.0) {
                        printf("\t\t\tdampStep: current trial step and damping"
                               " coefficient accepted because test step < 1\n");
                        printf("\t\t\t          s1 = %g, s0 = %g\n", s1, s0);
                    }
                }
            }
            break;
        } else {
            if (loglevel > 1) {
                printf("\t\t\tdampStep: current step rejected: (s1 = %g > "
                       "s0 = %g)", s1, s0);
                if (m < (NDAMP-1)) {
                    printf(" Decreasing damping factor and retrying");
                } else {
                    printf(" Giving up!!!");
                }
                printf("\n");
            }
        }
        num_backtracks++;
        damp /= DampFactor;
    }

    // If a damping coefficient was found, return 1 if the
    // solution after stepping by the damped step would represent
    // a converged solution, and return 0 otherwise. If no damping
    // coefficient could be found, return -2.
    if (m < NDAMP) {
        if (s1 > 1.0) {
            return 0;
        } else {
            return 1;
        }
    } else {
        if (s1 < 0.5 && (s0 < 0.5)) {
            return 1;
        }
        if (s1 < 1.0) {
            return 0;
        }
        return -2;
    }
}

int BEulerInt::solve_nonlinear_problem(double* const y_comm,
                                       double* const ydot_comm, double CJ,
                                       double time_curr,
                                       GeneralMatrix& jac,
                                       int& num_newt_its,
                                       int& num_linear_solves,
                                       int& num_backtracks,
                                       int loglevel)
{
    int m = 0;
    bool forceNewJac = false;
    double s1=1.e30;

    vector_fp y_curr(y_comm, y_comm + m_neq);
    vector_fp ydot_curr(ydot_comm, ydot_comm + m_neq);
    vector_fp stp(m_neq, 0.0);
    vector_fp stp1(m_neq, 0.0);
    vector_fp y_new(m_neq, 0.0);
    vector_fp ydot_new(m_neq, 0.0);

    bool frst = true;
    num_newt_its = 0;
    num_linear_solves = - m_numTotalLinearSolves;
    num_backtracks = 0;
    int i_backtracks;

    while (1 > 0) {

        /*
         * Increment Newton Solve counter
         */
        m_numTotalNewtIts++;
        num_newt_its++;


        if (loglevel > 1) {
            printf("\t\tSolve_Nonlinear_Problem: iteration %d:\n",
                   num_newt_its);
        }

        // Check whether the Jacobian should be re-evaluated.

        forceNewJac = true;

        if (forceNewJac) {
            if (loglevel > 1) {
                printf("\t\t\tGetting a new Jacobian and solving system\n");
            }
            beuler_jac(jac, &m_resid[0], time_curr, CJ, &y_curr[0], &ydot_curr[0],
                       num_newt_its);
        } else {
            if (loglevel > 1) {
                printf("\t\t\tSolving system with old Jacobian\n");
            }
        }

        // compute the undamped Newton step
        doNewtonSolve(time_curr, &y_curr[0], &ydot_curr[0], &stp[0], jac, loglevel);

        // damp the Newton step
        m = dampStep(time_curr, &y_curr[0], &ydot_curr[0], &stp[0], &y_new[0], &ydot_new[0],
                     &stp1[0], s1, jac, loglevel, frst, i_backtracks);
        frst = false;
        num_backtracks += i_backtracks;

        /*
         * Impose the minimum number of Newton iterations critera
         */
        if (num_newt_its < m_min_newt_its) {
            if (m == 1) {
                m = 0;
            }
        }
        /*
         * Impose max Newton iteration
         */
        if (num_newt_its > 20) {
            m = -1;
            if (loglevel > 1) {
                printf("\t\t\tDampnewton unsuccessful (max newts exceeded) sfinal = %g\n", s1);
            }
        }

        if (loglevel > 1) {
            if (m == 1) {
                printf("\t\t\tDampNewton iteration successful, nonlin "
                       "converged sfinal = %g\n", s1);
            } else if (m == 0) {
                printf("\t\t\tDampNewton iteration successful, get new"
                       "direction, sfinal = %g\n", s1);
            } else {
                printf("\t\t\tDampnewton unsuccessful sfinal = %g\n", s1);
            }
        }

        bool m_filterIntermediate = false;
        if (m_filterIntermediate) {
            if (m == 0) {
                (void) filterNewStep(time_n, &y_new[0], &ydot_new[0]);
            }
        }
        // Exchange new for curr solutions
        if (m == 0 || m == 1) {
            y_curr = y_new;
            calc_ydot(m_order, &y_curr[0], &ydot_curr[0]);
        }

        // convergence
        if (m == 1) {
            goto done;
        }

        // If dampStep fails, first try a new Jacobian if an old
        // one was being used. If it was a new Jacobian, then
        // return -1 to signify failure.
        else if (m < 0) {
            goto done;
        }
    }

done:
    // Copy into the return vectors
    copy(y_curr.begin(), y_curr.end(), y_comm);
    copy(ydot_curr.begin(), ydot_curr.end(), ydot_comm);
    // Increment counters
    num_linear_solves += m_numTotalLinearSolves;

    double time_elapsed = 0.0;
    if (loglevel > 1) {
        if (m == 1) {
            printf("\t\tNonlinear problem solved successfully in "
                   "%d its, time elapsed = %g sec\n",
                   num_newt_its, time_elapsed);
        }
    }
    return m;
}

void BEulerInt::print_solnDelta_norm_contrib(const double* const solnDelta0,
                                             const char* const s0,
                                             const double* const solnDelta1,
                                             const char* const s1,
                                             const char* const title,
                                             const double* const y0,
                                             const double* const y1,
                                             double damp,
                                             int num_entries)
{
    int i, j, jnum;
    bool used;
    double dmax0, dmax1, error, rel_norm;
    printf("\t\t%s currentDamp = %g\n", title, damp);
    printf("\t\t         I  ysoln %10s ysolnTrial "
           "%10s weight relSoln0 relSoln1\n", s0, s1);
    vector_int imax(num_entries, -1);
    printf("\t\t   ");
    writeline('-', 90);
    for (jnum = 0; jnum < num_entries; jnum++) {
        dmax1 = -1.0;
        for (i = 0; i < m_neq; i++) {
            used = false;
            for (j = 0; j < jnum; j++) {
                if (imax[j] == i) {
                    used = true;
                }
            }
            if (!used) {
                error     = solnDelta0[i] /  m_ewt[i];
                rel_norm = sqrt(error * error);
                error     = solnDelta1[i] /  m_ewt[i];
                rel_norm += sqrt(error * error);
                if (rel_norm > dmax1) {
                    imax[jnum] = i;
                    dmax1 = rel_norm;
                }
            }
        }
        if (imax[jnum] >= 0) {
            i = imax[jnum];
            error = solnDelta0[i] /  m_ewt[i];
            dmax0 = sqrt(error * error);
            error = solnDelta1[i] /  m_ewt[i];
            dmax1 = sqrt(error * error);
            printf("\t\t   %4d %12.4e %12.4e %12.4e  %12.4e "
                   "%12.4e %12.4e %12.4e\n",
                   i, y0[i], solnDelta0[i], y1[i],
                   solnDelta1[i], m_ewt[i], dmax0, dmax1);
        }
    }
    printf("\t\t   ");
    writeline('-', 90);
}

} // End of namespace Cantera
