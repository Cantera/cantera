/**
 *  @file BEulerInt.h
 */

/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */
#ifndef CT_BEULERINT_H
#define CT_BEULERINT_H

#include "cantera/base/ct_defs.h"
#include "cantera/numerics/ctlapack.h"
#include "cantera/base/utilities.h"
#include "cantera/base/ctexceptions.h"

#include "cantera/numerics/Integrator.h"
#include "cantera/numerics/ResidJacEval.h"

#include "cantera/numerics/GeneralMatrix.h"
#include "cantera/numerics/NonlinearSolver.h"

#define OPT_SIZE 10

#define SUCCESS 0
#define FAILURE 1

#define STEADY 0
#define TRANSIENT 1

namespace Cantera
{

enum BEulerMethodType {
    BEulerFixedStep,
    BEulerVarStep
};

/**
 * Exception class thrown when a BEuler error is encountered.
 */
class BEulerErr : public CanteraError
{
public:
    /**
     * Exception thrown when a BEuler error is encountered. We just call the
     * Cantera Error handler in the initialization list.
     */
    explicit BEulerErr(const std::string& msg);
};


#define BEULER_JAC_ANAL 2
#define BEULER_JAC_NUM  1

/*!
 *  Wrapper class for 'beuler' integrator
 *  We derive the class from the class Integrator
 */
class BEulerInt : public Integrator
{
public:
    /*!
     *  Constructor. Default settings: dense jacobian, no user-supplied
     *  Jacobian function, Newton iteration.
     */
    BEulerInt();

    virtual ~BEulerInt();

    virtual void setTolerances(double reltol, size_t n, double* abstol);
    virtual void setTolerances(double reltol, double abstol);
    virtual void setProblemType(int probtype);

    //! Find the initial conditions for y and ydot.
    virtual void initializeRJE(double t0, ResidJacEval& func);
    virtual void reinitializeRJE(double t0, ResidJacEval& func);
    virtual double integrateRJE(double tout, double tinit = 0.0);

    // This routine advances the calculations one step using a predictor
    // corrector approach. We use an implicit algorithm here.
    virtual doublereal step(double tout);

    //! Set the solution weights. This is a very important routine as it affects
    //! quite a few operations involving convergence.
    virtual void setSolnWeights();

    virtual double& solution(size_t k) {
        return m_y_n[k];
    }
    double* solution() {
        return &m_y_n[0];
    }
    int nEquations() const {
        return m_neq;
    }

    //! Return the total number of function evaluations
    virtual int nEvals() const;
    virtual void setMethodBEMT(BEulerMethodType t);
    virtual void setIterator(IterType t);
    virtual void setMaxStep(double hmax);
    virtual void setMaxNumTimeSteps(int);
    virtual void setNumInitialConstantDeltaTSteps(int);

    void  print_solnDelta_norm_contrib(const double* const soln0,
                                       const char* const s0,
                                       const double* const soln1,
                                       const char* const s1,
                                       const char* const title,
                                       const double* const y0,
                                       const double* const y1,
                                       double damp,
                                       int num_entries);

    //! This routine controls when the solution is printed
    /*!
     * @param printSolnStepInterval If greater than 0, then the soln is
     *                     printed every printSolnStepInterval steps.
     * @param printSolnNumberToTout The solution is printed at regular
     *                  invervals a total of "printSolnNumberToTout" times.
     * @param printSolnFirstSteps The solution is printed out the first
     *                   "printSolnFirstSteps" steps. After these steps the
     *                   other parameters determine the printing. default = 0
     * @param dumpJacobians Dump jacobians to disk.
     */
    virtual void setPrintSolnOptions(int printSolnStepInterval,
                                     int printSolnNumberToTout,
                                     int printSolnFirstSteps = 0,
                                     bool dumpJacobians = false);

    //! Set the options for the nonlinear method
    /*!
     *  Defaults are set in the .h file. These are the defaults:
     *     min_newt_its = 0
     *     matrixConditioning = false
     *     colScaling = false
     *     rowScaling = true
     */
    void setNonLinOptions(int min_newt_its = 0,
                          bool matrixConditioning = false,
                          bool colScaling = false,
                          bool rowScaling = true);
    virtual void setPrintFlag(int print_flag);

    //! Set the column scaling vector at the current time
    virtual void setColumnScales();

    /**
     * Calculate the solution error norm. if printLargest is true, then a table
     * of the largest values is printed to standard output.
     */
    virtual double soln_error_norm(const double* const,
                                   bool printLargest = false);
    virtual void setInitialTimeStep(double delta_t);

    /*!
     *  Function called by to evaluate the Jacobian matrix and the current
     *  residual at the current time step.
     *  @param J = Jacobian matrix to be filled in
     *  @param f = Right hand side. This routine returns the current
     *             value of the rhs (output), so that it does
     *             not have to be computed again.
     */
    void beuler_jac(GeneralMatrix& J, double* const f,
                    double, double, double* const, double* const, int);

protected:
    //!  Internal routine that sets up the fixed length storage based on
    //!  the size of the problem to solve.
    void internalMalloc();

    /*!
     * Function to calculate the predicted solution vector, m_y_pred_n for the
     * (n+1)th time step.  This routine can be used by a first order - forward
     * Euler / backward Euler predictor / corrector method or for a second order
     * Adams-Bashforth / Trapezoidal Rule predictor / corrector method.  See
     * Nachos documentation Sand86-1816 and Gresho, Lee, Sani LLNL report UCRL -
     * 83282 for more information.
     *
     * on input:
     *
     *     N          - number of unknowns
     *     order      - indicates order of method
     *                  = 1 -> first order forward Euler/backward Euler
     *                         predictor/corrector
     *                  = 2 -> second order Adams-Bashforth/Trapezoidal Rule
     *                         predictor/corrector
     *
     *    delta_t_n   - magnitude of time step at time n     (i.e., = t_n+1 - t_n)
     *    delta_t_nm1 - magnitude of time step at time n - 1 (i.e., = t_n - t_n-1)
     *    y_n[]       - solution vector at time n
     *    y_dot_n[]   - acceleration vector from the predictor at time n
     *    y_dot_nm1[] - acceleration vector from the predictor at time n - 1
     *
     * on output:
     *
     *    m_y_pred_n[]    - predicted solution vector at time n + 1
     */
    void calc_y_pred(int);

    /*!
     * Function to calculate the acceleration vector ydot for the first or
     * second order predictor/corrector time integrator.  This routine can be
     * called by a first order - forward Euler / backward Euler predictor /
     * corrector or for a second order Adams - Bashforth / Trapezoidal Rule
     * predictor / corrector.  See Nachos documentation Sand86-1816 and Gresho,
     * Lee, Sani LLNL report UCRL - 83282 for more information.
     *
     *    on input:
     *
     *       N          - number of local unknowns on the processor
     *                    This is equal to internal plus border unknowns.
     *       order      - indicates order of method
     *                    = 1 -> first order forward Euler/backward Euler
     *                           predictor/corrector
     *                    = 2 -> second order Adams-Bashforth/Trapezoidal Rule
     *                           predictor/corrector
     *
     *      delta_t_n   - Magnitude of the current time step at time n
     *                    (i.e., = t_n - t_n-1)
     *      y_curr[]    - Current Solution vector at time n
     *      y_nm1[]     - Solution vector at time n-1
     *      ydot_nm1[] - Acceleration vector at time n-1
     *
     *   on output:
     *
     *      ydot_curr[]   - Current acceleration vector at time n
     *
     * Note we use the current attribute to denote the possibility that
     * y_curr[] may not be equal to m_y_n[] during the nonlinear solve
     * because we may be using a look-ahead scheme.
     */
    void calc_ydot(int, double*, double*);

    /*!
     * Calculates the time step truncation error estimate from a very simple
     * formula based on Gresho et al.  This routine can be called for a first
     * order - forward Euler/backward Euler predictor/ corrector and for a
     * second order Adams- Bashforth/Trapezoidal Rule predictor/corrector. See
     * Nachos documentation Sand86-1816 and Gresho, Lee, LLNL report UCRL -
     * 83282 for more information.
     *
     *    on input:
     *
     *      abs_error   - Generic absolute error tolerance
     *      rel_error   - Generic realtive error tolerance
     *      x_coor[]    - Solution vector from the implicit corrector
     *      x_pred_n[]    - Solution vector from the explicit predictor
     *
     *   on output:
     *
     *      delta_t_n   - Magnitude of next time step at time t_n+1
     *      delta_t_nm1 - Magnitude of previous time step at time t_n
     */
    double time_error_norm();

    /*!
     * Time step control function for the selection of the time step size based on
     * a desired accuracy of time integration and on an estimate of the relative
     * error of the time integration process. This routine can be called for a
     * first order - forward Euler/backward Euler predictor/ corrector and for a
     * second order Adams- Bashforth/Trapezoidal Rule predictor/corrector. See
     * Nachos documentation Sand86-1816 and Gresho, Lee, Sani LLNL report UCRL -
     * 83282 for more information.
     *
     *    on input:
     *
     *       order      - indicates order of method
     *                    = 1 -> first order forward Euler/backward Euler
     *                           predictor/corrector
     *                    = 2 -> second order forward Adams-Bashforth/Trapezoidal
     *                          rule predictor/corrector
     *
     *      delta_t_n   - Magnitude of time step at time t_n
     *      delta_t_nm1 - Magnitude of time step at time t_n-1
     *      rel_error   - Generic realtive error tolerance
     *      time_error_factor   - Estimated value of the time step truncation error
     *                           factor. This value is a ratio of the computed
     *                           error norms. The premultiplying constants
     *                           and the power are not yet applied to normalize the
     *                           predictor/corrector ratio. (see output value)
     *
     *   on output:
     *
     *      return - delta_t for the next time step
     *               If delta_t is negative, then the current time step is
     *               rejected because the time-step truncation error is
     *               too large.  The return value will contain the negative
     *               of the recommended next time step.
     *
     *      time_error_factor  - This output value is normalized so that
     *                           values greater than one indicate the current time
     *                           integration error is greater than the user
     *                           specified magnitude.
     */
    double time_step_control(int m_order, double time_error_factor);

    //! Solve a nonlinear system
    /*!
     * Find the solution to F(X, xprime) = 0 by damped Newton iteration.  On
     * entry, y_comm[] contains an initial estimate of the solution and
     * ydot_comm[] contains an estimate of the derivative.
     *   On  successful return, y_comm[] contains the converged solution
     * and ydot_comm[] contains the derivative
     *
     * @param y_comm[] Contains the input solution. On output y_comm[] contains
     *                 the converged solution
     * @param ydot_comm  Contains the input derivative solution. On output
     *                 y_comm[] contains the converged derivative solution
     * @param CJ       Inverse of the time step
     * @param time_curr  Current value of the time
     * @param jac      Jacobian
     * @param num_newt_its  number of newton iterations
     * @param num_linear_solves number of linear solves
     * @param num_backtracks number of backtracs
     * @param loglevel  Log level
     */
    int solve_nonlinear_problem(double* const y_comm,
                                double* const ydot_comm, double CJ,
                                double time_curr,
                                GeneralMatrix& jac,
                                int& num_newt_its,
                                int& num_linear_solves,
                                int& num_backtracks,
                                int loglevel);

    /**
     * Compute the undamped Newton step. The residual function is
     * evaluated at the current time, t_n, at the current values of the
     * solution vector, m_y_n, and the solution time derivative, m_ydot_n,
     * but the Jacobian is not recomputed.
     */
    void doNewtonSolve(double, double*, double*, double*,
                       GeneralMatrix&, int);


    //!  Bound the Newton step while relaxing the solution
    /*!
     * Return the factor by which the undamped Newton step 'step0'
     * must be multiplied in order to keep all solution components in
     * all domains between their specified lower and upper bounds.
     * Other bounds may be applied here as well.
     *
     * Currently the bounds are hard coded into this routine:
     *
     *  Minimum value for all variables: - 0.01 * m_ewt[i]
     *  Maximum value = none.
     *
     * Thus, this means that all solution components are expected
     * to be numerical greater than zero in the limit of time step
     * truncation errors going to zero.
     *
     * Delta bounds: The idea behind these is that the Jacobian
     *               couldn't possibly be representative if the
     *               variable is changed by a lot. (true for
     *               nonlinear systems, false for linear systems)
     *  Maximum increase in variable in any one newton iteration: factor of 2
     *  Maximum decrease in variable in any one newton iteration: factor of 5
     *
     *   @param y       Current value of the solution
     *   @param step0   Current raw step change in y[]
     *   @param loglevel  Log level. This routine produces output if loglevel
     *                    is greater than one
     *
     *   @return        Returns the damping coefficient
     */
    double boundStep(const double* const y, const double* const step0, int loglevel);

    /*!
     * On entry, step0 must contain an undamped Newton step for the
     * solution x0. This method attempts to find a damping coefficient
     * such that the next undamped step would have a norm smaller than
     * that of step0. If successful, the new solution after taking the
     * damped step is returned in y1, and the undamped step at y1 is
     * returned in step1.
    */
    int dampStep(double, const double*, const double*,
                 const double*, double*, double*,
                 double*, double&, GeneralMatrix&, int&, bool, int&);

    //! Compute Residual Weights
    void computeResidWts(GeneralMatrix& jac);

    //! Filter a new step
    double filterNewStep(double, double*, double*);

    //! Get the next time to print out
    double getPrintTime(double time_current);

    /********************** Member data ***************************/
    /*********************
     * METHOD FLAGS
     *********************/

    //! IterType is used to specify how the nonlinear equations are
    //! to be relaxed at each time step.
    IterType m_iter;
    /**
     * MethodType is used to specify how the time step is to be
     * chosen. Currently, there are two choices, one is a fixed
     * step method while the other is based on a predictor-corrector
     * algorithm and a time-step truncation error tolerance.
     */
    BEulerMethodType m_method;
    /**
     * m_jacFormMethod determines how a matrix is formed.
     */
    int m_jacFormMethod;
    /**
     * m_rowScaling is a boolean. If true then row sum scaling
     * of the Jacobian matrix is carried out when solving the
     * linear systems.
     */
    bool m_rowScaling;
    /**
     * m_colScaling is a boolean. If true, then column scaling
     * is performed on each solution of the linear system.
     */
    bool m_colScaling;
    /**
     * m_matrixConditioning is a boolean. If true, then the
     * Jacobian and every rhs is multiplied by the inverse
     * of a matrix that is suppose to reduce the condition
     * number of the matrix. This is done before row scaling.
     */
    bool m_matrixConditioning;
    /**
     *  If m_itol =1 then each component has an individual
     *  value of atol. If m_itol = 0, the all atols are equal.
     */
    int m_itol;

    //! Relative time truncation error tolerances
    double m_reltol;

    /**
     *  Absolute time truncation error tolerances, when uniform
     *  for all variables.
     */
    double m_abstols;
    /**
     *  Vector of absolute time truncation error tolerance
     *  when not uniform for all variables.
     */
    vector_fp m_abstol;

    //! Error Weights. This is a surprisingly important quantity.
    vector_fp m_ewt;

    //! Maximum step size
    double m_hmax;

    //! Maximum integration order
    int m_maxord;

    //! Current integration order
    int m_order;

    //! Time step number
    int m_time_step_num;
    int m_time_step_attempts;

    //! Max time steps allowed
    int m_max_time_step_attempts;
    /**
     * Number of initial time steps to take where the time truncation error
     * tolerances are not checked. Instead the delta T is uniform
     */
    int m_numInitialConstantDeltaTSteps;

    //! Failure Counter -> keeps track of the number of consecutive failures
    int m_failure_counter;

    //! Minimum Number of Newton Iterations per nonlinear step. default = 0
    int m_min_newt_its;
    /************************
     *  PRINTING OPTIONS
     ************************/
    /**
     * Step Interval at which to print out the solution
     * default = 1;
     * If set to zero, there is no printout
     */
    int m_printSolnStepInterval;
    /**
     * Number of evenly spaced printouts of the solution
     * If zero, there is no printout from this option
     *  default 1
     * If set to zero there is no printout.
     */
    int m_printSolnNumberToTout;

    //! Number of initial steps that the solution is printed out. default = 0
    int m_printSolnFirstSteps;

    //! Dump Jacobians to disk. default false
    bool m_dumpJacobians;

    /*********************
     * INTERNAL SOLUTION VALUES
     *********************/

    //! Number of equations in the ode integrator
    int m_neq;
    vector_fp m_y_n;
    vector_fp m_y_nm1;
    vector_fp m_y_pred_n;
    vector_fp m_ydot_n;
    vector_fp m_ydot_nm1;
    /************************
     * TIME VARIABLES
     ************************/

    //! Initial time at the start of the integration
    double m_t0;

    //! Final time
    double m_time_final;

    double time_n;
    double time_nm1;
    double time_nm2;
    double delta_t_n;
    double delta_t_nm1;
    double delta_t_nm2;
    double delta_t_np1;

    //! Maximum permissible time step
    double delta_t_max;

    vector_fp m_resid;
    vector_fp m_residWts;
    vector_fp m_wksp;
    ResidJacEval* m_func;
    vector_fp m_rowScales;
    vector_fp m_colScales;

    //! Pointer to the jacobian representing the time dependent problem
    GeneralMatrix* tdjac_ptr;

    /**
     * Determines the level of printing for each time
     * step.
     *   0 -> absolutely nothing is printed for
     *        a single time step.
     *   1 -> One line summary per time step
     *   2 -> short description, points of interest
     *   3 -> Lots printed per time step (default)
     */
    int m_print_flag;
    /***************************************************************************
     * COUNTERS OF VARIOUS KINDS
     ***************************************************************************/

    //! Number of function evaluations
    int m_nfe;

    /**
     * Number of Jacobian Evaluations and
     * factorization steps (they are the same)
     */
    int m_nJacEval;

    //! Number of total newton iterations
    int m_numTotalNewtIts;

    //! Total number of linear iterations
    int m_numTotalLinearSolves;

    //! Total number of convergence failures.
    int m_numTotalConvFails;

    //! Total Number of time truncation error failures
    int m_numTotalTruncFails;

    int num_failures;
};

}    // namespace

#endif // CT_BEULER
