/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2011/08/01 21:20:16 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * This simple example problem for IDA, due to Robertson, 
 * is from chemical kinetics, and consists of the following three 
 * equations:
 *
 *      dy1/dt = -.04*y1 + 1.e4*y2*y3
 *      dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
 *         0   = y1 + y2 + y3 - 1
 *
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1, y2 = y3 = 0.
 *
 * While integrating the system, we also use the rootfinding
 * feature to find the points at which y1 = 1e-4 or at which
 * y3 = 0.01.
 *
 * The problem is solved with IDA using IDADENSE for the linear
 * solver, with a user-supplied Jacobian. Output is printed at
 * t = .4, 4, 40, ..., 4e10.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <math.h>

#include <cantera/Cantera.h>
#include <cantera/numerics.h>
#include <cantera/kernel/ResidJacEval.h>    // defines class Water
#include <cantera/kernel/NonlinearSolver.h>
#include "cantera/kernel/DAE_Solver.h"
#include "cantera/kernel/IDA_Solver.h"

#include <ida/ida.h>
#include <ida/ida_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

using namespace Cantera;
using namespace std;

/* Problem Constants */

#define NEQ   3
#define NOUT  12

#define ZERO RCONST(0.0);
#define ONE  RCONST(1.0);

/* Macro to define dense matrix elements, indexed from 1. */

#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

/* Prototypes of functions called by IDA */

int resrob(realtype tres, N_Vector yy, N_Vector yp, 
           N_Vector resval, void *user_data);

static int grob(realtype t, N_Vector yy, N_Vector yp,
                realtype *gout, void *user_data);

int jacrob(int Neq, realtype tt,  realtype cj, 
           N_Vector yy, N_Vector yp, N_Vector resvec,
           DlsMat JJ, void *user_data,
           N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);

/* Prototypes of private functions */
static void PrintHeader(realtype rtol, N_Vector avtol, N_Vector y);
static void PrintOutput(void *mem, realtype t, N_Vector y);
static void PrintRootInfo(int root_f1, int root_f2);
static void PrintFinalStats(void *mem);
static int check_flag(void *flagvalue, const char *funcname, int opt);


class RobertsFunc : public Cantera::ResidJacEval {
public:
  RobertsFunc(doublereal atol = 1.0e-13) :
    ResidJacEval(atol)
  {
    neq_ = 3;
  }
  virtual int  evalResidNJ(doublereal t, const doublereal deltaT,
                           const doublereal * const y,
                           const doublereal * const ydot,
                           doublereal * const resid,
                           ResidEval_Type_Enum evalType = Base_ResidEval,
                           int id_x = 0,
                           doublereal delta_x = 0.0)
  {
    resid[0]  = RCONST(-0.04) * y[0] + RCONST(1.0e4) * y[1] * y[2];
    resid[1]  = -resid[0] - RCONST(3.0e7) * y[1] * y[1] - ydot[1];
    resid[0] -=  ydot[0];
    resid[2]  =  y[0] + y[1] + y[2] - ONE;
    return(0);
  }

  //! Get initial conditions
  virtual int getInitialConditions(const doublereal t0, doublereal * const y,
				   doublereal * const ydot) {
    y[0] = 1.0;
    y[1] = 0.0;
    y[2] = 0.0;
    ydot[0]  = -0.04;
    ydot[1]  =  0.04;
    ydot[2]  = 0.0;  
    return 0;
  }


  virtual int evalJacobianDP(const doublereal t, const doublereal delta_t, doublereal cj,
			     const doublereal* const y,
			     const doublereal* const ydot,
			     doublereal * const *jacDP,
			     doublereal * const resid) {
    
    const double *yval = y;
    
    double *jacCol = jacDP[0];
    jacCol[0] = RCONST(-0.04) - cj;
    jacCol[1] = RCONST(0.04);
    jacCol[2] = 1.0;

    jacCol = jacDP[1];
    jacCol[0] = RCONST(1.0e4)*yval[2];
    jacCol[1] = RCONST(-1.0e4)*yval[2] - RCONST(6.0e7)*yval[1] - cj;
    jacCol[2] = 1.0;

    jacCol = jacDP[2];
    jacCol[0] = RCONST(1.0e4)*yval[1];
    jacCol[1] = RCONST(-1.0e4)*yval[1];
    jacCol[2] = 1.0;

    return 0;
  }


};



static void PrintHeaderp(realtype rtol, double *atval, double *yval)
{
 

  printf("\nidadenx: Robertson kinetics DAE serial example problem for IDA\n");
  printf("         Three equation chemical kinetics problem.\n\n");
  printf("Linear solver: IDADENSE, with user-supplied Jacobian.\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Tolerance parameters:  rtol = %Lg   atol = %Lg %Lg %Lg \n",
         rtol, atval[0],atval[1],atval[2]);
  printf("Initial conditions y0 = (%Lg %Lg %Lg)\n",
         yval[0], yval[1], yval[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Tolerance parameters:  rtol = %lg   atol = %lg %lg %lg \n",
         rtol, atval[0],atval[1],atval[2]);
  printf("Initial conditions y0 = (%lg %lg %lg)\n",
         yval[0], yval[1], yval[2]);
#else
  printf("Tolerance parameters:  rtol = %g   atol = %g %g %g \n",
         rtol, atval[0],atval[1],atval[2]);
  printf("Initial conditions y0 = (%g %g %g)\n",
         yval[0], yval[1], yval[2]);
#endif
  printf("Constraints and id not used.\n\n");
  printf("-----------------------------------------------------------------------\n");
  printf("  t             y1           y2           y3");
  printf("      | nst  k      h\n");
  printf("-----------------------------------------------------------------------\n");
}

static void PrintOutputp(void *mem, realtype t, const double *yval)
{
  int retval, kused;
  long int nst;
  realtype hused;

  retval = IDAGetLastOrder(mem, &kused);
  check_flag(&retval, "IDAGetLastOrder", 1);
  retval = IDAGetNumSteps(mem, &nst);
  check_flag(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetLastStep(mem, &hused);
  check_flag(&retval, "IDAGetLastStep", 1);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%10.4Le %12.4Le %12.4Le %12.4Le | %3ld  %1d %12.4Le\n", 
         t, yval[0], yval[1], yval[2], nst, kused, hused);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%10.4le %12.4le %12.4le %12.4le | %3ld  %1d %12.4le\n", 
         t, yval[0], yval[1], yval[2], nst, kused, hused);
#else
  printf("%10.4e %12.4e %12.4e %12.4e | %3ld  %1d %12.4e\n", 
         t, yval[0], yval[1], yval[2], nst, kused, hused);
#endif
}



/*
 *--------------------------------------------------------------------
 * Main Program
 *--------------------------------------------------------------------
 */

int main(void)
{
  void *mem;
  N_Vector yy, yp, avtol;
  realtype rtol, *yval, *ypval, *atval;
  realtype t0,  tret;
  int iout, retval, retvalr;
  int rootsfound[2];
  double tout;


  RobertsFunc f;
  
  try {


    double atval[10];
    double yval[10];
    double ypval[10];
    yval[0] = 1.0;
    yval[1] = 0.0;
    yval[2] = 0.0;
    ypval[0]  = RCONST(-0.04);
    ypval[1]  = RCONST(0.04);
    ypval[2]  = ZERO;  
    double rtol = RCONST(1.0e-4);
    atval[0] = RCONST(1.0e-8);
    atval[1] = RCONST(1.0e-14);
    atval[2] = RCONST(1.0e-6);
    
    double t0 = 0.0;
    double tout1 = 0.4;
    tout = tout1;

    Cantera::DAE_Solver *solver = newDAE_Solver("IDA", f);

    Cantera::IDA_Solver *idasolver =  dynamic_cast<Cantera::IDA_Solver *>(solver);
    solver->setMaxOrder(5);
    solver->setTolerances(rtol, atval);
    idasolver->setJacobianType(1);
    solver->init(t0);

    solver->setTolerances(rtol, atval);
  

    void *mem = idasolver->IDAMemory();
    iout = 0;

    PrintHeaderp(rtol, atval, yval);

    while(1) {
      
      // retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);


      retval = solver->solve(tout);
      if (retval != IDA_SUCCESS && retval != IDA_TSTOP_RETURN) {
	printf("Error return = %d\n", retval);
      } else {
	tret = tout;
      }

      const double  *yretn = solver->solutionVector();

      PrintOutputp(mem, tret, yretn);




      if (retval == IDA_SUCCESS || retval == IDA_TSTOP_RETURN) {
	iout++;
	tout *= RCONST(10.0);
      }

      if (iout == 12) break;
    }

    PrintFinalStats(mem);

    delete solver;
    appdelete();


  }
  catch (CanteraError) {
    showErrors(cout);
  }


  

  mem = NULL;
  yy = yp = avtol = NULL;
  yval = ypval = atval = NULL;

  /* Allocate N-vectors. */
  yy = N_VNew_Serial(NEQ);
  if(check_flag((void *)yy, "N_VNew_Serial", 0)) return(1);
  yp = N_VNew_Serial(NEQ);
  if(check_flag((void *)yp, "N_VNew_Serial", 0)) return(1);
  avtol = N_VNew_Serial(NEQ);
  if(check_flag((void *)avtol, "N_VNew_Serial", 0)) return(1);

  /* Create and initialize  y, y', and absolute tolerance vectors. */
  yval  = NV_DATA_S(yy);
  yval[0] = ONE;
  yval[1] = ZERO;
  yval[2] = ZERO;

  ypval = NV_DATA_S(yp);
  ypval[0]  = RCONST(-0.04);
  ypval[1]  = RCONST(0.04);
  ypval[2]  = ZERO;  

  rtol = RCONST(1.0e-4);

  atval = NV_DATA_S(avtol);
  atval[0] = RCONST(1.0e-8);
  atval[1] = RCONST(1.0e-14);
  atval[2] = RCONST(1.0e-6);

  /* Integration limits */
  t0 = ZERO;
  double tout1 = RCONST(0.4);

  PrintHeader(rtol, avtol, yy);

  /* Call IDACreate and IDAMalloc to initialize IDA memory */
  mem = IDACreate();
  if(check_flag((void *)mem, "IDACreate", 0)) return(1);
  retval = IDAInit(mem, resrob, t0, yy, yp);
  if(check_flag(&retval, "IDAInit", 1)) return(1);
  retval = IDASVtolerances(mem, rtol, avtol);
  if(check_flag(&retval, "IDASVtolerances", 1)) return(1);

  /* Free avtol */
  N_VDestroy_Serial(avtol);

  /* Call IDARootInit to specify the root function grob with 2 components */
  retval = IDARootInit(mem, 2, grob);
  if (check_flag(&retval, "IDARootInit", 1)) return(1);

  /* Call IDADense and set up the linear solver. */
  retval = IDADense(mem, NEQ);
  if(check_flag(&retval, "IDADense", 1)) return(1);
  retval = IDADlsSetDenseJacFn(mem, jacrob);
  if(check_flag(&retval, "IDADlsSetDenseJacFn", 1)) return(1);

  /* In loop, call IDASolve, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached. */

  iout = 0; tout = tout1;
  while(1) {

    retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);

    PrintOutput(mem,tret,yy);

    if(check_flag(&retval, "IDASolve", 1)) return(1);

    if (retval == IDA_ROOT_RETURN) {
      retvalr = IDAGetRootInfo(mem, rootsfound);
      check_flag(&retvalr, "IDAGetRootInfo", 1);
      PrintRootInfo(rootsfound[0],rootsfound[1]);
    }

    if (retval == IDA_SUCCESS) {
      iout++;
      tout *= RCONST(10.0);
    }

    if (iout == NOUT) break;
  }

  PrintFinalStats(mem);

  /* Free memory */

  IDAFree(&mem);
  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);

  return(0);
  
}

/*
 *--------------------------------------------------------------------
 * Functions called by IDA
 *--------------------------------------------------------------------
 */

/*
 * Define the system residual function. 
 */

int resrob(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)
{
  realtype *yval, *ypval, *rval;

  yval = NV_DATA_S(yy); 
  ypval = NV_DATA_S(yp); 
  rval = NV_DATA_S(rr);

  rval[0]  = RCONST(-0.04)*yval[0] + RCONST(1.0e4)*yval[1]*yval[2];
  rval[1]  = -rval[0] - RCONST(3.0e7)*yval[1]*yval[1] - ypval[1];
  rval[0] -=  ypval[0];
  rval[2]  =  yval[0] + yval[1] + yval[2] - ONE;

  return(0);
}

/*
 * Root function routine. Compute functions g_i(t,y) for i = 0,1. 
 */

static int grob(realtype t, N_Vector yy, N_Vector yp, realtype *gout,
                void *user_data)
{
  realtype *yval, y1, y3;

  yval = NV_DATA_S(yy); 
  y1 = yval[0]; y3 = yval[2];
  gout[0] = y1 - RCONST(0.0001);
  gout[1] = y3 - RCONST(0.01);

  return(0);
}

/*
 * Define the Jacobian function. 
 */

int jacrob(int Neq, realtype tt,  realtype cj, 
           N_Vector yy, N_Vector yp, N_Vector resvec,
           DlsMat JJ, void *user_data,
           N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{
  realtype *yval;
  
  yval = NV_DATA_S(yy);

  IJth(JJ,1,1) = RCONST(-0.04) - cj;
  IJth(JJ,2,1) = RCONST(0.04);
  IJth(JJ,3,1) = ONE;
  IJth(JJ,1,2) = RCONST(1.0e4)*yval[2];
  IJth(JJ,2,2) = RCONST(-1.0e4)*yval[2] - RCONST(6.0e7)*yval[1] - cj;
  IJth(JJ,3,2) = ONE;
  IJth(JJ,1,3) = RCONST(1.0e4)*yval[1];
  IJth(JJ,2,3) = RCONST(-1.0e4)*yval[1];
  IJth(JJ,3,3) = ONE;

  return(0);
}

/*
 *--------------------------------------------------------------------
 * Private functions
 *--------------------------------------------------------------------
 */

/* 
 * Print first lines of output (problem description)
 */

static void PrintHeader(realtype rtol, N_Vector avtol, N_Vector y)
{
  realtype *atval, *yval;

  atval  = NV_DATA_S(avtol);
  yval  = NV_DATA_S(y);

  printf("\nidadenx: Robertson kinetics DAE serial example problem for IDA\n");
  printf("         Three equation chemical kinetics problem.\n\n");
  printf("Linear solver: IDADENSE, with user-supplied Jacobian.\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Tolerance parameters:  rtol = %Lg   atol = %Lg %Lg %Lg \n",
         rtol, atval[0],atval[1],atval[2]);
  printf("Initial conditions y0 = (%Lg %Lg %Lg)\n",
         yval[0], yval[1], yval[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Tolerance parameters:  rtol = %lg   atol = %lg %lg %lg \n",
         rtol, atval[0],atval[1],atval[2]);
  printf("Initial conditions y0 = (%lg %lg %lg)\n",
         yval[0], yval[1], yval[2]);
#else
  printf("Tolerance parameters:  rtol = %g   atol = %g %g %g \n",
         rtol, atval[0],atval[1],atval[2]);
  printf("Initial conditions y0 = (%g %g %g)\n",
         yval[0], yval[1], yval[2]);
#endif
  printf("Constraints and id not used.\n\n");
  printf("-----------------------------------------------------------------------\n");
  printf("  t             y1           y2           y3");
  printf("      | nst  k      h\n");
  printf("-----------------------------------------------------------------------\n");
}

/*
 * Print Output
 */

static void PrintOutput(void *mem, realtype t, N_Vector y)
{
  realtype *yval;
  int retval, kused;
  long int nst;
  realtype hused;

  yval  = NV_DATA_S(y);

  retval = IDAGetLastOrder(mem, &kused);
  check_flag(&retval, "IDAGetLastOrder", 1);
  retval = IDAGetNumSteps(mem, &nst);
  check_flag(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetLastStep(mem, &hused);
  check_flag(&retval, "IDAGetLastStep", 1);
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("%10.4Le %12.4Le %12.4Le %12.4Le | %3ld  %1d %12.4Le\n", 
         t, yval[0], yval[1], yval[2], nst, kused, hused);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("%10.4le %12.4le %12.4le %12.4le | %3ld  %1d %12.4le\n", 
         t, yval[0], yval[1], yval[2], nst, kused, hused);
#else
  printf("%10.4e %12.4e %12.4e %12.4e | %3ld  %1d %12.4e\n", 
         t, yval[0], yval[1], yval[2], nst, kused, hused);
#endif
}

static void PrintRootInfo(int root_f1, int root_f2)
{
  printf("    rootsfound[] = %3d %3d\n", root_f1, root_f2);
  return;
}

/*
 * Print final integrator statistics
 */

static void PrintFinalStats(void *mem)
{
  int retval;
  long int nst, nni, nje, nre, nreLS, netf, ncfn, nge;

  retval = IDAGetNumSteps(mem, &nst);
  check_flag(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetNumResEvals(mem, &nre);
  check_flag(&retval, "IDAGetNumResEvals", 1);
  retval = IDADlsGetNumJacEvals(mem, &nje);
  check_flag(&retval, "IDADlsGetNumJacEvals", 1);
  retval = IDAGetNumNonlinSolvIters(mem, &nni);
  check_flag(&retval, "IDAGetNumNonlinSolvIters", 1);
  retval = IDAGetNumErrTestFails(mem, &netf);
  check_flag(&retval, "IDAGetNumErrTestFails", 1);
  retval = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
  check_flag(&retval, "IDAGetNumNonlinSolvConvFails", 1);
  retval = IDADlsGetNumResEvals(mem, &nreLS);
  check_flag(&retval, "IDADlsGetNumResEvals", 1);
  retval = IDAGetNumGEvals(mem, &nge);
  check_flag(&retval, "IDAGetNumGEvals", 1);

  printf("\nFinal Run Statistics: \n\n");
  printf("Number of steps                    = %ld\n", nst);
  printf("Number of residual evaluations     = %ld\n", nre+nreLS);
  printf("Number of Jacobian evaluations     = %ld\n", nje);
  printf("Number of nonlinear iterations     = %ld\n", nni);
  printf("Number of error test failures      = %ld\n", netf);
  printf("Number of nonlinear conv. failures = %ld\n", ncfn);
  printf("Number of root fn. evaluations     = %ld\n", nge);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;
  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, 
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1);
  } else if (opt == 1) {
    /* Check if flag < 0 */
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, 
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", 
              funcname, *errflag);
      return(1); 
    }
  } else if (opt == 2 && flagvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr, 
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1);
  }

  return(0);
}
