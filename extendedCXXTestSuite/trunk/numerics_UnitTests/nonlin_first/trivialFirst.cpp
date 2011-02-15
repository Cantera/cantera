/**
 *  @file mixGasTransport.cpp
 *       test problem for mixture transport
 */

//  Example 
//
// Test case for mixture transport in a gas
// The basic idea is to set up a gradient of some kind.
// Then the resulting transport coefficients out.
// Essentially all of the interface routines should be
// exercised and the results dumped out.
//
// A blessed solution test will make sure that the actual
// solution doesn't change as a function of time or 
// further development. 

// perhaps, later, an analytical solution could be added

// An open Rankine cycle

#include <string>
#include <map>

#include <cantera/Cantera.h>
#include <cantera/numerics.h>
#include <cantera/kernel/ResidJacEval.h>    // defines class Water
#include <cantera/kernel/NonlinearSolver.h>

using namespace Cantera;
using namespace std;


void printDbl(double val) {
  if (fabs(val) < 5.0E-17) {
    cout << " nil";
  } else {
    cout << val;
  }
}



class RoboFunc : public ResidJacEval {

public:
  RoboFunc(doublereal atol = 1.0e-13) :
    ResidJacEval(atol)
  {
    neq_ = 8;
  }

  int nEquations() const 
  {
    return neq_;
  }

  virtual int  evalResidNJ(doublereal t, const doublereal deltaT,
			   const doublereal * const y,
			   const doublereal * const ydot,
			   doublereal * const resid,
			   ResidEval_Type_Enum evalType = Base_ResidEval,
			   int id_x = 0,
			   doublereal delta_x = 0.0)
  {

    double x1, x2, x3, x4, x5, x6, x7, x8;
    double eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8;

    x1 = y[0]; 
    x2 = y[1]; 
    x3 = y[2]; 
    x4 = y[3]; 
    x5 = y[4]; 
    x6 = y[5]; 
    x7 = y[6]; 
    x8 = y[7]; 

    eq1 = - 0.1238*x1 + x7 - 0.001637*x2 
      - 0.9338*x4 + 0.004731*x1*x3 - 0.3578*x2*x3 - 0.3571;
    eq2 = 0.2638*x1 - x7 - 0.07745*x2 
      - 0.6734*x4 + 0.2238*x1*x3 + 0.7623*x2*x3 - 0.6022;
    eq3 = 0.3578*x1 + 0.004731*x2 + x6*x8;
    eq4 = - 0.7623*x1 + 0.2238*x2 + 0.3461;
    eq5 = x1*x1 + x2*x2 - 1;
    eq6 = x3*x3 + x4*x4 - 1;
    eq7 = x5*x5 + x6*x6 - 1;
    eq8 = x7*x7 + x8*x8 - 1;
  
    resid[0] = eq1;
    resid[1] = eq2; 
    resid[2] = eq3; 
    resid[3] = eq4; 
    resid[4] = eq5; 
    resid[5] = eq6;
    resid[6] = eq7; 
    resid[7] = eq8; 

    return 1; 
  }





};



int main(int argc, char** argv) {
  int k;
  string infile = "diamond.xml";

  try {
    int num_newt_its = 0;
    int num_linear_solves = 0;
    int numBacktracks = 0;
    int loglevelInput = 9;
    RoboFunc r1;
    int neq = r1.nEquations();
    double y_comm[8];
    double deltaBounds[8];
    double ydot_comm[8];;
    double time_curr = 0.0;
    SquareMatrix jac(8);
    double CJ = 0.0;
    for (k = 0; k < neq; k++) {
      y_comm[k] = 0.1;
      deltaBounds[k] = 0.1;
    }
  
    NonlinearSolver *nls = new NonlinearSolver(&r1);

    int solnType =       NSOLN_TYPE_STEADY_STATE ;

    nls->setDeltaBoundsMagnitudes(deltaBounds);

    nls->solve_nonlinear_problem(solnType, y_comm,
				 ydot_comm, CJ,
				 time_curr,
				 jac,
				 num_newt_its,
				 num_linear_solves,
				 numBacktracks,
				 loglevelInput);
 

    double res[10];

    r1.evalResidNJ(0.0, 0.0, y_comm, ydot_comm, res);
 
    printf("\nSOLN:\n");

    printf(" Unk    Value        residual if GE 1.0E-15\n");
 
    for (k = 0; k < neq; k++) {
   
      double tmp = res[k];
      if (fabs(tmp) < 1.0E-15) {
	tmp = 0.0;
      }
      printf("    %d %13.5E %13.5E\n", k, y_comm[k], tmp);
 
    }

    delete nls;
    appdelete();
 

  }
  catch (CanteraError) {
    showErrors(cout);
  }

  return 0;
}
/***********************************************************/
