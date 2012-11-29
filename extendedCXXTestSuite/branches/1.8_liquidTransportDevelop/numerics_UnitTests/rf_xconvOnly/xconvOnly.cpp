//  Example 
//
// This is an example of a discontinuous function. The root finder should still
// find the solution given hints at the monotonicity of the solution.
//

#include <string>
#include <map>

#include <cantera/Cantera.h>
#include <cantera/numerics.h>
#include <cantera/kernel/ResidEval.h>    
#include <cantera/kernel/RootFind.h>

using namespace Cantera;
using namespace std;


void printDbl(double val) {
  if (fabs(val) < 5.0E-17) {
    cout << " nil";
  } else {
    cout << val;
  }
}



class funcDiscontinuous : public ResidEval {

public:
  funcDiscontinuous() :
    ResidEval(),
    neq_(1)
  {
  
  }

  int nEquations() const 
  {
    return neq_;
  }

  int  evalSS(const doublereal t, const doublereal * const y,
	      doublereal * const r)
  {
    double x = y[0];
    
    if (x > 0.5) {
      r[0] = 1.0;
    } else {
      r[0] = -1.0;
    }
    return 1; 
  }

  int nEquations() {
    return 1;
  }

  int neq_;

};



int main(int argc, char** argv) {
  string infile = "diamond.xml";

  try {


    // Define a residual. The definition of a residual involves a lot more work than is shown here.
    funcDiscontinuous * func_ptr = new funcDiscontinuous();
    // Instantiate the root finder with the residual to be solved, ec.
    RootFind rf(func_ptr);
    // Set the relative and absolute tolerancess for f and x.
    rf.setTol(1.0E-5, 1.0E-10, 1.0E-5, 1.0E-11);
    // Give a hint about the function's dependence on x. This is needed, for example, if the function has
    // flat regions.
    rf.setFuncIsGenerallyIncreasing(true);
    rf.setDeltaX(0.1);
    // Supply an initial guess for the solution
    double xbest = 0.12;

    // Set the print level for the solver. Zero produces no output. Two produces a summary table of each iteration.
    rf.setPrintLvl(2);
    // Define a minimum and maximum for the independent variable.
    double xmin = -2.0;
    double xmax = 5.0;
    // Define a maximum iteration number
    int itmax = 100;
    // Define the f_0 value, and on return will contain the actual value of f(x) obtained
    double f_value = 0.1;
    // Call the solver
    int status = rf.solve(xmin, xmax, itmax, f_value, &xbest);
    if (status >= 0) {
      printf(" f (%g) = %g\n", xbest, f_value);
      if (status == 1) {
	printf("  function returned correct status of Xonly converged\n");
      }
    } else {  
      printf(" bad status = %d f (%g) = %g\n",  status, xbest, f_value);
    }
    appdelete();

  }
  catch (CanteraError) {
    showErrors(cout);
  }

  return 0;
}
/***********************************************************/
