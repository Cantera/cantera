/**
 *  @file IDA_Solver.h
 *
 */

// Copyright 2006  California Institute of Technology

using namespace std;

namespace Cantera {
  
  /**
   * A simple class to hold an array of parameter values and a pointer to 
   * an instance of a subclass of ResidEval.
   */
  class ResidData {
    
  public:
    ResidData(ResidEval* f, int npar = 0) ;
    virtual ~ResidData();
    ResidEval* func_;
  };

}


extern "C" {
  
  /**
   *  Function called by IDA to evaluate the residual, given y and
   *  ydot.  IDA allows passing in a void* pointer to access
   *  external data. Instead of requiring the user to provide a
   *  residual function directly to IDA (which would require using
   *  the sundials data types N_Vector, etc.), we define this
   *  function as the single function that IDA always calls. The
   *  real evaluation of the residual is done by an instance of a
   *  subclass of ResidEval, passed in to this function as a pointer
   *  in the parameters.
   */
  static int ida_resid(realtype t, N_Vector y, N_Vector ydot, 
		       N_Vector r, void *f_data) {
    double* ydata = NV_DATA_S(y);
    double* ydotdata = NV_DATA_S(ydot);
    double* rdata = NV_DATA_S(r);
    Cantera::ResidData* d = (Cantera::ResidData*)f_data;
    Cantera::ResidEval* f = d->func;
    f->eval(t, ydata, ydotdata, rdata);
    return 0;
  }
  
}

namespace Cantera {
  
  /**
   *  Constructor. Default settings: dense jacobian, no user-supplied
   *  Jacobian function, Newton iteration.
   */
  IDA_Solver::IDA_Solver(ResidEval& f);
  
  /// Destructor.
  IDA_Solver::~IDA_Solver();
  
  Real IDA_Solver::solution(int k) const;
  
  const Real* IDA_Solver::solutionVector() const;
  
  Real IDA_Solver::derivative(int k) const;
  
  const Real* IDA_Solver::derivativeVector() const;
  
  
  void IDA_Solver::setTolerances(double reltol, double* abstol) ;
  void IDA_Solver::setTolerances(double reltol, double abstol);

  void IDA_Solver::setLinearSolverType(int solverType);

  void IDA_Solver::init(double );

  void IDA_Solver::solve(double tout);

  double IDA_Solver::step(double tout);

  Real IDA_Solver::getOutputParameter(int flag);
}

