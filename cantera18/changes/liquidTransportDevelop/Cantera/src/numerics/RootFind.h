/**
 * @file RootFind.h
 *       Header file for implicit nonlinear solver of a one dimensional function
 *  (see \ref numerics and class \link Cantera::RootFind RootFind\endlink).
 */
/*
 * $Id: solveSP.h 381 2010-01-15 21:20:41Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef ROOTFIND_H
#define ROOTFIND_H
/**
 * @defgroup solverGroup Solvers for Equation Systems
 */

#include <vector>
#include "ResidEval.h"

namespace Cantera {

#define ROOTFIND_SUCCESS             0
#define ROOTFIND_FAILEDCONVERGENCE  -1
#define ROOTFIND_BADINPUT           -2

  class RootFind {

  public:

    //! Constructor for the object
    /*!
     */ 
    RootFind(ResidEval* resid);
   
    //! Destructor. Deletes the integrator.
    ~RootFind();

  private:

    //! Unimplemented private copy constructor
    RootFind(const RootFind &right);

    //! Unimplemented private assignment operator
    RootFind& operator=(const RootFind &right);

  public:
  
    int solve(double xmin, double xmax, int itmax, double funcTargetValue, double *xbest) ;

    double func(double x) ;

    void setTol(double rtol, double atol);

    void setPrintLvl(int printLvl) ;

  public:
    ResidEval *m_residFunc;
    double m_funcTargetValue;
    double m_atol;
    double m_rtol;
    double m_maxstep;
    int    printLvl;

  };
}
#endif
