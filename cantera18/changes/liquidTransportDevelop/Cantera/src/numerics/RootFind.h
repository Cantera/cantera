/**
 * @file RootFind.h
 *       Header file for implicit nonlinear solver of a one dimensional function
 *  (see \ref numerics and class \link Cantera::RootFind RootFind\endlink).
 */
/*
 * $Id$
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef CT_ROOTFIND_H
#define CT_ROOTFIND_H
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
  
    int solve(doublereal xmin, doublereal xmax, int itmax, doublereal funcTargetValue, doublereal *xbest) ;

    doublereal func(doublereal x);

    void setTol(doublereal rtol, doublereal atol);

    void setPrintLvl(int printLvl);
    void setFuncIsGenerallyIncreasing(bool value);
    void setFuncIsGenerallyDecreasing(bool value);
    void setDeltaX(doublereal deltaXNorm); 

  public:
    ResidEval *m_residFunc;
    doublereal m_funcTargetValue;
    doublereal m_atol;
    doublereal m_rtol;
    doublereal m_maxstep;
  protected:
    int    printLvl;
    doublereal DeltaXnorm_;
    bool FuncIsGenerallyIncreasing_;
    bool FuncIsGenerallyDecreasing_;

  };
}
#endif
