/**
 *  @file GeneralMatrix.cpp
 *
 */
/*
 * $Revision: 725 $
 * $Date: 2011-05-16 18:45:08 -0600 (Mon, 16 May 2011) $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "GeneralMatrix.h"
using namespace std;

namespace Cantera {
  //====================================================================================================================
  GeneralMatrix::GeneralMatrix(int matType) :
    matrixType_(matType)
  {
  }
  //====================================================================================================================
  GeneralMatrix::GeneralMatrix(const GeneralMatrix &y) :
    matrixType_(y.matrixType_)
  {
  }
  //====================================================================================================================
  GeneralMatrix&  GeneralMatrix::operator=(const GeneralMatrix &y)
  { 
    if (&y == this) return *this;
    matrixType_ = y.matrixType_;
    return *this;
  }
  //====================================================================================================================
  GeneralMatrix::~GeneralMatrix()
  {
  }
  //====================================================================================================================
}
