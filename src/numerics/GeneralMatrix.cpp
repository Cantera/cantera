/**
 *  @file GeneralMatrix.cpp
 */

/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "cantera/numerics/GeneralMatrix.h"
using namespace std;

namespace Cantera
{

GeneralMatrix::GeneralMatrix(int matType) :
    matrixType_(matType),
    m_factored(0)
{
}

GeneralMatrix::GeneralMatrix(const GeneralMatrix& y) :
    matrixType_(y.matrixType_),
    m_factored(y.m_factored)
{
}

GeneralMatrix&  GeneralMatrix::operator=(const GeneralMatrix& y)
{
    if (&y == this) {
        return *this;
    }
    matrixType_ = y.matrixType_;
    m_factored = y.m_factored;
    return *this;
}

}
