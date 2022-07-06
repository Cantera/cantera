// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PY_WRAPPERS_H
#define CT_PY_WRAPPERS_H

#include "cantera/base/Array.h"

// Function for assigning elements of Array2D, since Cython has trouble
// with assigning to the reference returned by operator()
inline void CxxArray2D_set(Cantera::Array2D& array, size_t i, size_t j, double value)
{
    array(i,j) = value;
}

// Function which populates a 1D array
#define ARRAY_FUNC(PREFIX, CLASS_NAME, FUNC_NAME) \
    inline void PREFIX ## _ ## FUNC_NAME(Cantera::CLASS_NAME* object, double* data) \
    { object->FUNC_NAME(data); }

// function which takes a stride as the first argument and populates a 2D array
#define ARRAY_FUNC2(PREFIX, CLASS_NAME, FUNC_NAME) \
    inline void PREFIX ## _ ## FUNC_NAME(Cantera::CLASS_NAME* object, size_t dim, double* data) \
    { object->FUNC_NAME(dim, data); }

// Function which populates a 1D array, extra arguments
#define ARRAY_POLY(PREFIX, CLASS_NAME, FUNC_NAME) \
    inline void PREFIX ## _ ## FUNC_NAME(Cantera::CLASS_NAME* object, size_t i, double* data) \
    { object->FUNC_NAME(i, data); }

#define ARRAY_POLY_BINARY(PREFIX, CLASS_NAME, FUNC_NAME) \
    inline void PREFIX ## _ ## FUNC_NAME(Cantera::CLASS_NAME* object, size_t i, size_t j, double* data) \
    { object->FUNC_NAME(i, j, data); }

#endif
