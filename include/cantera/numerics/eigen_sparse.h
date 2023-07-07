// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_EIGEN_SPARSE_H
#define CT_EIGEN_SPARSE_H

#include "cantera/base/config.h"
#if CT_USE_SYSTEM_EIGEN
    #if CT_USE_SYSTEM_EIGEN_PREFIXED
    #include <eigen3/Eigen/Sparse>
    #else
    #include <Eigen/Sparse>
    #endif
#else
#include "cantera/ext/Eigen/Sparse"
#endif

namespace Cantera
{
//! @ingroup matrices
typedef std::vector<Eigen::Triplet<double>> SparseTriplets;
}

#endif
