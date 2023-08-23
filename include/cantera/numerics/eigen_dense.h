// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_EIGEN_DENSE_H
#define CT_EIGEN_DENSE_H

#include "cantera/base/config.h"
#if CT_USE_SYSTEM_EIGEN
    #if CT_USE_SYSTEM_EIGEN_PREFIXED
    #include <eigen3/Eigen/Dense>
    #else
    #include <Eigen/Dense>
    #endif
#else
#include "cantera/ext/Eigen/Dense"
#endif

namespace Cantera
{

//! @addtogroup matrices
//! @{

typedef Eigen::Map<Eigen::MatrixXd> MappedMatrix;
typedef Eigen::Map<const Eigen::MatrixXd> ConstMappedMatrix;
typedef Eigen::Map<Eigen::VectorXd> MappedVector;
typedef Eigen::Map<const Eigen::VectorXd> ConstMappedVector;
typedef Eigen::Map<Eigen::RowVectorXd> MappedRowVector;
typedef Eigen::Map<const Eigen::RowVectorXd> ConstMappedRowVector;

//! @}

}

#endif
