// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_EIGEN_DEFS_H
#define CT_EIGEN_DEFS_H

#include "cantera/base/config.h"
#if CT_USE_SYSTEM_EIGEN
#include <Eigen/Dense>
#include <Eigen/Sparse>
#else
#include "cantera/ext/Eigen/Dense"
#include "cantera/ext/Eigen/Sparse"
#endif

namespace Cantera
{

typedef Eigen::Map<Eigen::MatrixXd> MappedMatrix;
typedef Eigen::Map<Eigen::VectorXd> MappedVector;
typedef Eigen::Map<const Eigen::VectorXd> ConstMappedVector;
typedef Eigen::Map<Eigen::RowVectorXd> MappedRowVector;
typedef Eigen::Map<const Eigen::RowVectorXd> ConstMappedRowVector;

typedef Eigen::Map<Eigen::SparseMatrix<double>> MappedSparseMatrix;
typedef Eigen::Map<const Eigen::SparseMatrix<double>> ConstMappedSparseMatrix;
}

#endif
