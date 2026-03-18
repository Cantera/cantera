// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_EIGEN_DENSE_H
#define CT_EIGEN_DENSE_H

#include "cantera/base/ct_defs.h"

#if CT_USE_SYSTEM_EIGEN
    #if CT_USE_SYSTEM_EIGEN_PREFIXED
    #include <eigen3/Eigen/Dense>
    #else
    #include <Eigen/Dense>
    #endif
#else
#include "cantera/ext/Eigen/Dense"
#endif
#include <concepts>

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

template<class Derived>
concept EigenDenseDouble = std::same_as<typename Derived::Scalar, double>;

// Require a 1D type with contiguous storage
template<class Derived>
concept EigenVectorLike = EigenDenseDouble<Derived>
    && (Derived::IsVectorAtCompileTime == 1);

//! Convenience wrapper for accessing Eigen vector/array/map data as a span
//! @todo Remove once %Cantera requires Eigen 5.0 or newer.
template<EigenVectorLike Derived>
inline span<double> asSpan(Eigen::DenseBase<Derived>& v)
{
    return span<double>(v.derived().data(), v.size());
}

//! Convenience wrapper for accessing Eigen vector/array/map data as a span
//! @todo Remove once %Cantera requires Eigen 5.0 or newer.
template<EigenVectorLike Derived>
inline span<const double> asSpan(const Eigen::DenseBase<Derived>& v)
{
    return span<const double>(v.derived().data(), v.size());
}

}

#endif
