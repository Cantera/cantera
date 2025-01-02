//! @file EigenSparseDirectJacobian.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/EigenSparseDirectJacobian.h"
#include "cantera/numerics/eigen_dense.h"

namespace Cantera
{

void EigenSparseDirectJacobian::factorize()
{
    m_matrix.makeCompressed();
    // analyze and factorize
    m_solver.compute(m_matrix);
    // check for errors
    if (m_solver.info() != Eigen::Success) {
        throw CanteraError("EigenSparseDirectJacobian::factorize",
                           "error code: {}", static_cast<int>(m_solver.info()));
    }
}

void EigenSparseDirectJacobian::solve(const size_t stateSize, double* b, double* x)
{
    MappedVector(x, m_dim) = m_solver.solve(MappedVector(b, m_dim));
    // check for errors
    if (m_solver.info() != Eigen::Success) {
        throw CanteraError("EigenSparseDirectJacobian::solve",
                           "error code: {}", static_cast<int>(m_solver.info()));
    }
}

}
