//! @file AdaptivePreconditioner.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/global.h"
#include "cantera/numerics/AdaptivePreconditioner.h"

namespace Cantera
{

AdaptivePreconditioner::AdaptivePreconditioner()
{
    setPreconditionerSide("right");
}

void AdaptivePreconditioner::stateAdjustment(vector<double>& state) {
    // Only keep positive composition based on given tol
    for (size_t i = 0; i < state.size(); i++) {
        state[i] = std::max(state[i], m_atol);
    }
}

void AdaptivePreconditioner::initialize(size_t networkSize)
{
    EigenSparseJacobian::initialize(networkSize);
    // don't use legacy rate constants
    use_legacy_rate_constants(false);
    // setting default ILUT parameters
    if (m_drop_tol == 0) {
        setIlutDropTol(1e-12);
    }
    if (m_drop_tol == 0) {
        setIlutFillFactor(static_cast<int>(m_dim) / 2);
    }
    // update initialized status
    m_init = true;
}

void AdaptivePreconditioner::factorize()
{
    if (m_prune_precon) {
        prunePreconditioner();
    }
    // compress sparse matrix structure
    m_matrix.makeCompressed();
    // analyze and factorize
    m_solver.compute(m_matrix);
    // check for errors
    if (m_solver.info() != Eigen::Success) {
        throw CanteraError("AdaptivePreconditioner::factorize",
                           "error code: {}", static_cast<int>(m_solver.info()));
    }
}

void AdaptivePreconditioner::prunePreconditioner()
{
    for (int k=0; k<m_matrix.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(m_matrix, k); it;
            ++it) {
            if (std::abs(it.value()) < m_threshold && it.row() != it.col()) {
                it.valueRef() = 0;
            }
        }
    }
}

void AdaptivePreconditioner::solve(const size_t stateSize, double* rhs_vector, double*
    output)
{
    // creating vectors in the form of Ax=b
    Eigen::Map<Eigen::VectorXd> bVector(rhs_vector, stateSize);
    Eigen::Map<Eigen::VectorXd> xVector(output, stateSize);
    // solve for xVector
    xVector = m_solver.solve(bVector);
    if (m_solver.info() != Eigen::Success) {
        throw CanteraError("AdaptivePreconditioner::solve",
                           "error code: {}", static_cast<int>(m_solver.info()));
    }
}

}
