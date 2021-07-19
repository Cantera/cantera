//! @file AdaptivePreconditioner.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/global.h"
#include "cantera/numerics/AdaptivePreconditioner.h"

namespace Cantera
{

void AdaptivePreconditioner::setValue(size_t row, size_t col, double value)
{
    m_jac_trips.emplace_back(row, col, value);
}

void AdaptivePreconditioner::stateAdjustment(vector_fp& state) {
    // Only keep positive composition based on given tol
    for (size_t i = 0; i < state.size(); i++) {
        state[i] = std::max(state[i], m_atol);
    }
}

void AdaptivePreconditioner::initialize(size_t networkSize)
{
    // don't use legacy rate constants
    use_legacy_rate_constants(false);
    // reset arrays in case of re-initialization
    m_jac_trips.clear();
    // set dimensions of preconditioner from network
    m_dim = networkSize;
    // reserve some space for vectors making up SparseMatrix
    m_jac_trips.reserve(3 * networkSize);
    // reserve space for preconditioner
    m_precon_matrix.resize(m_dim, m_dim);
    // creating sparse identity matrix
    m_identity.resize(m_dim, m_dim);
    m_identity.setIdentity();
    m_identity.makeCompressed();
    // setting default ILUT parameters
    if (m_drop_tol == 0) {
        setIlutDropTol(1e-10);
    }
    if (m_drop_tol == 0) {
        setIlutFillFactor(m_dim/4);
    }
    // update initialized status
    m_init = true;
}


void AdaptivePreconditioner::setup()
{
    // make into preconditioner as P = (I - gamma * J_bar)
    transformJacobianToPreconditioner();
    // compressing sparse matrix structure
    m_precon_matrix.makeCompressed();
    // analyze and factorize
    m_solver.compute(m_precon_matrix);
    // check for errors
    if (m_solver.info() != Eigen::Success) {
        throw CanteraError("AdaptivePreconditioner::setup",
                           "error code: {}", m_solver.info());
    }
}

void AdaptivePreconditioner::transformJacobianToPreconditioner()
{
    // set precon to jacobian
    m_precon_matrix.setFromTriplets(m_jac_trips.begin(), m_jac_trips.end());
    // convert to preconditioner
    m_precon_matrix = m_identity - m_gamma * m_precon_matrix;
    // prune by threshold if desired
    if (m_prune_precon) {
        prunePreconditioner();
    }
}

void AdaptivePreconditioner::prunePreconditioner()
{
    for (int k=0; k<m_precon_matrix.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(m_precon_matrix, k); it; ++it) {
            if (std::abs(it.value()) < m_threshold && it.row() != it.col()) {
                m_precon_matrix.coeffRef(it.row(), it.col()) = 0;
            }
        }
    }
}

void AdaptivePreconditioner::solve(const size_t stateSize, double *rhs_vector, double* output)
{
    // creating vectors in the form of Ax=b
    Eigen::Map<Eigen::VectorXd> bVector(rhs_vector, stateSize);
    Eigen::Map<Eigen::VectorXd> xVector(output, stateSize);
    // solve for xVector
    xVector = m_solver.solve(bVector);
    if (m_solver.info() != Eigen::Success) {
        throw CanteraError("AdaptivePreconditioner::solve",
                           "error code: {}", m_solver.info());
    }
}

}
