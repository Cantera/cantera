//! @file EigenSparseJacobian.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/EigenSparseJacobian.h"
#include "cantera/base/global.h"

namespace Cantera
{

void EigenSparseJacobian::initialize(size_t networkSize)
{
    // reset arrays in case of re-initialization
    m_jac_trips.clear();
    // set dimensions of preconditioner from network
    m_dim = networkSize;
    // reserve some space for vectors making up SparseMatrix
    m_jac_trips.reserve(3 * networkSize);
    // reserve space for preconditioner
    m_matrix.resize(m_dim, m_dim);
    // creating sparse identity matrix
    m_identity.resize(m_dim, m_dim);
    m_identity.setIdentity();
    m_identity.makeCompressed();
    // update initialized status
    m_init = true;
}

void EigenSparseJacobian::reset()
{
    m_matrix.setZero();
    m_jac_trips.clear();
}

void EigenSparseJacobian::setValue(size_t row, size_t col, double value)
{
    m_jac_trips.emplace_back(static_cast<int>(row), static_cast<int>(col), value);
}

void EigenSparseJacobian::updatePreconditioner()
{
    // set precon to jacobian
    m_matrix.setFromTriplets(m_jac_trips.begin(), m_jac_trips.end());
    // make into preconditioner as P = (I - gamma * J_bar)
    m_matrix = m_identity - m_gamma * m_matrix;
    // prune by threshold if desired
    factorize();
}

void EigenSparseJacobian::updateTransient(double rdt, span<const int> mask)
{
    // set matrix to steady Jacobian
    m_matrix.setFromTriplets(m_jac_trips.begin(), m_jac_trips.end());
    // update transient diagonal terms
    Eigen::VectorXd diag = Eigen::Map<const Eigen::VectorXi>(mask.data(), m_dim).cast<double>();
    m_matrix -= rdt * diag.matrix().asDiagonal();
    factorize();
}

Eigen::SparseMatrix<double> EigenSparseJacobian::jacobian()
{
    Eigen::SparseMatrix<double> jacobian_mat(m_dim, m_dim);
    jacobian_mat.setFromTriplets(m_jac_trips.begin(), m_jac_trips.end());
    return jacobian_mat;
}

void EigenSparseJacobian::printPreconditioner() {
    std::stringstream ss;
    Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
    ss << Eigen::MatrixXd(m_matrix).format(HeavyFmt);
    writelog(ss.str());
}

void EigenSparseJacobian::printJacobian() {
    std::stringstream ss;
    Eigen::SparseMatrix<double> jacobian(m_dim, m_dim);
    jacobian.setFromTriplets(m_jac_trips.begin(), m_jac_trips.end());
    ss << Eigen::MatrixXd(jacobian);
    writelog(ss.str());
}

}
