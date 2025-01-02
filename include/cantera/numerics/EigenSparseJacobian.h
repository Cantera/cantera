//! @file EigenSparseJacobian.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef EIGENSPARSEJACOBIAN_H
#define EIGENSPARSEJACOBIAN_H

#include "cantera/numerics/SystemJacobian.h"
#include "cantera/numerics/eigen_sparse.h"

namespace Cantera
{

//! System Jacobians that use Eigen sparse matrices for storage
class EigenSparseJacobian : public SystemJacobian
{
public:
    EigenSparseJacobian() = default;
    void initialize(size_t networkSize) override;
    void setValue(size_t row, size_t col, double value) override;
    void updatePreconditioner() override;
    void updateTransient(double rdt, int* mask) override;

    //! Return underlying Jacobian matrix
    //! @ingroup derivGroup
    Eigen::SparseMatrix<double> jacobian();

    //! Return the internal preconditioner matrix
    Eigen::SparseMatrix<double> matrix() {
        updatePreconditioner();
        return m_matrix;
    }

    //! Print preconditioner contents
    void printPreconditioner() override;

    //! Print jacobian contents
    void printJacobian();

protected:
    //! Vector of triples representing the jacobian used in preconditioning
    vector<Eigen::Triplet<double>> m_jac_trips;

    //! Storage of appropriately sized identity matrix for making the preconditioner
    Eigen::SparseMatrix<double> m_identity;

    //! Container that is the sparse preconditioner
    Eigen::SparseMatrix<double> m_matrix;
};

}

#endif
