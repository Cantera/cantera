//! @file EigenSparseDirectJacobian.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef EIGENSPARSEDIRECTJACOBIAN_H
#define EIGENSPARSEDIRECTJACOBIAN_H

#include "cantera/numerics/EigenSparseJacobian.h"

namespace Cantera
{

//! A system matrix solver that uses Eigen's sparse direct (LU) algorithm
class EigenSparseDirectJacobian : public EigenSparseJacobian
{
public:
    EigenSparseDirectJacobian() = default;
    const string type() const override { return "eigen-sparse-direct"; }
    void factorize() override;
    void solve(span<const double> rhs_vector, span<double> output) override;

protected:
    Eigen::SparseLU<Eigen::SparseMatrix<double>> m_solver;
};

}

#endif
