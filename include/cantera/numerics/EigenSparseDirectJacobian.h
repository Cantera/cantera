//! @file EigenSparseDirectJacobian.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef EIGENSPARSEDIRECTJACOBIAN_H
#define EIGENSPARSEDIRECTJACOBIAN_H

#include "cantera/numerics/EigenSparseJacobian.h"

namespace Cantera
{

class EigenSparseDirectJacobian : public EigenSparseJacobian
{
public:
    EigenSparseDirectJacobian() = default;
    const string type() const override { return "eigen-sparse-direct"; }
    void factorize() override;
    void solve(const size_t stateSize, double* rhs_vector, double* output) override;

protected:
    Eigen::SparseLU<Eigen::SparseMatrix<double>> m_solver;
};

}

#endif
