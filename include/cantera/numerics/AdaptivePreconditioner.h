/**
 *  @file AdaptivePreconditioner.h Declarations for the class
 *   AdaptivePreconditioner which is a child class of PreconditionerBase
 *   for preconditioners used by sundials
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef ADAPTIVEPRECONDITIONER_H
#define ADAPTIVEPRECONDITIONER_H

#include "cantera/numerics/PreconditionerBase.h"
#include "cantera/numerics/eigen_sparse.h"
#include "cantera/base/global.h"
#include <iostream>

namespace Cantera
{

//! AdaptivePreconditioner a preconditioner designed for use with large
//! mechanisms that leverages sparse solvers. It does this by pruning
//! the preconditioner by a threshold value. It also neglects pressure
//! dependence and thirdbody contributions in its formation and has a
//! finite difference approximation for temperature.
class AdaptivePreconditioner : public PreconditionerBase
{
public:
    AdaptivePreconditioner() {}

    void initialize(size_t networkSize);

    void reset() {
        m_precon_matrix.setZero();
        m_jac_trips.clear();
    };

    void setup();

    void solve(const size_t stateSize, double* rhs_vector, double* output);

    PreconditionerType preconditionerType() {
        return PreconditionerType::LEFT_PRECONDITION;
    }

    void setValue(size_t row, size_t col, double value);

    virtual void stateAdjustment(vector_fp& state);

    virtual void updatePreconditioner();

    //! Prune preconditioner elements
    void prunePreconditioner();

    //! Function used to return semi-analytical jacobian matrix
    Eigen::SparseMatrix<double> jacobian() {
        Eigen::SparseMatrix<double> jacobian_mat(m_dim, m_dim);
        jacobian_mat.setFromTriplets(m_jac_trips.begin(), m_jac_trips.end());
        return jacobian_mat;
    }

    //! Return the internal preconditioner matrix
    Eigen::SparseMatrix<double> matrix() {
        updatePreconditioner();
        return m_precon_matrix;
    }

    //! Get the threshold value for setting elements
    double threshold() { return m_threshold; }

    //! Get ilut fill factor
    double ilutFillFactor() { return m_fill_factor; }

    //! Get ilut drop tolerance
    double ilutDropTol() { return m_drop_tol; }

    //! Set the threshold value to compare elements against
    //! @param threshold double value used in setting by threshold
    void setThreshold(double threshold) {
        m_threshold = threshold;
        m_prune_precon = (threshold <= 0) ? false : true;
    }

    //! Set drop tolerance for ILUT
    //! @param droptol double value used in setting solver drop tolerance
    void setIlutDropTol(double droptol) {
        m_drop_tol = droptol;
        m_solver.setDroptol(droptol);
        }

    //! Set the fill factor for ILUT solver
    //! @param fillFactor fill in factor for ILUT solver
    void setIlutFillFactor(int fillFactor) {
        m_fill_factor = fillFactor;
        m_solver.setFillfactor(fillFactor);
    }

    //! Print preconditioner contents
    void printPreconditioner();

    //! Print jacobian contents
    void printJacobian();

protected:
    //! ilut fill factor
    double m_fill_factor = 0;

    //! ilut drop tolerance
    double m_drop_tol = 0;

    //! Vector of triples representing the jacobian used in preconditioning
    std::vector<Eigen::Triplet<double>> m_jac_trips;

    //! Storage of appropriately sized identity matrix for making the preconditioner
    Eigen::SparseMatrix<double> m_identity;

    //! Container that is the sparse preconditioner
    Eigen::SparseMatrix<double> m_precon_matrix;

    //! Solver used in solving the linear system
    Eigen::IncompleteLUT<double> m_solver;

    //! Minimum value a non-diagonal element must be to be included in
    //! the preconditioner
    double m_threshold = 1e-8;

    //! Bool set whether to prune the matrix or not
    double m_prune_precon = true;
};

}

#endif
