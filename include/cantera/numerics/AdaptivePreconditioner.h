/**
 *  @file AdaptivePreconditioner.h Declarations for the class
 *   AdaptivePreconditioner which is a child class of SystemJacobian
 *   for preconditioners used by sundials
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef ADAPTIVEPRECONDITIONER_H
#define ADAPTIVEPRECONDITIONER_H

#include "cantera/numerics/EigenSparseJacobian.h"

namespace Cantera
{

//! AdaptivePreconditioner a preconditioner designed for use with large
//! mechanisms that leverages sparse solvers. It does this by pruning
//! the preconditioner by a threshold value. It also neglects pressure
//! dependence and third body contributions in its formation and has a
//! finite difference approximation for temperature.
class AdaptivePreconditioner : public EigenSparseJacobian
{
public:
    AdaptivePreconditioner();

    void initialize(size_t networkSize) override;
    const string type() const override { return "Adaptive"; }
    void factorize() override;
    void solve(const size_t stateSize, double* rhs_vector, double* output) override;
    void stateAdjustment(vector<double>& state) override;

    int info() const override {
        return static_cast<int>(m_solver.info());
    }

    //! Prune preconditioner elements
    void prunePreconditioner();

    //! Get the threshold value for setting elements
    double threshold() { return m_threshold; }

    //! Get ILUT fill factor
    double ilutFillFactor() { return m_fill_factor; }

    //! Get ILUT drop tolerance
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

protected:
    //! ILUT fill factor
    double m_fill_factor = 0;

    //! ILUT drop tolerance
    double m_drop_tol = 0;

    //! Solver used in solving the linear system
    Eigen::IncompleteLUT<double> m_solver;

    //! Minimum value a non-diagonal element must be to be included in
    //! the preconditioner
    double m_threshold = 0.0;

    //! Bool set whether to prune the matrix or not
    double m_prune_precon = true;
};

}

#endif
