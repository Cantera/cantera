//! @file IdealGasConstPressureMoleReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IDEALGASCONSTPRESSMOLE_REACTOR_H
#define CT_IDEALGASCONSTPRESSMOLE_REACTOR_H

#include "cantera/zeroD/ConstPressureMoleReactor.h"

namespace Cantera
{

/**
 * IdealGasConstPressureMoleReactor is a class for ideal gas constant-pressure reactors
 * which use a state of moles.
 * @since New in %Cantera 3.0
 * @ingroup reactorGroup
 */
class IdealGasConstPressureMoleReactor : public ConstPressureMoleReactor
{
public:
    using ConstPressureMoleReactor::ConstPressureMoleReactor; // inherit constructors

    string type() const override {
        return "IdealGasConstPressureMoleReactor";
    };

    size_t componentIndex(const string& nm) const override;

    string componentName(size_t k) override;

    void getState(double* y) override;

    void initialize(double t0=0.0) override;

    void eval(double t, double* LHS, double* RHS) override;

    void updateState(double* y) override;

    //! Calculate an approximate Jacobian to accelerate preconditioned solvers

    //! Neglects derivatives with respect to mole fractions that would generate a
    //! fully-dense Jacobian. Currently also neglects terms related to interactions
    //! between reactors, for example via inlets and outlets.
    Eigen::SparseMatrix<double> jacobian() override;

    bool preconditionerSupported() const override { return true; };

protected:
    void setThermo(ThermoPhase& thermo) override;

    vector<double> m_hk; //!< Species molar enthalpies
};

}

#endif
