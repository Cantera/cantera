//! @file IdealGasMoleReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IDEALGASMOLE_REACTOR_H
#define CT_IDEALGASMOLE_REACTOR_H

#include "cantera/zeroD/MoleReactor.h"

namespace Cantera
{

/**
 * IdealGasMoleReactor is a class for ideal gas constant-volume reactors which use a
 * state of moles.
 * @since New in %Cantera 3.0
 * @ingroup reactorGroup
 */
class IdealGasMoleReactor : public MoleReactor
{
public:
    using MoleReactor::MoleReactor; // inherit constructors

    string type() const override {
        return "IdealGasMoleReactor";
    }

    size_t componentIndex(const string& nm) const override;
    string componentName(size_t k) override;
    double upperBound(size_t k) const override;
    double lowerBound(size_t k) const override;
    vector<size_t> steadyConstraints() const override;

    void getState(double* y) override;

    void initialize(double t0=0.0) override;

    void eval(double t, double* LHS, double* RHS) override;

    void updateState(double* y) override;

    //! Calculate the Jacobian to accelerate preconditioned solvers
    void buildJacobian(vector<Eigen::Triplet<double>>& jacVector) override;

    bool preconditionerSupported() const override {return true;};

    double temperatureDerivative() override { return 1; };

    double temperatureRadiationDerivative() override {
        return 4 * std::pow(temperature(), 3);
    }

    double moleDerivative(size_t index) override;

    double moleRadiationDerivative(size_t index) override;

    size_t speciesOffset() const override { return m_sidx; };

protected:
    vector<double> m_uk; //!< Species molar internal energies
};

}

#endif
