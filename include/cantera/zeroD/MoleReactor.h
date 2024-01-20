//! @file MoleReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MOLEREACTOR_H
#define CT_MOLEREACTOR_H

#include "Reactor.h"

namespace Cantera
{

/**
 * MoleReactor is meant to serve the same purpose as the reactor class but with a state
 * vector composed of moles. It also serves as the base class for other mole reactors.
 * @since New in %Cantera 3.0
 * @ingroup reactorGroup
 */
class MoleReactor : public Reactor
{
public:
    using Reactor::Reactor; // inherit constructors

    string type() const override {
        return "MoleReactor";
    }

    void initialize(double t0=0.0) override;

    void getState(double* y) override;

    void updateState(double* y) override;

    void eval(double t, double* LHS, double* RHS) override;

    size_t componentIndex(const string& nm) const override;

    string componentName(size_t k) override;

protected:
    //! For each surface in the reactor, update vector of triplets with all relevant
    //! surface jacobian derivatives of species with respect to species
    //! which are appropriately offset to align with the reactor's state vector.
    virtual void addSurfaceJacobian(vector<Eigen::Triplet<double>> &triplets);

    //! Get moles of the system from mass fractions stored by thermo object
    //! @param y vector for moles to be put into
    void getMoles(double* y);

    //! Set internal mass variable based on moles given
    //! @param y vector of moles of the system
    void setMassFromMoles(double* y);

    void evalSurfaces(double* LHS, double* RHS, double* sdot) override;

    void updateSurfaceState(double* y) override;

    void getSurfaceInitialConditions(double* y) override;

    //! const value for the species start index
    const size_t m_sidx = 2;
};

}

#endif
