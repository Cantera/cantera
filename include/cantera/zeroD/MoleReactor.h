//! @file MoleReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MOLEREACTOR_H
#define CT_MOLEREACTOR_H

#include "Reactor.h"

namespace Cantera
{

/*!
 * MoleReactor is meant to serve the same purpose as the reactor class but with a state
 * vector composed of moles. It also serves as the base class for other mole reactors.
 */
class MoleReactor : public Reactor
{
public:
    MoleReactor() {}

    virtual std::string type() const {
        return "MoleReactor";
    }

    virtual void initialize(double t0 = 0.0);

    virtual void getState(double* y);

    virtual void updateState(double* y);

    virtual void eval(double t, double* LHS, double* RHS);

    size_t componentIndex(const std::string& nm) const;

    std::string componentName(size_t k);

    //! Add the surface chemistry Jacobian values to m_jac_trips
    virtual void addSurfJacobian();

protected:
    //! Get moles of the system from mass fractions stored by thermo object
    //! @param y vector for moles to be put into
    virtual void getMoles(double* y);

    //! Set internal mass variable based on moles given
    //! @param y vector of moles of the system
    virtual void setMassFromMoles(double* y);

    virtual void evalSurfaces(double* LHS, double* RHS, double* sdot);

    virtual void updateSurfaceState(double* y);

    virtual void getSurfaceInitialConditions(double* y);

    //! const value for the species start index
    const size_t m_sidx = 2;
};

}

#endif
