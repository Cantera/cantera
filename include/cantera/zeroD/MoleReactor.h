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
 * vector composed of moles. It is currently not functional and only serves as a base
 * class for other reactors.
 */
class MoleReactor : public Reactor
{
public:
    MoleReactor() {}

    virtual std::string type() const {
        return "MoleReactor";
    }

    virtual void initialize(double t0 = 0.0);

    virtual void eval(double t, double* LHS, double* RHS) {
        throw NotImplementedError("MoleReactor::eval()");
    }

    size_t componentIndex(const std::string& nm) const {
        throw NotImplementedError("MoleReactor::componentIndex");
    }

    std::string componentName(size_t k) {
        throw NotImplementedError("MoleReactor::componentName");
    }

protected:
    virtual void evalSurfaces(double* LHS, double* RHS, double* sdot);

    virtual void updateSurfaceState(double* y);

    virtual void getSurfaceInitialConditions(double* y);

    //! const value for the species start index
    const int m_sidx = 2;
};

}

#endif
