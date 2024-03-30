//! @file IdealGasReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IDEALGASREACTOR_H
#define CT_IDEALGASREACTOR_H

#include "Reactor.h"

namespace Cantera
{

/**
 * Class IdealGasReactor is a class for stirred reactors that is specifically
 * optimized for ideal gases. In this formulation, temperature replaces the
 * total internal energy as a state variable.
 * @ingroup reactorGroup
 */
class IdealGasReactor : public Reactor
{
public:
    using Reactor::Reactor; // inherit constructors

    string type() const override {
        return "IdealGasReactor";
    }

    void getState(double* y) override;

    void initialize(double t0=0.0) override;

    void eval(double t, double* LHS, double* RHS) override;

    void updateState(double* y) override;

    //! Return the index in the solution vector for this reactor of the
    //! component named *nm*. Possible values for *nm* are "mass",
    //! "volume", "temperature", the name of a homogeneous phase species, or the
    //! name of a surface species.
    size_t componentIndex(const string& nm) const override;
    string componentName(size_t k) override;

protected:
    void setThermo(ThermoPhase& thermo) override;

    vector<double> m_uk; //!< Species molar internal energies
};

}

#endif
