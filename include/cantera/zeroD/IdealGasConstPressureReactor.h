//! @file ConstPressureReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IDEALGASCONSTP_REACTOR_H
#define CT_IDEALGASCONSTP_REACTOR_H

#include "ConstPressureReactor.h"

namespace Cantera
{

/**
 * Class ConstPressureReactor is a class for constant-pressure reactors. The
 * reactor may have an arbitrary number of inlets and outlets, each of which may
 * be connected to a "flow device" such as a mass flow controller, a pressure
 * regulator, etc. Additional reactors may be connected to the other end of the
 * flow device, allowing construction of arbitrary reactor networks.
 * @ingroup reactorGroup
 */
class IdealGasConstPressureReactor : public ConstPressureReactor
{
public:
    IdealGasConstPressureReactor() {}

    string type() const override {
        return "IdealGasConstPressureReactor";
    }

    void setThermoMgr(ThermoPhase& thermo) override;

    void getState(double* y) override;

    void initialize(double t0=0.0) override;
    void eval(double t, double* LHS, double* RHS) override;

    void updateState(double* y) override;

    //! Return the index in the solution vector for this reactor of the
    //! component named *nm*. Possible values for *nm* are "mass",
    //! "temperature", the name of a homogeneous phase species, or the name of a
    //! surface species.
    size_t componentIndex(const string& nm) const override;
    string componentName(size_t k) override;

protected:
    vector<double> m_hk; //!< Species molar enthalpies
};
}

#endif
