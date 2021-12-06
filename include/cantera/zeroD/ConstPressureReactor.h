//! @file ConstPressureReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_CONSTP_REACTOR_H
#define CT_CONSTP_REACTOR_H

#include "Reactor.h"

namespace Cantera
{

/**
 * Class ConstPressureReactor is a class for constant-pressure reactors. The
 * reactor may have an arbitrary number of inlets and outlets, each of which
 * may be connected to a "flow device" such as a mass flow controller, a
 * pressure regulator, etc. Additional reactors may be connected to the other
 * end of the flow device, allowing construction of arbitrary reactor
 * networks.
 */
class ConstPressureReactor : public Reactor
{
public:
    ConstPressureReactor() {}

    virtual std::string typeStr() const {
        warn_deprecated("ConstPressureReactor::typeStr",
                        "To be removed after Cantera 2.6. Use type() instead.");
        return "ConstPressureReactor";
    }

    virtual std::string type() const {
        return "ConstPressureReactor";
    }

    virtual void getState(doublereal* y);

    virtual void initialize(doublereal t0 = 0.0);
    virtual void eval(double t, double* LHS, double* RHS);

    virtual void updateState(doublereal* y);

    //! Return the index in the solution vector for this reactor of the
    //! component named *nm*. Possible values for *nm* are "mass", "enthalpy",
    //! the name of a homogeneous phase species, or the name of a surface
    //! species.
    virtual size_t componentIndex(const std::string& nm) const;
    std::string componentName(size_t k);
};

}

#endif
