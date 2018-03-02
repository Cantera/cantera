//! @file ConstPressureReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

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

    virtual int type() const {
        return ConstPressureReactorType;
    }

    virtual void getState(double* y);

    virtual void initialize(double t0 = 0.0);
    virtual void evalEqs(double t, double* y,
                         double* ydot, double* params);

    virtual void updateState(double* y);

    //! Return the index in the solution vector for this reactor of the
    //! component named *nm*. Possible values for *nm* are "mass", "enthalpy",
    //! the name of a homogeneous phase species, or the name of a surface
    //! species.
    virtual size_t componentIndex(const std::string& nm) const;
    std::string componentName(size_t k);
};

}

#endif
