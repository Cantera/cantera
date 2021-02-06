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
 */
class IdealGasConstPressureReactor : public ConstPressureReactor
{
public:
    IdealGasConstPressureReactor() {}

    virtual std::string typeStr() const {
        return "IdealGasConstPressureReactor";
    }

    /*!
     * @deprecated To be changed after Cantera 2.5.
     */
    virtual int type() const {
        warn_deprecated("IdealGasConstPressureReactor::type",
                        "To be changed after Cantera 2.5. "
                        "Return string instead of magic number; use "
                        "IdealGasConstPressureReactor::typeStr during "
                        "transition");
        return IdealGasConstPressureReactorType;
    }

    virtual void setThermoMgr(ThermoPhase& thermo);

    virtual void getState(doublereal* y);

    virtual void initialize(doublereal t0 = 0.0);
    virtual void evalEqs(doublereal t, doublereal* y,
                         doublereal* ydot, doublereal* params);

    virtual void updateState(doublereal* y);

    //! Return the index in the solution vector for this reactor of the
    //! component named *nm*. Possible values for *nm* are "mass",
    //! "temperature", the name of a homogeneous phase species, or the name of a
    //! surface species.
    virtual size_t componentIndex(const std::string& nm) const;
    std::string componentName(size_t k);

protected:
    vector_fp m_hk; //!< Species molar enthalpies
};
}

#endif
