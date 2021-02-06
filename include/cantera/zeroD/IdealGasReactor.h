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
 */
class IdealGasReactor : public Reactor
{
public:
    IdealGasReactor() {}

    virtual std::string typeStr() const {
        return "IdealGasReactor";
    }

    /*!
     * @deprecated To be changed after Cantera 2.5.
     */
    virtual int type() const {
        warn_deprecated("IdealGasReactor::type",
                        "To be changed after Cantera 2.5. "
                        "Return string instead of magic number; use "
                        "IdealGasReactor::typeStr during transition");
        return IdealGasReactorType;
    }

    virtual void setThermoMgr(ThermoPhase& thermo);

    virtual void getState(doublereal* y);

    virtual void initialize(doublereal t0 = 0.0);

    virtual void evalEqs(doublereal t, doublereal* y,
                         doublereal* ydot, doublereal* params);

    virtual void updateState(doublereal* y);

    //! Return the index in the solution vector for this reactor of the
    //! component named *nm*. Possible values for *nm* are "mass",
    //! "volume", "temperature", the name of a homogeneous phase species, or the
    //! name of a surface species.
    virtual size_t componentIndex(const std::string& nm) const;
    std::string componentName(size_t k);

protected:
    vector_fp m_uk; //!< Species molar internal energies
};

}

#endif
