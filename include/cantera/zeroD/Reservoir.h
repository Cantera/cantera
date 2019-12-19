//! @file Reservoir.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_RESERVOIR_H
#define CT_RESERVOIR_H

#include "ReactorBase.h"

namespace Cantera
{

class Reservoir : public ReactorBase
{
public:
    Reservoir() {}

    virtual std::string typeStr() const {
        return "Reservoir";
    }

    /*!
     * @deprecated To be changed after Cantera 2.5.
     */
    virtual int type() const {
        warn_deprecated("Reservoir::type",
                        "To be changed after Cantera 2.5. "
                        "Return string instead of magic number; use "
                        "Reservoir::typeStr during transition");
        return ReservoirType;
    }
    virtual void initialize(doublereal t0 = 0.0) {}

    void insert(ThermoPhase& contents) {
        setThermoMgr(contents);
    }

    void insert(shared_ptr<Solution> sol) {
        setThermoMgr(*sol->thermo());
    }
};

}

#endif
