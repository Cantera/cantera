//! @file Reservoir.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_RESERVOIR_H
#define CT_RESERVOIR_H

#include "ReactorBase.h"

namespace Cantera
{

class Reservoir : public ReactorBase
{
public:
    Reservoir() {}
    virtual int type() const {
        return ReservoirType;
    }
    virtual void initialize(doublereal t0 = 0.0) {}

    void insert(ThermoPhase& contents) {
        setThermoMgr(contents);
    }
};

}

#endif
