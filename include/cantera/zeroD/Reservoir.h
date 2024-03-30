//! @file Reservoir.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_RESERVOIR_H
#define CT_RESERVOIR_H

#include "ReactorBase.h"
#include "cantera/base/Solution.h"

namespace Cantera
{

//! A source or sink whose state remains constant regardless of any flows or other
//! interactions with other Reactor objects.
//! @ingroup reactorGroup
class Reservoir : public ReactorBase
{
public:
    using ReactorBase::ReactorBase; // inherit constructors

    string type() const override {
        return "Reservoir";
    }

    void initialize(double t0=0.0) override {}

    //! @deprecated Unused; to be removed after %Cantera 3.1.
    void insert(ThermoPhase& contents) {
        warn_deprecated("Reservoir::insert",
            "Unused; to be removed after Cantera 3.1.");
        setThermo(contents);
    }

    using ReactorBase::insert;
};

}

#endif
