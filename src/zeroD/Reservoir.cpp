//! @file Reservoir.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/Reservoir.h"

namespace Cantera
{

shared_ptr<Reservoir> newReservoir(shared_ptr<Solution> contents, const string& name)
{
    return make_shared<Reservoir>(contents, name);
}

}
