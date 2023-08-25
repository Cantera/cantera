//! @file utils.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/tpx/utils.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/global.h"

#include "CarbonDioxide.h"
#include "Heptane.h"
#include "HFC134a.h"
#include "Hydrogen.h"
#include "Methane.h"
#include "Nitrogen.h"
#include "Oxygen.h"
#include "Water.h"

namespace tpx
{
Substance* newSubstance(const std::string& name)
{
    std::string lcname = Cantera::toLowerCopy(name);
    if (lcname == "water") {
        return new water;
    } else if (lcname == "nitrogen") {
        return new nitrogen;
    } else if (lcname == "methane") {
        return new methane;
    } else if (lcname == "hydrogen") {
        return new hydrogen;
    } else if (lcname == "oxygen") {
        return new oxygen;
    } else if (lcname == "hfc-134a" || lcname == "hfc134a") {
        return new HFC134a;
    } else if (lcname == "carbon-dioxide" || lcname == "carbondioxide") {
        return new CarbonDioxide;
    } else if (lcname == "heptane") {
        return new Heptane;
    } else {
        throw Cantera::CanteraError("tpx::newSubstance", "No Substance"
            " definition known for '{}'.", name);
    }
}

}
