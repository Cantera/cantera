//! @file utils.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/tpx/utils.h"
#include "cantera/base/stringUtils.h"

#include "CarbonDioxide.h"
#include "Heptane.h"
#include "HFC134a.h"
#include "Hydrogen.h"
#include "Methane.h"
#include "Nitrogen.h"
#include "Oxygen.h"
#include "Water.h"
#include "RedlichKwong.h"

namespace tpx
{
Substance* GetSubstanceByName(std::string name)
{
    std::string lcname = boost::algorithm::to_lower_copy(name);
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
    } else if (lcname == "hfc134a") {
        return new HFC134a;
    } else if (lcname == "rk") {
        return new RedlichKwong;
    } else if (lcname == "carbondioxide") {
        return new CarbonDioxide;
    } else if (lcname == "heptane") {
        return new Heptane;
    } else {
        return 0;
    }
}

Substance* GetSub(int isub)
{
    if (isub == 0) {
        return new water;
    } else if (isub == 1) {
        return new nitrogen;
    } else if (isub == 2) {
        return new methane;
    } else if (isub == 3) {
        return new hydrogen;
    } else if (isub == 4) {
        return new oxygen;
    } else if (isub == 5) {
        return new HFC134a;
    } else if (isub == 6) {
        return new RedlichKwong;
    } else if (isub == 7) {
        return new CarbonDioxide;
    } else if (isub == 8) {
        return new Heptane;
    } else {
        return 0;
    }
}

}
