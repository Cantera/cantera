//! @file SoaveRedlichKwong.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/SoaveRedlichKwong.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Species.h"
#include "cantera/base/stringUtils.h"

#include <boost/math/tools/roots.hpp>

namespace bmt = boost::math::tools;

namespace Cantera
{
const double SoaveRedlichKwong::omega_a = 4.27480233540E-01;
const double SoaveRedlichKwong::omega_b = 8.66403499650E-02;
const double SoaveRedlichKwong::omega_vc = 3.33333333333333E-01;
}
