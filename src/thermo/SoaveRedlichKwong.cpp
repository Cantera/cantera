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



int SoaveRedlichKwong::solveCubic(double T, double pres, double a, double b, double aAlpha,
                                  double Vroot[3]) const
{
    double an = 1.0;
    double bn = - GasConstant * T / pres;
    double cn = (aAlpha - b * GasConstant * T) / pres - b * b;
    double dn = aAlpha * b / pres;

    double tc = a * omega_b / (b * omega_a * GasConstant);
    double pc = omega_b * R * tc / b;
    double vc = omega_vc * GasConstant * tc / pc

    return MixtureFugacityTP::solveCubic(T, pres, a, b, aAlpha, Vroot,
                                         an, bn, cn, dn, tc, vc);
}

}
