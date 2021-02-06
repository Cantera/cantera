//! @file Water.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef TPX_WATER_H
#define TPX_WATER_H

#include "cantera/tpx/Sub.h"

namespace tpx
{
//! Pure species representation of water. Values and functions are from
//! "Thermodynamic Properties in SI" by W.C. Reynolds
class water : public Substance
{
public:
    water() {
        m_name = "water";
        m_formula = "H2O";
    }

    double MolWt();
    double Tcrit();
    double Pcrit();
    double Vcrit();
    double Tmin();
    double Tmax();

    double Pp();
    double up();
    double sp();
    double Psat();
    double dPsatdT();

private:
    double ldens();
    double C(int i);
    double Cprime(int i);
    double I(int i);
    double H(int i);
};

}
#endif // ! WATER_H
