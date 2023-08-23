//! @file Water.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef TPX_WATER_H
#define TPX_WATER_H

#include "cantera/tpx/Sub.h"

namespace tpx
{
//! Pure species representation of water. Values and functions are from
//! from Reynolds @cite reynolds1979.
class water : public Substance
{
public:
    water() {
        m_name = "water";
        m_formula = "H2O";
    }

    double MolWt() override;
    double Tcrit() override;
    double Pcrit() override;
    double Vcrit() override;
    double Tmin() override;
    double Tmax() override;

    double Pp() override;
    double up() override;
    double sp() override;
    double Psat() override;
    double dPsatdT();

protected:
    double ldens() override;

private:
    double C(int i);
    double Cprime(int i);
    double I(int i);
    double H(int i);
};

}
#endif // ! WATER_H
