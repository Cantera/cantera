//! @file Hydrogen.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef TPX_HYDROGEN_H
#define TPX_HYDROGEN_H

#include "cantera/tpx/Sub.h"

namespace tpx
{

//! Pure species representation of hydrogen. Values and functions are
//! from "Thermodynamic Properties in SI" by W.C. Reynolds
class hydrogen : public Substance
{
public:
    hydrogen() {
        m_name = "hydrogen";
        m_formula = "H2";
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

    //! Saturation pressure. Equation s3 in Reynolds TPSI.
    double Psat();

private:
    //! Liquid density. Equation D4 in Reynolds TPSI.
    double ldens();
    double C(int i, double rt, double rt2);
    double Cprime(int i, double rt, double rt2, double rt3);
    double I(int i, double egrho);
    double H(int i, double egrho);
    double W(int i, double egrho);
    double icv(int i, double x, double xlg);
};

}

#endif // ! HYDROGEN_H
