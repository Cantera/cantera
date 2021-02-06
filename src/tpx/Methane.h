//! @file Methane.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef TPX_METHANE_H
#define TPX_METHANE_H

#include "cantera/tpx/Sub.h"

namespace tpx
{

//! Pure species representation of methane. Values and functions are
//! from "Thermodynamic Properties in SI" by W.C. Reynolds
class methane : public Substance
{
public:
    methane() {
        m_name = "methane";
        m_formula = "CH4";
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

    //! Saturation pressure. Equation S3 from Reynolds TPSI.
    double Psat();

private:
    //! Liquid density. Equation D3 from Reynolds TPSI.
    double ldens();

    double C(int i, double rt, double rt2);
    double Cprime(int i, double rt, double rt2, double rt3);
    double I(int i, double egrho);
    double H(int i, double egrho);
    double W(int i, double egrho);
};
}
#endif // ! METHANE_H
