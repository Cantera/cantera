//! @file Nitrogen.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef TPX_NITROGEN_H
#define TPX_NITROGEN_H

#include "cantera/tpx/Sub.h"

namespace tpx
{

//! Pure species representation of nitrogen. Values and functions are
//! from "Thermodynamic Properties in SI" by W.C. Reynolds
class nitrogen : public Substance
{
public:
    nitrogen() {
        m_name = "nitrogen";
        m_formula = "N2";
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

    //! Saturation pressure. Equation S4 from Reynolds TPSI.
    double Psat();

private:
    //! Liquid density. Equation D2 from Reynolds TPSI.
    double ldens();

    //! Equation P4 from Reynolds TPSI.
    double C(int i, double rt, double rt2);
    double Cprime(int i, double rt, double rt2, double rt3);
    double I(int i, double egrho);
    double H(int i, double egrho);
    double W(int i, double egrho);
};

}

#endif // ! TPX_NITROGEN_H

