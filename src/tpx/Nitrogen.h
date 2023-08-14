//! @file Nitrogen.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef TPX_NITROGEN_H
#define TPX_NITROGEN_H

#include "cantera/tpx/Sub.h"

namespace tpx
{

//! Pure species representation of nitrogen. Values and functions are
//! from Reynolds @cite reynolds1979.
class nitrogen : public Substance
{
public:
    nitrogen() {
        m_name = "nitrogen";
        m_formula = "N2";
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

    //! Saturation pressure. Equation S4 from Reynolds TPSI.
    double Psat() override;

protected:
    //! Liquid density. Equation D2 from Reynolds TPSI.
    double ldens() override;

private:
    //! Equation P4 from Reynolds TPSI.
    double C(int i, double rt, double rt2);
    double Cprime(int i, double rt, double rt2, double rt3);
    double I(int i, double egrho);
    double H(int i, double egrho);
    double W(int i, double egrho);
};

}

#endif // ! TPX_NITROGEN_H
