//! @file CarbonDioxide.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef TPX_CARBONDIOXIDE_H
#define TPX_CARBONDIOXIDE_H

#include "cantera/tpx/Sub.h"

namespace tpx
{

//! Pure species representation of carbon dioxide. Values and functions are
//! from "Thermodynamic Properties in SI" by W.C. Reynolds
class CarbonDioxide : public Substance
{
public:
    CarbonDioxide() {
        m_name="carbon-dioxide";
        m_formula="CO2";
    }

    double MolWt();
    double Tcrit();
    double Pcrit();
    double Vcrit();
    double Tmin();
    double Tmax();

    //! Pressure. Equation P-3 in Reynolds. P(rho, T).
    double Pp();

    /*!
     * internal energy. See Reynolds eqn (15) section 2
     *
     *  u = (the integral from T to To of co(T)dT) +
     *         sum from i to N ([C(i) - T*Cprime(i)] + uo
     */
    double up();

    //! entropy. See Reynolds eqn (16) section 2
    double sp();

    //! Pressure at Saturation. Equation S-2 in Reynolds.
    double Psat();

private:
    //! Liquid density. Equation D2 in Reynolds.
    double ldens();

    /*!
     * C returns a multiplier in each term of the sum in P-3, used in
     * conjunction with C in the function Pp
     * - j is used to represent which of the values in the summation to calculate
     * - j=0 is the second additive in the formula in reynolds
     * - j=1 is the third...
     * (this part does not include the multiplier rho^n)
     */
    double C(int jm, double, double, double, double);

    //! Derivative of C(i)
    double Cprime(int i, double, double, double);

    /*!
     * I = integral from o-rho { 1/(rho^2) * H(i, rho) d rho }
     * ( see section 2 of Reynolds TPSI )
     */
    double I(int i, double, double);

    /*!
     * H returns a multiplier in each term of the sum in P-3. This is used in
     * conjunction with C in the function Pp this represents the product
     * rho^n
     * - i=0 is the second additive in the formula in reynolds
     * - i=1 is the third ...
     */
    double H(int i, double egrho);
};

}

#endif // ! TPX_CARBONDIOXIDE_H
