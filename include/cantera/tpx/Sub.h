//! @file Sub.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef TPX_SUB_H
#define TPX_SUB_H

#include "cantera/base/ctexceptions.h"
#include <algorithm>

namespace tpx
{

namespace PropertyPair
{
enum type {
    TV = 12, HP = 34, SP = 54, PV = 42, TP = 14, UV = 62, ST = 51,
    SV = 52, UP = 64, VH = 23, TH = 13, SH = 53, PX = 47, TX = 17,
    VT = -12, PH = -34, PS = -54, VP = -42, PT = -14, VU = -62, TS = -51,
    VS = -52, PU = -64, HV = -23, HT = -13, HS = -53, XP = -47, XT = -17
};
}

const int Pgiven = 0, Tgiven = 1;

namespace propertyFlag
{
enum type { H, S, U, V, P, T };
}

const double Undef = 999.1234;

/*!
 * Base class from which all pure substances are derived
 */
class Substance
{
public:
    Substance();

    virtual ~Substance() {}

    void setStdState(double h0 = 0.0, double s0 = 0.0,
                     double t0 = 298.15, double p0 = 1.01325e5);

    //! @name Information about a substance
    //! @{

    //! Molecular weight [kg/kmol]
    virtual double MolWt()=0;

    //! Critical temperature [K]
    virtual double Tcrit()=0;

    //! Critical pressure [Pa]
    virtual double Pcrit()=0;

    //! Critical specific volume [m^3/kg]
    virtual double Vcrit()=0;

    //! Minimum temperature for which the equation of state is valid
    virtual double Tmin()=0;

    //! Maximum temperature for which the equation of state is valid
    virtual double Tmax()=0;

    //! Name of the substance
    const char* name() {
        return m_name.c_str();
    }

    //! Chemical formula for the substance
     const char* formula() {
        return m_formula.c_str();
    }

    //! @}

    //! @name Properties
    //! @{

    //! Pressure [Pa]. If two phases are present, return the saturation
    //! pressure; otherwise return the pressure computed directly from the
    //! underlying eos.
    double P();

    //! Temperature [K]
    double Temp() {
        return T;
    }

    //! Specific volume [m^3/kg]
    double v() {
        return prop(propertyFlag::V);
    }

    //! Internal energy [J/kg]
    double u() {
        return prop(propertyFlag::U);
    }

    //! Enthalpy [J/kg]
    double h() {
        return prop(propertyFlag::H);
    }

    //! Entropy [J/kg/K]
    double s() {
        return prop(propertyFlag::S);
    }

    //! Helmholtz function [J/kg]
    double f() {
        return u() - T*s();
    }

    //! Gibbs function [J/kg]
    double g() {
        return h() - T*s();
    }

    //! Specific heat at constant volume [J/kg/K]
    virtual double cv();

    //! Specific heat at constant pressure [J/kg/K]
    virtual double cp();

    virtual double thermalExpansionCoeff();

    virtual double isothermalCompressibility();

    //! @}
    //! @name Saturation Properties
    //! @{

    double Ps();

    //! The derivative of the saturation pressure with respect to temperature.
    virtual double dPsdT();

    //! Saturation temperature at pressure *p*.
    double Tsat(double p);

    //! Vapor mass fraction. If T >= Tcrit, 0 is returned for v < Vcrit, and 1
    //! is returned if v > Vcrit.
    double x();

    //! Returns 1 if the current state is a liquid/vapor mixture, 0 otherwise
    int TwoPhase();
    //! @}

    virtual double Pp()=0;

    //! Enthaply of a single-phase state
    double hp() {
        return up() + Pp()/Rho;
    }

    //! Gibbs function of a single-phase state
    double gp() {
        return hp() - T*sp();
    }

    double prop(propertyFlag::type ijob);

    //! set T and P
    void set_TPp(double t0, double p0);

    //! Function to set or change the state for a property pair *XY* where
    //! *x0* is the value of first property and *y0* is the value of the
    //! second property.
    void Set(PropertyPair::type XY, double x0, double y0);

protected:
    double T, Rho;
    double Tslast, Rhf, Rhv;
    double Pst;
    double m_energy_offset;
    double m_entropy_offset;
    std::string m_name;
    std::string m_formula;

    virtual double ldens()=0;

    //! Saturation pressure, Pa
    virtual double Psat()=0;

    //! Internal energy of a single-phase state
    virtual double up()=0;

    //! Entropy of a single-phase state
    virtual double sp()=0;

    virtual int ideal() {
        return 0;
    }

    double vp() {
        return 1.0/Rho;
    }

    //! Uses the lever rule to set state in the dome. Returns 1 if in dome,
    //! 0 if not, in which case state not set.
    int Lever(int itp, double sat, double val, propertyFlag::type ifunc);

    //! Update saturated liquid and vapor densities and saturation pressure
    void update_sat();

private:
    void set_Rho(double r0);
    void set_T(double t0);
    void set_v(double v0);
    void BracketSlope(double p);
    double vprop(propertyFlag::type ijob);
    void set_xy(propertyFlag::type if1, propertyFlag::type if2,
                double X, double Y,
                double atx, double aty, double rtx, double rty);

    int kbr;
    double Vmin, Vmax;
    double Pmin, Pmax;
    double dvbf, dv;
    double v_here, P_here;
};

}

#endif
