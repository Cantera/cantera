/*
 * This is the base substance class from which all substances are derived
 *
 * Kate Talmazan: SURF -- July, 1995
 *   original implementation of this class and all derived classes from
 *   formulas given in TPSI. Implementation of P(Rho, T), cv0(T), ldens(T),
 *   and Psat(T) for all substances in TPSI.f
 *
 * Dave Goodwin: Fall, 1996
 *   functions for u, h, s, f, g;
 *   functions to set state
 *   error handling
 *   documentation
 *
 *   Sept., 2001: minor modifications to use with Cantera
 *
 */

#ifndef TPX_SUB_H
#define TPX_SUB_H

#include "cantera/base/ctexceptions.h"

#include <iostream>
#include <string>

namespace tpx
{

class TPX_Error : public Cantera::CanteraError
{
public:
    TPX_Error(std::string p, std::string e) : CanteraError(p, e) { }
};


const double OneAtm = 1.01325e5;
const double Liquid = 0.0;
const double Vapor = 1.0;

const int TV = 12, HP = 34, SP = 54, PV = 42, TP = 14, UV = 62,
          ST = 51, SV = 52, UP = 64, VH = 23, TH = 13, SH = 53,
          PX = 47, TX = 17;

const int VT = -12, PH = -34, PS = -54, VP = -42, PT = -14, VU = -62,
          TS = -51, VS = -52, PU = -64, HV = -23, HT = -13, HS = -53,
          XP = -47, XT = -17;

const int Pgiven = 0, Tgiven = 1;

const int EvalH = 3;
const int EvalS = 5;
const int EvalU = 6;
const int EvalV = 2;
const int EvalP = 4;
const int EvalT = 1;
const int EvalX = 7;
const int EvalF = 8;
const int EvalG = 9;
const int EvalC = 10;
const int EvalM = 11;
const int EvalN = 12;
const int EvalMW = 13;
const int EvalEA = 14;
const int EvalCdot = 15;
const int EvalDdot = 16;
const int EvalWdot = 17;
const int EvalTchem = 18;
const int EvalRgas = 19;

const double Undef = 999.1234;

class Substance
{
public:
    Substance();

    virtual ~Substance() {}

    void setStdState(double h0 = 0.0, double s0 = 0.0,
                     double t0 = 298.15, double p0 = 1.01325e5) {
        Set(TP, t0, p0);
        double hh = h();
        double ss = s();
        double hoff = h0 - hh;
        double soff = s0 - ss;
        m_entropy_offset += soff;
        m_energy_offset += hoff;
    }

    // information about a substance:

    virtual double MolWt()=0;          // molecular weight, kg/kmol
    virtual double Tcrit()=0;          // critical temperature, K
    virtual double Pcrit()=0;          // critical pressure, Pa
    virtual double Vcrit()=0;          // critical specific vol, m^3/kg
    virtual double Tmin()=0;           // min. temp for which equations valid
    virtual double Tmax()=0;           // max. temp for which equations valid
    virtual char* name() = 0;          // name
    virtual char* formula() = 0;       // chemical formula

    // properties:

    double P();                        // pressure, Pa
    double Temp() {
        return T;   // temperature, K
    }
    double v() {                       // specific vol, m^3/kg
        return prop(EvalV);
    }
    double u() {                       // int. energy, J/kg
        return prop(EvalU);
    }
    double h() {                       // enthalpy, J/kg
        return prop(EvalH);
    }
    double s() {                       // entropy, J/kg/K
        return prop(EvalS);
    }
    double f() {                       // Helmholtz function, J/kg
        return u() - T*s();
    }
    double g() {                       // Gibbs function, J/kg
        return h() - T*s();
    }

    virtual double cv() {
        double Tsave = T, dt = 1.e-4*T;
        set_T(Tsave - dt);
        double s1 = s();
        set_T(Tsave + dt);
        double s2 = s();
        set_T(Tsave);
        return T*(s2 - s1)/(2.0*dt);
    }

    virtual double cp() {
        double Tsave = T, dt = 1.e-4*T;
        double p0 = P();
        Set(TP, Tsave - dt, p0);
        double s1 = s();
        Set(TP, Tsave + dt, p0);
        double s2 = s();
        Set(TP, Tsave, p0);
        return T*(s2 - s1)/(2.0*dt);
    }

    virtual double thermalExpansionCoeff() {
        double Tsave = T, dt = 1.e-4*T;
        double p0 = P();
        Set(TP, Tsave - dt, p0);
        double v1 = v();
        Set(TP, Tsave + dt, p0);
        double v2 = v();
        Set(TP, Tsave, p0);
        return (v2 - v1)/((v2 + v1)*dt);
    }

    virtual double isothermalCompressibility() {
        double Psave = P(), dp = 1.e-4*Psave;
        Set(TP, T, Psave - dp);
        double v1 = v();
        Set(TP, T, Psave + dp);
        double v2 = v();
        Set(TP, T, Psave);
        return -(v2 - v1)/((v2 + v1)*dp);
    }


    // saturation properties

    double Ps();
    virtual double dPsdT();          // d(Psat)/dT, Pa/K
    double Tsat(double p);             // saturation temp at p
    double x();                        // vapor mass fraction
    int TwoPhase();                    // =1 if vapor/liquid, 0 otherwise
    virtual double Pp()=0;
    double hp() {
        return up() + Pp()/Rho;
    }
    double gp() {
        return hp() - T*sp();
    }

    double prop(int ijob);
    void set_TPp(double t0, double p0);    // set T and P


    // functions to set or change state:

    void Set(int XY, double x0, double y0);
    void Set_meta(double phase, double pp);

protected:

    double T, Rho;
    double Tslast, Rhf, Rhv;
    double Pst;
    double m_energy_offset;
    double m_entropy_offset;
    std::string m_name;
    std::string m_formula;

    //virtual double Xm(int k) { return 1.0;}
    //virtual int Species() { return 1;}

    virtual double ldens()=0;
    virtual double Psat()=0;           // saturation pressure, Pa
    virtual double up()=0;
    virtual double sp()=0;
    virtual int ideal() {
        return 0;   // added 9/2/98; default is false
    }
    double vp() {
        return 1.0/Rho;
    }
    int Lever(int itp, double sat, double val, int ifunc);
    void update_sat();

private:
    void set_Rho(double r0);
    void set_T(double t0);
    void set_v(double v0);
    void BracketSlope(double p);
    double lprop(int ijob);
    double vprop(int ijob);
    void set_xy(int if1, int if2, double X, double Y,
                double atx, double aty, double rtx, double rty);

    int kbr;
    double Vmin, Vmax;
    double Pmin, Pmax;
    double dvbf, dv;
    double v_here, P_here;
};

}


#endif
