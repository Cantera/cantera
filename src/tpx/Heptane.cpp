/**
 * @file Heptane.cpp representation of substance Heptane.
 *
 * Values and functions are from "Thermodynamic Properties in SI" by W.C.
 * Reynolds. AUTHOR: jrh@stanford.edu: GCEP, Stanford University
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "Heptane.h"
#include "cantera/base/stringUtils.h"

using namespace Cantera;

namespace tpx
{

// Heptane constants
static const double Tmn = 182.56; // [K] minimum temperature for which calculations are valid
static const double Tmx = 1000.0; // [K] maximum temperature for which calculations are valid
static const double Tc=537.68; // [K] critical temperature
static const double Roc=197.60; // [kg/m^3] critical density
static const double To=300; // [K] reference Temperature
static const double R=82.99504; // [J/(kg*K)] gas constant (for this substance)
static const double Gamma=9.611604E-6; // [??]
static const double u0=3.4058439E5; // [] internal energy at To
static const double s0=1.1080254E3; // [] entropy at To
static const double Tp=400; // [K] ??
static const double Pc=2.6199E6; // [Pa] critical pressure
static const double M=100.20; // [kg/kmol] molar density

// array Ahept is used by the function Pp
static const double Ahept[]= {
    2.246032E-3,
    2.082990E2,
    5.085746E7,
    3.566396E9,
    1.622168E9,
    1.065237E-5,
    5.987922E-1,
    7.736602,
    1.929386E5,
    5.291379E-9
};

// array F is used by Psat
static const double F[]= {
    -7.2298764,
    3.8607475E-1,
    -3.4216472,
    4.6274432E-1,
    -9.7926124,
    -4.2058094E1,
    7.5468678E1,
    3.1758992E2
};

// array D is used by the function ldens
static const double D[]= {
    1.9760405E2,
    8.9451237E2,
    -1.1462908E3,
    1.7996947E3,
    -1.7250843E3,
    9.7088329E2
};

// array G is used by the function sp
static const double G[]= {
    1.1925213E5,
    -7.7231363E2,
    7.4463527,
    -3.0888167E-3,
    0.0,
    0.0
};

double Heptane::C(int j,double Tinverse, double T2inverse, double T3inverse, double T4inverse)
{
    switch (j) {
    case 0:
        return Ahept[0] * R * T -
               Ahept[1] -
               Ahept[2] * T2inverse +
               Ahept[3] * T3inverse -
               Ahept[4] * T4inverse;
    case 1:
        return Ahept[5] * R * T -
               Ahept[6] -
               Ahept[7] * Tinverse;
    case 2:
        return Ahept[9] * (Ahept[6] + Ahept[7] * Tinverse);
    case 3:
        return Ahept[8] * T2inverse;
    default:
        return 0.0;
    }
}

double Heptane::Cprime(int j, double T2inverse, double T3inverse, double T4inverse)
{
    switch (j) {
    case 0:
        return Ahept[0] * R -
               -2 * Ahept[2] * T3inverse +
               -3 * Ahept[3] * T4inverse -
               -4 * Ahept[4] * pow(T, -5.0);
    case 1:
        return Ahept[5] * R -
               -1 * Ahept[7] * T2inverse;
    case 2:
        return Ahept[9] * (-1 * Ahept[7] * T2inverse);
    case 3:
        return -2 * Ahept[8] * T3inverse;
    default:
        return 0.0;
    }
}

double Heptane::I(int j, double ergho, double Gamma)
{
    switch (j) {
    case 0:
        return Rho;
    case 1:
        return Rho * Rho / 2;
    case 2:
        return pow(Rho, 5.0)/ 5;
    case 3:
        return 1 / Gamma - (Gamma * Rho * Rho + 2) * ergho / (2 * Gamma);
    default:
        return 0.0;
    }
}

double Heptane::H(int i, double egrho)
{
    if (i < 2) {
        return pow(Rho,i+2);
    } else if (i == 2) {
        return pow(Rho,6.0);
    } else if (i == 3) {
        return pow(Rho,3) * (1 + Gamma * Rho * Rho) * egrho;
    } else {
        return 0;
    }
}

double Heptane::up()
{
    double Tinverse = 1.0/T;
    double T2inverse = pow(T, -2);
    double T3inverse = pow(T, -3);
    double T4inverse = pow(T, -4);
    double egrho = exp(-Gamma*Rho*Rho);

    double sum = 0.0;
    int i;
    for (i=1; i<=5; i++) {
        sum += G[i]*(pow(T,i) - pow(To,i))/double(i);
    }
    sum += G[0]*log(T/To);
    for (i=0; i<=6; i++) {
        sum += (C(i, Tinverse, T2inverse, T3inverse, T4inverse) - T*Cprime(i,T2inverse, T3inverse, T4inverse))*I(i,egrho, Gamma);
    }
    sum += u0;
    return sum + m_energy_offset;
}

double Heptane::sp()
{
    double T2inverse = pow(T, -2);
    double T3inverse = pow(T, -3);
    double T4inverse = pow(T, -4);
    double egrho = exp(-Gamma*Rho*Rho);

    double sum = 0.0;
    for (int i=2; i<=5; i++) {
        sum += G[i]*(pow(T,i-1) - pow(To,i-1))/double(i-1);
    }
    sum += G[1]*log(T/To);
    sum -= G[0]*(1.0/T - 1.0/To);
    for (int i=0; i<=6; i++) {
        sum -= Cprime(i,T2inverse, T3inverse, T4inverse)*I(i,egrho, Gamma);
    }
    sum += s0 - R*log(Rho);
    return sum + m_entropy_offset;
}

double Heptane::Pp()
{
    double Tinverse = pow(T,-1);
    double T2inverse = pow(T, -2);
    double T3inverse = pow(T, -3);
    double T4inverse = pow(T, -4);
    double egrho = exp(-Gamma*Rho*Rho);

    double P = Rho*R*T;
    for (int i=0; i<=3; i++) {
        P += C(i,Tinverse, T2inverse, T3inverse, T4inverse)*H(i,egrho);
    }
    return P;
}

double Heptane::Psat()
{
    double log, sum=0;
    if ((T < Tmn) || (T > Tc)) {
        throw CanteraError("Heptane::Psat",
                           "Temperature out of range. T = {}", T);
    }
    for (int i=1; i<=8; i++) {
        sum += F[i-1] * pow((T/Tp -1),double(i-1));
    }

    log = ((Tc/T)-1)*sum;
    return exp(log)*Pc;
}

double Heptane::ldens()
{
    double xx=1-(T/Tc), sum=0;
    if ((T < Tmn) || (T > Tc)) {
        throw CanteraError("Heptane::ldens",
                           "Temperature out of range. T = {}", T);
    }
    for (int i=1; i<=6; i++) {
        sum+=D[i-1]*pow(xx,double(i-1)/3.0);
    }

    return sum;
}

// The following functions allow users to get the properties of Heptane that
// are not dependent on the state

double Heptane::Tcrit()
{
    return Tc;
}
double Heptane::Pcrit()
{
    return Pc;
}
double Heptane::Vcrit()
{
    return 1.0/Roc;
}
double Heptane::Tmin()
{
    return Tmn;
}
double Heptane::Tmax()
{
    return Tmx;
}
double Heptane::MolWt()
{
    return M;
}
}
