//! @file Hydrogen.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "Hydrogen.h"
#include "cantera/base/stringUtils.h"

using namespace Cantera;

namespace tpx
{
static const double
M = 2.0159,
Tmn = 13.8,
Tmx = 5000.0,
Tc = 32.938,
Pc = 1.2838e6,
Roc= 31.36,
To = 13.8,
Tt = 13.8,
Pt = 7042.09,
R = 4124.299539,
Gamma = 1.008854772e-3,
u0 = 3.9275114e5,
s0 = 2.3900333e4,
T1 = 35,
T2 = 400,
alpha = 1.5814454428, //to be used with psat
alpha1 = .3479; //to be used with ldens

static const double Ahydro[] = {
    1.150470519352900e1, 1.055427998826072e3, -1.270685949968568e4,
    7.287844527295619e4, -7.448780703363973e5, 2.328994151810363e-1,
    -1.635308393739296e1, 3.730678064960389e3, 6.299667723184813e5,
    1.210920358305697e-3, 1.753651095884817, -1.367022988058101e2,
    -6.869936641299885e-3, 3.644494201750974e-2, -2.559784772600182,
    -4.038855202905836e-4, 1.485396303520942e-6, 4.243613981060742e-4,
    -2.307910113586888e-6, -6.082192173879582e5, -1.961080967486886e6,
    -5.786932854076408e2, 2.799129504191752e4, -2.381566558300913e-1,
    8.918796032452872e-1, -6.985739539036644e-5, -7.339554179182899e-3,
    -5.597033440289980e-9, 8.842130160884514e-8, -2.655507264539047e-12,
    -4.544474518140164e-12, 9.818775257001922e-11
};
static const double Dhydro[]= {
    4.8645813003e1, -3.4779278180e1, 4.0776538192e2,
    -1.1719787304e3, 1.6213924400e3, -1.1531096683e3, 3.3825492039e2
};
static const double Fhydro[]=
{ 3.05300134164, 2.80810925813, -6.55461216567e-1, 1.59514439374 };
static const double Ghydro[]= {
    6.1934792e3, 2.9490437e2, -1.5401979e3, -4.9176101e3,
    6.8957165e4, -2.2282185e5, 3.7990059e5, -3.7094216e5, 2.1326792e5,
    -7.1519411e4, 1.2971743e4, -9.8533014e2, 1.0434776e4,
    -3.9144179e2, 5.8277696e2, 6.5409163e2, -1.8728847e2
};

double hydrogen::C(int i, double rt, double rt2)
{
    switch (i) {
    case 0:
        return Ahydro[0] * T + Ahydro[1] * sqrt(T) + Ahydro[2] + (Ahydro[3] + Ahydro[4] * rt) * rt;
    case 1:
        return Ahydro[5] * T + Ahydro[6] + rt * (Ahydro[7] + Ahydro[8] * rt);
    case 2:
        return Ahydro[9] * T + Ahydro[10] + Ahydro[11] * rt;
    case 3:
        return Ahydro[12];
    case 4:
        return rt*(Ahydro[13] + Ahydro[14]*rt);
    case 5:
        return Ahydro[15]*rt;
    case 6:
        return rt*(Ahydro[16] + Ahydro[17]*rt);
    case 7:
        return Ahydro[18]*rt2;
    case 8:
        return rt2*(Ahydro[19] + Ahydro[20]*rt);
    case 9:
        return rt2*(Ahydro[21] + Ahydro[22]*rt2);
    case 10:
        return rt2*(Ahydro[23] + Ahydro[24]*rt);
    case 11:
        return rt2*(Ahydro[25] + Ahydro[26]*rt2);
    case 12:
        return rt2*(Ahydro[27] + Ahydro[28]*rt);
    case 13:
        return rt2*(Ahydro[29] + Ahydro[30]*rt + Ahydro[31]*rt2);
    default:
        return 0.0;
    }
}

double hydrogen::Cprime(int i, double rt, double rt2, double rt3)
{
    switch (i) {
    case 0:
        return Ahydro[0] + 0.5*Ahydro[1]/sqrt(T) - (Ahydro[3] + 2.0*Ahydro[4]*rt)*rt2;
    case 1:
        return Ahydro[5] - rt2*(Ahydro[7] + 2.0*Ahydro[8]*rt);
    case 2:
        return Ahydro[9] - Ahydro[11]*rt2;
    case 3:
        return 0.0;
    case 4:
        return -rt2*(Ahydro[13] + 2.0*Ahydro[14]*rt);
    case 5:
        return -Ahydro[15]*rt2;
    case 6:
        return -rt2*(Ahydro[16] + 2.0*Ahydro[17]*rt);
    case 7:
        return -2.0*Ahydro[18]*rt3;
    case 8:
        return -rt3*(2.0*Ahydro[19] + 3.0*Ahydro[20]*rt);
    case 9:
        return -rt3*(2.0*Ahydro[21] + 4.0*Ahydro[22]*rt2);
    case 10:
        return -rt3*(2.0*Ahydro[23] + 3.0*Ahydro[24]*rt);
    case 11:
        return -rt3*(2.0*Ahydro[25] + 4.0*Ahydro[26]*rt2);
    case 12:
        return -rt3*(2.0*Ahydro[27] + 3.0*Ahydro[28]*rt);
    case 13:
        return -rt3*(2.0*Ahydro[29] + 3.0*Ahydro[30]*rt + 4.0*Ahydro[31]*rt2);
    default:
        return 0.0;
    }
}

double hydrogen::W(int n, double egrho)
{
    return (n == 0 ? (1.0 - egrho)/(2.0*Gamma) :
            (n*W(n-1, egrho) - 0.5*pow(Rho,2*n)*egrho)/Gamma);
}

double hydrogen::H(int i, double egrho)
{
    return (i < 8 ? pow(Rho,i+2) : pow(Rho,2*i-13)*egrho);
}

double hydrogen::I(int i, double egrho)
{
    return (i < 8 ? pow(Rho,i+1)/double(i+1) : W(i-8, egrho));
}

double hydrogen::icv(int i, double x, double xlg)
{
    return (i == 0 ? x - 1 : x*pow(xlg,i) - i*icv(i-1,x,xlg));
}

double hydrogen::up()
{
    double rt = 1.0/T;
    double rt2 = rt*rt;
    double rt3 = rt*rt2;
    double egrho = exp(-Gamma*Rho*Rho);
    double sum = u0;
    for (int i=0; i<14; i++) {
        sum += (C(i, rt, rt2) - T*Cprime(i, rt, rt2, rt3))*I(i, egrho);
    }

    //   add \int c_{v,0} term
    sum += Ghydro[0] * (std::min(T, T1) - To);
    if (T > T1) {
        double x = std::min(T, T2) / T1;
        for (int i = 0; i < 12; i++) {
            sum += Ghydro[i] * T1 * icv(i, x, log(x));
        }
    }
    if (T > T2) {
        double x = T/T2;
        for (int i = 0; i < 5; i++) {
            sum += Ghydro[i+12] * T2 * icv(i, x, log(x));
        }
    }
    return sum + m_energy_offset;
}

double hydrogen::sp()
{
    double rt = 1.0/T;
    double rt2 = rt*rt;
    double rt3 = rt*rt2;
    double egrho = exp(-Gamma*Rho*Rho);
    double sum = s0 - R*log(Rho);
    for (int i=0; i<14; i++) {
        sum -= Cprime(i, rt, rt2, rt3)*I(i, egrho);
    }

    //   add \int c_{v,0}/T term
    sum += Ghydro[0] * log(std::min(T, T1)/ To);
    if (T > T1) {
        double xlg = log(std::min(T, T2)/T1);
        for (int i = 0; i < 12; i++) {
            sum += Ghydro[i] / (i + 1) * pow(xlg, i+1);
        }
    }
    if (T > T2) {
        double xlg = log(T/T2);
        for (int i = 0; i < 5; i++) {
            sum += Ghydro[i+12] / (i + 1) * pow(xlg, i+1);
        }
    }
    return sum + m_entropy_offset;
}

double hydrogen::Pp()
{
    double rt = 1.0/T;
    double rt2 = rt*rt;
    double egrho = exp(-Gamma*Rho*Rho);

    double P = Rho*R*T;
    for (int i=0; i<14; i++) {
        P += C(i, rt, rt2)*H(i, egrho);
    }
    return P;
}

double hydrogen::ldens()
{
    if ((T < Tmn) || (T > Tc)) {
        throw CanteraError("hydrogen::ldens",
                           "Temperature out of range. T = {}", T);
    }
    double x=1-T/Tc;
    double sum;
    int i;
    for (i=1, sum=0; i<=6; i++) {
        sum+=Dhydro[i]*pow(x, 1+double(i-1)/3.0);
    }
    return sum+Roc+Dhydro[0]*pow(x,alpha1);
}

double hydrogen::Psat()
{
    double x = (1.0 - Tt/T)/(1.0 - Tt/Tc);
    double result;
    if ((T < Tmn) || (T > Tc)) {
        throw CanteraError("hydrogen::Psat",
                           "Temperature out of range. T = {}", T);
    }
    result = Fhydro[0]*x + Fhydro[1]*x*x + Fhydro[2]*x*x*x +
             Fhydro[3]*x*pow(1-x, alpha);
    return exp(result)*Pt;
}

double hydrogen::Tcrit()
{
    return Tc;
}
double hydrogen::Pcrit()
{
    return Pc;
}
double hydrogen::Vcrit()
{
    return 1.0/Roc;
}
double hydrogen::Tmin()
{
    return Tmn;
}
double hydrogen::Tmax()
{
    return Tmx;
}
double hydrogen::MolWt()
{
    return M;
}

}
