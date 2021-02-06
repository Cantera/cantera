//! @file Methane.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "Methane.h"
#include "cantera/base/stringUtils.h"

using namespace std;
using namespace Cantera;

namespace tpx
{
static const double
M = 16.04996,
Tmn = 90.68,
Tmx = 1700.0,
Tc = 190.555,
Pc = 4.5988e6,
Roc= 160.43,
Tt = 90.68,
Pt=11743.5675,
R = 5.18253475866e2,
Gamma=3.72992471469e-5,
alpha = 1.5, //Used with Psat
alpha1 = .36, //used with ldens;
Rot=451.562,
beta=2009.152,
u0 = 357696.0858,
s0 = -1918.035071;

static const double Ameth[] = {
    -7.25929210183, 4.13766054566e2, -6.32167316855e3,
    3.34015577724e5, -1.68253379982e7, 1.87884851902e-2, -1.18673201223e1,
    2.09062618015e3, -4.07532656958e5, -5.73917603241e-5,4.37711441593e-2,
    -4.38766500673, 1.13524630779e-5, -5.07028240949e-5, 2.28002199522e-2,
    9.25611329590e-9, 1.33865662546e-10, -1.65439044196e-7, 1.81030980110e-10,
    5.45753645958e5, -3.63192281933e7, 4.81463473761, 1.56633022620e5,
    7.89977010972e-5, 1.39993881210e-2, -1.70656092212e-11, -4.55256623445e-5,
    -2.29314170748e-14,8.31548197665e-12, 6.84673626259e-20,
    -4.70845544152e-17, 5.21465091383e-16
};

static const double Dmeth[]=
{  -1.78860165e-1, 4.83847500e-2, -1.84898700e-2  };

static const double Fmeth[]=
{  4.77748580, 1.76065363, -5.67888940e-1, 1.32786231  };

static const double Gmeth[]=
{  1.34740610e3, 1.35512060e2, -2.93910458e1, 2.12774600, 2.44656600e3  };

double methane::C(int i, double rt, double rt2)
{
    switch (i) {
    case 0:
        return Ameth[0] * T + Ameth[1] * sqrt(T) + Ameth[2] + (Ameth[3] + Ameth[4] * rt) * rt;
    case 1:
        return Ameth[5] * T + Ameth[6] + rt * (Ameth[7] + Ameth[8] * rt);
    case 2:
        return Ameth[9] * T + Ameth[10] + Ameth[11] * rt;
    case 3:
        return Ameth[12];
    case 4:
        return rt*(Ameth[13] + Ameth[14]*rt);
    case 5:
        return Ameth[15]*rt;
    case 6:
        return rt*(Ameth[16] + Ameth[17]*rt);
    case 7:
        return Ameth[18]*rt2;
    case 8:
        return rt2*(Ameth[19] + Ameth[20]*rt);
    case 9:
        return rt2*(Ameth[21] + Ameth[22]*rt2);
    case 10:
        return rt2*(Ameth[23] + Ameth[24]*rt);
    case 11:
        return rt2*(Ameth[25] + Ameth[26]*rt2);
    case 12:
        return rt2*(Ameth[27] + Ameth[28]*rt);
    case 13:
        return rt2*(Ameth[29] + Ameth[30]*rt + Ameth[31]*rt2);
    default:
        return 0.0;
    }
}

double methane::Cprime(int i, double rt, double rt2, double rt3)
{
    switch (i) {
    case 0:
        return Ameth[0] + 0.5*Ameth[1]/sqrt(T) - (Ameth[3] + 2.0*Ameth[4]*rt)*rt2;
    case 1:
        return Ameth[5] - rt2*(Ameth[7] + 2.0*Ameth[8]*rt);
    case 2:
        return Ameth[9] - Ameth[11]*rt2;
    case 3:
        return 0.0;
    case 4:
        return -rt2*(Ameth[13] + 2.0*Ameth[14]*rt);
    case 5:
        return -Ameth[15]*rt2;
    case 6:
        return -rt2*(Ameth[16] + 2.0*Ameth[17]*rt);
    case 7:
        return -2.0*Ameth[18]*rt3;
    case 8:
        return -rt3*(2.0*Ameth[19] + 3.0*Ameth[20]*rt);
    case 9:
        return -rt3*(2.0*Ameth[21] + 4.0*Ameth[22]*rt2);
    case 10:
        return -rt3*(2.0*Ameth[23] + 3.0*Ameth[24]*rt);
    case 11:
        return -rt3*(2.0*Ameth[25] + 4.0*Ameth[26]*rt2);
    case 12:
        return -rt3*(2.0*Ameth[27] + 3.0*Ameth[28]*rt);
    case 13:
        return -rt3*(2.0*Ameth[29] + 3.0*Ameth[30]*rt + 4.0*Ameth[31]*rt2);
    default:
        return 0.0;
    }
}

double methane::W(int n, double egrho)
{
    return (n == 0 ? (1.0 - egrho)/(2.0*Gamma) :
            (n*W(n-1, egrho) - 0.5*pow(Rho,2*n)*egrho)/Gamma);
}

double methane::H(int i, double egrho)
{
    return (i < 8 ? pow(Rho,i+2) : pow(Rho,2*i-13)*egrho);
}

double methane::I(int i, double egrho)
{
    return (i < 8 ? pow(Rho,i+1)/double(i+1) : W(i-8, egrho));
}

double methane::up()
{
    double rt = 1.0/T;
    double rt2 = rt*rt;
    double rt3 = rt*rt2;
    double egrho = exp(-Gamma*Rho*Rho);
    double t3 = pow(T,1.0/3.0);

    double sum = 0.0;
    for (int i=0; i<14; i++) {
        sum += (C(i, rt, rt2) - T*Cprime(i, rt, rt2, rt3))*I(i, egrho);
    }
    sum += T*(Gmeth[0] + 0.75*Gmeth[1]*t3 + 0.6*Gmeth[2]*t3*t3 + 0.5*Gmeth[3]*T)
           + Gmeth[4]*beta/(exp(beta*rt) - 1.0) + u0;
    return sum + m_energy_offset;
}

double methane::sp()
{
    double rt = 1.0/T;
    double rt2 = rt*rt;
    double rt3 = rt*rt2;
    double egrho = exp(-Gamma*Rho*Rho);
    double t3 = pow(T,1.0/3.0);
    double sum = 0.0;
    sum = s0 - R*log(Rho);
    for (int i=0; i<14; i++) {
        sum -= Cprime(i, rt, rt2, rt3)*I(i, egrho);
    }
    sum += Gmeth[0]*log(T) + 3.0*Gmeth[1]*t3 + 1.5*Gmeth[2]*t3*t3 + Gmeth[3]*T
           + Gmeth[4]*(beta*rt + beta*rt/(exp(beta*rt) - 1.0)
                       - log(exp(beta*rt) - 1.0));
    return sum + m_entropy_offset;
}

double methane::Pp()
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

double methane::Psat()
{
    double x = (1.0 - Tt/T)/(1.0 - Tt/Tc);
    double result;
    if ((T < Tmn) || (T > Tc)) {
        throw CanteraError("methane::Psat",
                           "Temperature out of range. T = {}", T);
    }
    result = Fmeth[0]*x + Fmeth[1]*x*x + Fmeth[2]*x*x*x +
             Fmeth[3]*x*pow(1-x, alpha);
    return exp(result)*Pt;
}

double methane::ldens()
{
    double result;
    double sum;
    double w;
    if ((T < Tmn) || (T > Tc)) {
        throw CanteraError("methane::ldens",
                           "Temperature out of range. T = {}", T);
    }
    w = (Tc - T)/(Tc - Tt);
    sum = Dmeth[0]*(1.0 - pow(w, 2.0/3.0)) + Dmeth[1]*(1.0 - pow(w, 4.0/3.0))
          + Dmeth[2]*(1.0 - pow(w, 2));
    result = pow(w,alpha1)*exp(sum);
    result *= (Rot-Roc);
    result += Roc;
    return result;
}

double methane::Tcrit()
{
    return Tc;
}
double methane::Pcrit()
{
    return Pc;
}
double methane::Vcrit()
{
    return 1.0/Roc;
}
double methane::Tmin()
{
    return Tmn;
}
double methane::Tmax()
{
    return Tmx;
}
double methane::MolWt()
{
    return M;
}

}
