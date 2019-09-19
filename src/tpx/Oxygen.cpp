//! @file Oxygen.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "Oxygen.h"
#include "cantera/base/stringUtils.h"

using namespace Cantera;

namespace tpx
{
static const double
M = 31.9994,
Tmn = 54.34,
Tmx = 2000.0,
Tc = 154.581,
Pc = 5.0429e6,
Roc = 436.15,
R = 2.59820853437877e2,
Gamma = 5.46895508389297e-6,
alpha = 1.91576,
beta = 2239.18105,
u0 = 198884.2435,
s0 = 668.542976;

static const double Aoxy[] = {
    -4.26396872798684e-1, 3.48334938784107e1, -5.77516910418738e2,
    2.40961751553325e4, -1.23332307855543e6, 3.73585286319658e-4,
    -1.70178244046465e-1 ,-3.33226903068473e-4, 8.61334799901291e3,
    -6.80394661057309e-7, 7.09583347162704e-4, -5.73905688255053e-2,
    -1.92123080811409e-7, 3.11764722329504e-8, -8.09463854745591e-6,
    -2.22562296356501e-11, 9.18401045361994e-15, 5.75758417511114e-12,
    -2.10752269644774e-15, 3.62884761272184e3, -1.23317754317110e6,
    -5.03800414800672e-2, 3.30686173177055e2, -5.26259633964252e-8 ,
    5.53075442383100e-6, -2.71042853363688e-13, -1.65732450675251e-9 ,
    -5.82711196409204e-20, 4.42953322148281e-17 ,-2.95529679136244e-25,
    -1.92361786708846e-23, 9.43758410350413e-23
};

static const double Foxy[] = {
    -5.581932039e2, -1.0966262185e2, -8.3456211630e-2,
    2.6603644330e-3, 1.6875023830e-5, -2.1262477120e-7,
    9.5741096780e-10, -1.6617640450e-12, 2.7545605710e1
};

static const double Doxy[] =
{  4.3615175e2, 7.5897189e2, -4.2576866e2, 2.3487106e3, -3.0474660e3, 1.4850169e3  };

static const double Goxy[] = {
    -1.29442711174062e6, 5.98231747005341e4, -8.97850772730944e2,
    6.55236176900400e2, -1.13131252131570e-2,
    3.4981070244228e-6, 4.21065222886885e-9, 2.67997030050139e2
};

double oxygen::C(int i, double rt, double rt2)
{
    switch (i) {
    case 0:
        return Aoxy[0] * T + Aoxy[1] * sqrt(T) + Aoxy[2] + (Aoxy[3] + Aoxy[4] * rt) * rt;
    case 1:
        return Aoxy[5] * T + Aoxy[6] + rt * (Aoxy[7] + Aoxy[8] * rt);
    case 2:
        return Aoxy[9] * T + Aoxy[10] + Aoxy[11] * rt;
    case 3:
        return Aoxy[12];
    case 4:
        return rt*(Aoxy[13] + Aoxy[14]*rt);
    case 5:
        return Aoxy[15]*rt;
    case 6:
        return rt*(Aoxy[16] + Aoxy[17]*rt);
    case 7:
        return Aoxy[18]*rt2;
    case 8:
        return rt2*(Aoxy[19] + Aoxy[20]*rt);
    case 9:
        return rt2*(Aoxy[21] + Aoxy[22]*rt2);
    case 10:
        return rt2*(Aoxy[23] + Aoxy[24]*rt);
    case 11:
        return rt2*(Aoxy[25] + Aoxy[26]*rt2);
    case 12:
        return rt2*(Aoxy[27] + Aoxy[28]*rt);
    case 13:
        return rt2*(Aoxy[29] + Aoxy[30]*rt + Aoxy[31]*rt2);
    default:
        return 0.0;
    }
}

double oxygen::Cprime(int i, double rt, double rt2, double rt3)
{
    switch (i) {
    case 0:
        return Aoxy[0] + 0.5*Aoxy[1]/sqrt(T) - (Aoxy[3] + 2.0*Aoxy[4]*rt)*rt2;
    case 1:
        return Aoxy[5] - rt2*(Aoxy[7] + 2.0*Aoxy[8]*rt);
    case 2:
        return Aoxy[9] - Aoxy[11]*rt2;
    case 3:
        return 0.0;
    case 4:
        return -rt2*(Aoxy[13] + 2.0*Aoxy[14]*rt);
    case 5:
        return -Aoxy[15]*rt2;
    case 6:
        return -rt2*(Aoxy[16] + 2.0*Aoxy[17]*rt);
    case 7:
        return -2.0*Aoxy[18]*rt3;
    case 8:
        return -rt3*(2.0*Aoxy[19] + 3.0*Aoxy[20]*rt);
    case 9:
        return -rt3*(2.0*Aoxy[21] + 4.0*Aoxy[22]*rt2);
    case 10:
        return -rt3*(2.0*Aoxy[23] + 3.0*Aoxy[24]*rt);
    case 11:
        return -rt3*(2.0*Aoxy[25] + 4.0*Aoxy[26]*rt2);
    case 12:
        return -rt3*(2.0*Aoxy[27] + 3.0*Aoxy[28]*rt);
    case 13:
        return -rt3*(2.0*Aoxy[29] + 3.0*Aoxy[30]*rt + 4.0*Aoxy[31]*rt2);
    default:
        return 0.0;
    }
}

double oxygen::W(int n, double egrho)
{
    return (n == 0 ? (1.0 - egrho)/(2.0*Gamma) :
            (n*W(n-1, egrho) - 0.5*pow(Rho,2*n)*egrho)/Gamma);
}

double oxygen::H(int i, double egrho)
{
    return (i < 8 ? pow(Rho,i+2) : pow(Rho,2*i-13)*egrho);
}

double oxygen::I(int i, double egrho)
{
    return (i < 8 ? pow(Rho,i+1)/double(i+1) : W(i-8, egrho));
}

double oxygen::up()
{
    double rt = 1.0/T;
    double rt2 = rt*rt;
    double rt3 = rt*rt2;
    double egrho = exp(-Gamma*Rho*Rho);

    double sum = 0.0;
    for (int i=0; i<14; i++) {
        sum += (C(i,rt,rt2) - T*Cprime(i,rt,rt2,rt3))*I(i,egrho);
    }
    sum += (((0.25*Goxy[6]*T + Goxy[5]/3.0)*T + 0.5*Goxy[4])*T + Goxy[3])*T + Goxy[2]*log(T)
           - (Goxy[1] + 0.5*Goxy[0]*rt)*rt + Goxy[7]*beta/(exp(beta*rt) - 1.0) + u0;
    return sum + m_energy_offset;
}

double oxygen::sp()
{
    double rt = 1.0/T;
    double rt2 = rt*rt;
    double rt3 = rt*rt2;
    double egrho = exp(-Gamma*Rho*Rho);

    double sum = 0.0;
    sum = s0 - R*log(Rho);
    for (int i=0; i<14; i++) {
        sum -= Cprime(i,rt,rt2,rt3)*I(i,egrho);
    }
    sum += (((Goxy[6]/3.0)*T + 0.5*Goxy[5])*T + Goxy[4])*T + Goxy[3]*log(T)
           -((Goxy[0]*rt/3.0 + 0.5*Goxy[1])*rt + Goxy[2])*rt
           + Goxy[7]*(beta*rt + beta*rt/(exp(beta*rt) - 1.0)
                      - log(exp(beta*rt) - 1.0));
    return sum + m_entropy_offset;
}

double oxygen::Pp()
{
    double rt = 1.0/T;
    double rt2 = rt*rt;
    double egrho = exp(-Gamma*Rho*Rho);

    double P = Rho*R*T;
    for (int i=0; i<14; i++) {
        P += C(i,rt,rt2)*H(i,egrho);
    }
    return P;
}

double oxygen::Psat()
{
    double lnp;
    int i;
    if ((T < Tmn) || (T > Tc)) {
        throw CanteraError("oxygen::Psat",
                           "Temperature out of range. T = {}", T);
    }
    for (i=0, lnp=0; i<=7; i++) {
        if (i==3) {
            lnp+=Foxy[i]*pow(Tc-T, alpha);
        } else {
            lnp+=Foxy[i]*pow(T,i-1);
        }
    }
    lnp+=Foxy[8]*log(T);
    return exp(lnp);
}

double oxygen::ldens()
{
    double xx=1-T/Tc, sum=0;
    if ((T < Tmn) || (T > Tc)) {
        throw CanteraError("oxygen::ldens",
                           "Temperature out of range. T = {}", T);
    }
    for (int i=0; i<=5; i++) {
        sum+=Doxy[i]*pow(xx,double(i)/3.0);
    }
    return sum;
}

double oxygen::Tcrit()
{
    return Tc;
}
double oxygen::Pcrit()
{
    return Pc;
}
double oxygen::Vcrit()
{
    return 1.0/Roc;
}
double oxygen::Tmin()
{
    return Tmn;
}
double oxygen::Tmax()
{
    return Tmx;
}
double oxygen::MolWt()
{
    return M;
}

}
