//! @file HFC134a.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "HFC134a.h"
#include "cantera/base/stringUtils.h"

using namespace Cantera;

namespace tpx
{

const double
M = 102.032,
Tmn = 170.0,
Tmx = 455.0,
Tc = 374.18,
Pc = 4056290.0,
Roc = 508.0,
R = 81.48885644;

const double a134[] = {
    0.5586817e-1,
    0.4982230,
    0.2458698e-1,
    0.8570145e-3,
    0.4788584e-3,
    -0.1800808e1,
    0.2671641,
    -0.4781652e-1,
    0.1423987e-1,
    0.3324062,
    -0.7485907e-2,
    0.1017263e-3,
    -0.5184567,
    -0.8692288e-1,
    0.2057144,
    -0.5000457e-2,
    0.4603262e-3,
    -0.3497836e-2,
    0.6995038e-2,
    -0.1452184e-1,
    -0.1285458e-3
};

const double t134[] = {
    -0.5, 0.0, 0.0, 0.0, 1.5, 1.5, 2.0, 2.0, 1.0, 3.0, 5.0,
    1.0, 5.0, 5.0, 6.0, 10.0, 10.0, 10.0, 18.0, 22.0, 50.0
};

const int d134[] = {
    2, 1, 3, 6, 6, 1, 1, 2, 5, 2, 2,
    4, 1, 4, 1, 2, 4, 1, 5, 3, 10
};

const double b134[] = {
    -1.019535,
    9.047135,
    -1.629789,
    -9.723916,
    -3.927170
};

double HFC134a::fp()
{
    double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0,
           sum4 = 0.0, sum5 = 0.0;
    double tau = Tc/T;
    double delta = Rho/Roc;

    double phi0 = b134[0] + b134[1]*tau + b134[2]*log(tau)
                  + log(delta) + b134[3]/sqrt(tau) + b134[4]*pow(tau,-0.75);
    int i;
    for (i = 0; i<8; i++) {
        sum1 += a134[i]*pow(tau,t134[i])*pow(delta,d134[i]);
    }
    for (i = 8; i<11; i++) {
        sum2 += a134[i]*pow(tau,t134[i])*pow(delta,d134[i]);
    }
    for (i = 11; i<17; i++) {
        sum3 += a134[i]*pow(tau,t134[i])*pow(delta,d134[i]);
    }
    for (i = 17; i<20; i++) {
        sum4 += a134[i]*pow(tau,t134[i])*pow(delta,d134[i]);
    }
    sum5 = a134[20]*pow(tau,t134[20])*pow(delta,d134[20]);
    double phir = sum1 + exp(-delta)*sum2 + exp(-delta*delta)*sum3
                  + exp(-delta*delta*delta)*sum4
                  + exp(-delta*delta*delta*delta)*sum5;
    return R*T*(phir + phi0);
}

double HFC134a::up()
{
    double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0,
           sum4 = 0.0, sum5 = 0.0;
    double tau = Tc/T;
    double delta = Rho/Roc;

    double phi0t = b134[1]*tau + b134[2]
                   - 0.5*b134[3]*pow(tau,-0.5) - 0.75*b134[4]*pow(tau,-0.75);
    int i;
    for (i = 0; i<8; i++) {
        sum1 += a134[i]*t134[i]*pow(tau,t134[i])*pow(delta,d134[i]);
    }
    for (i = 8; i<11; i++) {
        sum2 += a134[i]*t134[i]*pow(tau,t134[i])*pow(delta,d134[i]);
    }
    for (i = 11; i<17; i++) {
        sum3 += a134[i]*t134[i]*pow(tau,t134[i])*pow(delta,d134[i]);
    }
    for (i = 17; i<20; i++) {
        sum4 += a134[i]*t134[i]*pow(tau,t134[i])*pow(delta,d134[i]);
    }
    sum5 = a134[20]*t134[20]*pow(tau,t134[20])*pow(delta,d134[20]);
    double phirt = sum1 + exp(-delta)*sum2 + exp(-delta*delta)*sum3
                   + exp(-delta*delta*delta)*sum4
                   + exp(-delta*delta*delta*delta)*sum5;
    return R*T*(phirt + phi0t) + m_energy_offset;
}

double HFC134a::Pp()
{
    double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0,
           sum4 = 0.0, sum5 = 0.0;
    double tau = Tc/T;
    double delta = Rho/Roc;
    double phi0d = 1.0/delta;

    int i;
    for (i = 0; i<8; i++) {
        sum1 += a134[i]*pow(tau,t134[i])*d134[i]*pow(delta,d134[i]-1);
    }
    for (i = 8; i<11; i++) {
        sum2 += a134[i]*pow(tau,t134[i])*(d134[i] - delta)*pow(delta,d134[i]-1);
    }
    sum2 *= exp(-delta);
    double dk = delta*delta;
    for (i = 11; i<17; i++) {
        sum3 += a134[i]*pow(tau,t134[i])*(d134[i] - 2.0*dk)*pow(delta,d134[i]-1);
    }
    sum3 *= exp(-dk);
    dk *= delta;
    for (i = 17; i<20; i++) {
        sum4 += a134[i]*pow(tau,t134[i])*(d134[i] - 3.0*dk)*pow(delta,d134[i]-1);
    }
    sum4 *= exp(-dk);
    dk *= delta;
    sum5 = a134[20]*pow(tau,t134[20])*(d134[20] - 4.0*dk)*pow(delta,d134[20]-1);
    sum5 *= exp(-dk);
    double phird = sum1 + sum2 + sum3 + sum4 + sum5;
    return R*T*delta*delta*Roc*(phird + phi0d);
}

double HFC134a::Psat()
{
    if ((T < Tmn) || (T > Tc)) {
        throw CanteraError("HFC134a::Psat",
                           "Temperature out of range. T = {}", T);
    }
    double x1 = T/Tc;
    double x2 = 1.0 - x1;
    double f = -7.686556*x2 + 2.311791*pow(x2,1.5)
               - 2.039554*x2*x2 - 3.583758*pow(x2,4);
    return Pc*exp(f/x1);
}

double HFC134a::ldens()
{
    if ((T < Tmn) || (T > Tc)) {
        throw CanteraError("HFC134a::ldens",
                           "Temperature out of range. T = {}", T);
    }
    double x1 = T/Tc;
    double x2 = 1.0 - x1;
    return 518.2 + 884.13*pow(x2,1.0/3.0) + 485.84*pow(x2,2.0/3.0)
           + 193.29*pow(x2,10.0/3.0);
}

double HFC134a::Tcrit()
{
    return 374.21;
}
double HFC134a::Pcrit()
{
    return 4059280.0;
}
double HFC134a::Vcrit()
{
    return 1.0/511.95;
}
double HFC134a::Tmin()
{
    return Tmn;
}
double HFC134a::Tmax()
{
    return Tmx;
}
double HFC134a::MolWt()
{
    return M;
}

}
