//! @file Nitrogen.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "Nitrogen.h"
#include "cantera/base/stringUtils.h"

using namespace Cantera;

namespace tpx
{
static const double M = 28.01348,
                    Tmn = 63.15,
                    Tmx = 2000.0,
                    Tc = 126.200,
                    Pc = 3.4e6,
                    Roc = 314.03,
                    R = 2.96790515164171e2,
                    Gamma = 7.13602531283233e-6,
                    alpha = 1.95,
                    beta = 3353.40610,
                    u0 = 150877.551,
                    s0 = 214.9352518;

static const double Ann[] = {
    1.75889959256970e-1,  1.38197604384933e1,  -3.14918412133921e2,
    4.40300150239380e3,  -5.45358971644916e5, 4.84413320182919e-4,
    -5.18964416491365e-2,  6.57265859197103e-4 ,8.51299771713314e4 ,
    1.33459405162578e-8,  3.83381319826746e-4, -8.35421151028455e-2,
    2.84874912286101e-7,
    -2.38296116270360e-7, -1.48321912935764e-4, 5.62605853190540e-10,
    -2.98201050924595e-13, 9.85319087685241e-11, -1.92002176056468e-14,
    -7.82250103373122e4,  -5.51801778744598e5, -5.72781957607352e-1,
    3.25760529488327e2,  -1.34659309828737e-6, -1.92036423064911e-5,
    -3.94564337674524e-12,-2.44388245328965e-9, -1.50970602460077e-18,
    1.25854885346038e-16,-8.34271144923969e-24, -1.17299202018417e-22,
    9.06544823455730e-22
};

static const double Fnn[]= {
    8.3944094440e3, -1.8785191705e3, -7.2822291650,
    1.0228509660e-2, 5.5560638250e-4,
    -5.9445446620e-6, 2.7154339320e-8,
    -4.8795359040e-11, 5.0953608240e2
};

static const double Dnn[] = {
    3.1402991e2, 4.4111015e2, 9.4622994e2 ,
    -2.9067111e3, 4.4785979e3, -2.2746914e3
};

static const double Gnn[] = {
    -2.18203473713518e5, 1.01573580096247e4, -1.65504721657240e2,
    7.43175999190430e2, -5.14605623546025e-3,
    5.18347156760489e-6, -1.05922170493616e-9, 2.98389393363817e2
};

double nitrogen::C(int i, double rt, double rt2)
{
    switch (i) {
    case 0:
        return Ann[0] * T + Ann[1] * sqrt(T)
               + Ann[2] + (Ann[3] + Ann[4] * rt) * rt;
    case 1:
        return Ann[5] * T + Ann[6] + rt * (Ann[7] + Ann[8] * rt);
    case 2:
        return Ann[9] * T + Ann[10] + Ann[11] * rt;
    case 3:
        return Ann[12];
    case 4:
        return rt*(Ann[13] + Ann[14]*rt);
    case 5:
        return Ann[15]*rt;
    case 6:
        return rt*(Ann[16] + Ann[17]*rt);
    case 7:
        return Ann[18]*rt2;
    case 8:
        return rt2*(Ann[19] + Ann[20]*rt);
    case 9:
        return rt2*(Ann[21] + Ann[22]*rt2);
    case 10:
        return rt2*(Ann[23] + Ann[24]*rt);
    case 11:
        return rt2*(Ann[25] + Ann[26]*rt2);
    case 12:
        return rt2*(Ann[27] + Ann[28]*rt);
    case 13:
        return rt2*(Ann[29] + Ann[30]*rt + Ann[31]*rt2);
    default:
        return 0.0;
    }
}

double nitrogen::Cprime(int i, double rt, double rt2, double rt3)
{
    switch (i) {
    case 0:
        return Ann[0] + 0.5*Ann[1]/sqrt(T) - (Ann[3] + 2.0*Ann[4]*rt)*rt2;
    case 1:
        return Ann[5] - rt2*(Ann[7] + 2.0*Ann[8]*rt);
    case 2:
        return Ann[9] - Ann[11]*rt2;
    case 3:
        return 0.0;
    case 4:
        return -rt2*(Ann[13] + 2.0*Ann[14]*rt);
    case 5:
        return -Ann[15]*rt2;
    case 6:
        return -rt2*(Ann[16] + 2.0*Ann[17]*rt);
    case 7:
        return -2.0*Ann[18]*rt3;
    case 8:
        return -rt3*(2.0*Ann[19] + 3.0*Ann[20]*rt);
    case 9:
        return -rt3*(2.0*Ann[21] + 4.0*Ann[22]*rt2);
    case 10:
        return -rt3*(2.0*Ann[23] + 3.0*Ann[24]*rt);
    case 11:
        return -rt3*(2.0*Ann[25] + 4.0*Ann[26]*rt2);
    case 12:
        return -rt3*(2.0*Ann[27] + 3.0*Ann[28]*rt);
    case 13:
        return -rt3*(2.0*Ann[29] + 3.0*Ann[30]*rt + 4.0*Ann[31]*rt2);
    default:
        return 0.0;
    }
}

double nitrogen::W(int n, double egrho)
{
    return (n == 0 ? (1.0 - egrho)/(2.0*Gamma) :
            (n*W(n-1, egrho) - 0.5*pow(Rho,2*n)*egrho)/Gamma);
}

double nitrogen::H(int i, double egrho)
{
    return (i < 8 ? pow(Rho,i+2) : pow(Rho,2*i-13)*egrho);
}

double nitrogen::I(int i, double egrho)
{
    return (i < 8 ? pow(Rho,i+1)/double(i+1) : W(i-8, egrho));
}

double nitrogen::up()
{
    double rt = 1.0/T;
    double rt2 = rt*rt;
    double rt3 = rt*rt2;
    double egrho = exp(-Gamma*Rho*Rho);

    double sum = 0.0;
    for (int i=0; i<14; i++) {
        sum += (C(i,rt,rt2) - T*Cprime(i,rt,rt2,rt3))*I(i,egrho);
    }
    sum += (((0.25*Gnn[6]*T + Gnn[5]/3.0)*T
             + 0.5*Gnn[4])*T + Gnn[3])*T + Gnn[2]*log(T)
           - (Gnn[1] + 0.5*Gnn[0]*rt)*rt
           + Gnn[7]*beta/(exp(beta*rt) - 1.0) + u0
           + m_energy_offset;
    return sum;
}

double nitrogen::sp()
{
    double rt = 1.0/T;
    double rt2 = rt*rt;
    double rt3 = rt*rt2;
    double egrho = exp(-Gamma*Rho*Rho);

    double sum = 0.0;
    sum = s0 + m_entropy_offset - R*log(Rho);
    for (int i=0; i<14; i++) {
        sum -= Cprime(i,rt,rt2,rt3)*I(i,egrho);
    }
    sum += (((Gnn[6]/3.0)*T + 0.5*Gnn[5])*T + Gnn[4])*T + Gnn[3]*log(T)
           -((Gnn[0]*rt/3.0 + 0.5*Gnn[1])*rt + Gnn[2])*rt
           + Gnn[7]*(beta*rt + beta*rt/(exp(beta*rt) - 1.0)
                     - log(exp(beta*rt) - 1.0));
    return sum;
}

double nitrogen::Pp()
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

double nitrogen::Psat()
{
    double lnp;
    int i;
    if ((T < Tmn) || (T > Tc)) {
        throw CanteraError("nitrogen::Psat",
                           "Temperature out of range. T = {}", T);
    }
    for (i=0, lnp=0; i<=7; i++) {
        if (i==3) {
            lnp+=Fnn[i]*pow(Tc-T, alpha);
        } else {
            lnp+=Fnn[i]*pow(T,i-1);
        }
    }
    lnp+=Fnn[8]*log(T);
    return exp(lnp);
}

double nitrogen::ldens()
{
    double xx=1-T/Tc, sum=0;
    if ((T < Tmn) || (T > Tc)) {
        throw CanteraError("nitrogen::ldens",
                           "Temperature out of range. T = {}", T);
    }
    for (int i=0; i<=5; i++) {
        sum+=Dnn[i]*pow(xx,double(i)/3.0);
    }
    return sum;
}

double nitrogen::Tcrit()
{
    return Tc;
}
double nitrogen::Pcrit()
{
    return Pc;
}
double nitrogen::Vcrit()
{
    return 1.0/Roc;
}
double nitrogen::Tmin()
{
    return Tmn;
}
double nitrogen::Tmax()
{
    return Tmx;
}
double nitrogen::MolWt()
{
    return M;
}
}
