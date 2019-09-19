//! @file Water.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "Water.h"
#include "cantera/base/stringUtils.h"

using namespace Cantera;

namespace tpx
{
static const double Tmn=273.16;
static const double Tmx=1600.0;
static const double M=18.016;
static const double Tc=647.286;
static const double Pc=22.089e6;
static const double Roc=317.0;
static const double To=273.16;
static const double R=461.51;
static const double E=4.8E-3;
static const double Ta=1000.0;
static const double tauc=1.544912;
static const double Tp=338.15;
static const double aww=0.01;
static const double Roa1=634.0;
static const double Roaj=1000.0;
static const double u0=2375470.875;
static const double s0=6697.356635;

static const double A[10][7]= {{
        2.9492937E-2,-5.1985860E-3,
        6.8335354E-3,-1.5641040E-4,
        -6.3972405E-3, -3.9661401E-3, -6.9048554E-4
    },
    {
        -1.3213917E-4,7.7779182E-6, -2.6149751E-5,-7.2546108E-7,
        2.6409282E-5, 1.5453061E-5,2.7407416E-6
    },
    {
        2.7464632E-7,-3.3301902E-8,6.5326396E-8,-9.2734289E-9,
        -4.7740374E-8,-2.9142470E-8,-5.1028070E-9
    },
    {
        -3.6093828E-10, -1.6254622E-11, -2.6181978E-11, 4.3125840E-12,
        5.6323130E-11, 2.9568796E-11,3.9636085E-12
    },
    {3.4218431E-13, -1.7731074E-13,0,0,0,0,0},
    {-2.4450042E-16, 1.2748742E-16,0,0,0,0,0},
    {1.5518535E-19, 1.3746153E-19,0,0,0,0,0},
    {5.9728487E-24,1.5597836E-22, 0,0,0,0,0},
    {
        -4.1030848E-1, 3.3731180E-1, -1.3746678E-1, 6.7874983E-3,
        1.3687317E-1, 7.984797E-2, 1.3041253E-2
    },
    {
        -4.1605860E-4, -2.0988866E-4,-7.3396848E-4,1.0401717E-5,
        6.4581880E-4, 3.9917570E-4, 7.1531353E-5
    }
};

static const double F[]= {-7.4192420, 2.9721E-1,-1.155286E-1,8.685635E-3,
                          1.0940980E-3, -4.39993E-3, 2.5206580E-3, -5.2186840E-4
                         };

static const double D[]= {3.6711257,-2.8512396E1,2.2265240E2,-8.8243852E2,
                          2.0002765E3,-2.6122557E3,1.8297674E3,-5.3350520E2
                         };

static const double G[]= {4.6E4,1.011249E3,8.3893E-1,-2.19989E-4,2.466619E-7,
                          -9.704700E-11
                         };

static const double taua[] = {1.544912, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5};

double water::C(int i)
{
    double tau = Ta/T;
    return (i == 0 ? R*T : R*T*(tau - tauc)*pow(tau - taua[i],i-1));
}

double water::Cprime(int i)
{
    double tau = Ta/T;
    return (i == 0 ? R : (i == 1 ? -R*tauc :
                          -R*pow(tau - taua[i],i-2)*(tauc*(tau - taua[i])
                                  + (i-1)*tau*(tau - tauc))));
}

double water::I(int j)
{
    double factor, sum, rho_aj;
    rho_aj = (j == 0 ? Roa1 : Roaj);
    sum = 0.0;
    factor = Rho - rho_aj;
    for (int i=7; i>0; i--) {
        sum += A[i][j];
        sum *= factor;
    }
    sum += A[0][j];
    sum += (exp(-E*Rho)*(A[8][j] + A[9][j]*Rho));
    return Rho*sum;
}

double water::H(int j)
{
    double factor, sum, rho_aj;
    rho_aj = (j == 0 ? Roa1 : Roaj);
    sum = 0.0;
    factor = Rho - rho_aj;
    for (int i=6; i>0; i--) {
        sum += (A[i][j] + Rho*(i+1)*A[i+1][j]);
        sum *= factor;
    }
    sum += (A[0][j] + Rho*A[1][j]);
    sum += (exp(-E*Rho)*((1.0 - Rho*E)*A[8][j]
                         + Rho*(2.0 - Rho*E)*A[9][j]));
    sum += A[7][j]*pow(factor,7);
    return Rho*Rho*sum;
}

double water::up()
{
    double sum = 0.0;
    int i;
    for (i=0; i<7; i++) {
        sum += (C(i) - T*Cprime(i))*I(i);
    }
    for (i=1; i<6; i++) {
        sum += G[i]*(pow(T,i) - pow(To,i))/double(i);
    }
    sum += G[0]*log(T/To) + u0;
    return sum + m_energy_offset;
}

double water::sp()
{
    double sum = 0.0;
    int i;
    for (i=2; i<6; i++) {
        sum += G[i]*(pow(T,i-1) - pow(To,i-1))/double(i-1);
    }
    sum += G[1]*log(T/To);
    sum -= G[0]*(1.0/T - 1.0/To);
    sum += s0 - R*log(Rho);
    for (i=0; i<7; i++) {
        sum -= Cprime(i)*I(i);
    }
    return sum + m_entropy_offset;
}

double water::Pp()
{
    double P = Rho*R*T;
    for (int i=0; i<7; i++) {
        P += C(i)*H(i);
    }
    return P;
}

double water::Psat()
{
    double log, sum=0;
    if ((T < Tmn) || (T > Tc)) {
        throw CanteraError("water::Psat",
                           "Temperature out of range. T = {}", T);
    }
    for (int i=1; i<=8; i++) {
        sum += F[i-1]*pow(aww*(T-Tp),double(i-1)); // DGG mod
    }
    log = (Tc/T-1)*sum;
    return exp(log)*Pc;
}

double water::ldens()
{
    double sum=0;
    int i;
    if ((T < Tmn) || (T >= Tc)) {
        throw CanteraError("water::ldens",
                           "Temperature out of range. T = {}", T);
    }
    for (i=0; i<8; i++) {
        sum+=D[i]*pow(1.0 - T/Tc, double(i+1)/3.0);
    }
    return Roc*(1+sum);
}

double water::Tcrit()
{
    return Tc;
}
double water::Pcrit()
{
    return Pc;
}
double water::Vcrit()
{
    return 1.0/Roc;
}
double water::Tmin()
{
    return Tmn;
}
double water::Tmax()
{
    return Tmx;
}
double water::MolWt()
{
    return M;
}

}
