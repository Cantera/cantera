/**
 *  @file WaterProps.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/WaterProps.h"
#include "cantera/base/ctml.h"
#include "cantera/thermo/PDSS_Water.h"
#include "cantera/base/stringUtils.h"

namespace Cantera
{
WaterProps::WaterProps():
    m_waterIAPWS(0),
    m_own_sub(false)
{
    // object owns its own water evaluator
    m_waterIAPWS = new WaterPropsIAPWS();
    m_own_sub = true;
}

WaterProps::WaterProps(PDSS_Water* wptr)  :
    m_waterIAPWS(0),
    m_own_sub(false)
{
    if (wptr) {
        // object in slave mode; it doesn't own its own water evaluator.
        m_waterIAPWS = wptr->getWater();
        m_own_sub = false;
    } else {
        m_waterIAPWS = new WaterPropsIAPWS();
        m_own_sub = true;
    }
}

WaterProps::WaterProps(WaterPropsIAPWS* waterIAPWS)  :
    m_waterIAPWS(0),
    m_own_sub(false)
{
    if (waterIAPWS) {
        m_waterIAPWS = waterIAPWS;
        m_own_sub = false;
    } else {
        m_waterIAPWS = new WaterPropsIAPWS();
        m_own_sub = true;
    }
}

WaterProps::~WaterProps()
{
    if (m_own_sub) {
        delete m_waterIAPWS;
    }
}

double WaterProps::density_T(double T, double P, int ifunc)
{
    static const double Tc = T - 273.15;
    static const double U1 = 288.9414;
    static const double U2 = 508929.2;
    static const double U3 = 68.12963;
    static const double U4 = -3.9863;

    double tmp1 = Tc + U1;
    double tmp4 = Tc + U4;
    double t4t4 = tmp4 * tmp4;
    double tmp3 = Tc + U3;
    double rho = 1000. * (1.0 - tmp1*t4t4/(U2 * tmp3));

    // Impose an ideal gas lower bound on rho. We need this to ensure positivity
    // of rho, even though it is grossly unrepresentative.
    double rhomin = P / (GasConstant * T);
    if (rho < rhomin) {
        rho = rhomin;
        if (ifunc == 1) {
            return - rhomin / T;
        } else if (ifunc == 3) {
            return rhomin / P;
        } else if (ifunc == 2) {
            return 2.0 * rhomin / (T * T);
        }
    }

    if (ifunc == 1) {
        double drhodT = 1000./U2 * (
                                - tmp4 * tmp4 / (tmp3)
                                - tmp1 * 2 * tmp4 / (tmp3)
                                + tmp1 * t4t4 / (tmp3*tmp3)
                            );
        return drhodT;
    } else if (ifunc == 3) {
        return 0.0;
    } else if (ifunc == 2) {
        double t3t3 = tmp3 * tmp3;
        double d2rhodT2 = 1000./U2 *
                               ((-4.0*tmp4-2.0*tmp1)/tmp3 +
                                (2.0*t4t4 + 4.0*tmp1*tmp4)/t3t3
                                - 2.0*tmp1 * t4t4/(t3t3*tmp3));
        return d2rhodT2;
    }
    return rho;
}

double WaterProps::relEpsilon(double T, double P_pascal,
                                  int ifunc)
{
    static const double U1 = 3.4279E2;
    static const double U2 = -5.0866E-3;
    static const double U3 = 9.4690E-7;
    static const double U4 = -2.0525;
    static const double U5 = 3.1159E3;
    static const double U6 = -1.8289E2;
    static const double U7 = -8.0325E3;
    static const double U8 = 4.2142E6;
    static const double U9 = 2.1417;
    double T2 = T * T;

    double eps1000 = U1 * exp(U2 * T + U3 * T2);
    double C = U4 + U5/(U6 + T);
    double B = U7 + U8/T + U9 * T;
    double Pbar = P_pascal * 1.0E-5;
    double tmpBpar = B + Pbar;
    double tmpB1000 = B + 1000.0;
    double ltmp = log(tmpBpar/tmpB1000);
    double epsRel = eps1000 + C * ltmp;

    if (ifunc == 1 || ifunc == 2) {
        double tmpC = U6 + T;
        double dCdT = - U5/(tmpC * tmpC);
        double dBdT = - U8/(T * T) + U9;
        double deps1000dT = eps1000 * (U2 + 2.0 * U3 * T);
        double dltmpdT = (dBdT/tmpBpar - dBdT/tmpB1000);
        if (ifunc == 1) {
            return deps1000dT + dCdT * ltmp + C * dltmpdT;
        }
        double T3 = T2 * T;
        double d2CdT2 = - 2.0 * dCdT / tmpC;
        double d2BdT2 = 2.0 * U8 / (T3);
        double d2ltmpdT2 = (d2BdT2*(1.0/tmpBpar - 1.0/tmpB1000) +
                                dBdT*dBdT*(1.0/(tmpB1000*tmpB1000) - 1.0/(tmpBpar*tmpBpar)));
        double d2eps1000dT2 = (deps1000dT * (U2 + 2.0 * U3 * T) + eps1000 * (2.0 * U3));

        if (ifunc == 2) {
            double d2epsReldT2 = (d2eps1000dT2 + d2CdT2 * ltmp + 2.0 * dCdT * dltmpdT
                                      + C * d2ltmpdT2);
            return d2epsReldT2;
        }
    }
    if (ifunc == 3) {
        double dltmpdP = 1.0E-5 / tmpBpar;
        return C * dltmpdP;
    }
    return epsRel;
}

double WaterProps::ADebye(double T, double P_input, int ifunc)
{
    double psat = satPressure(T);
    double P;
    if (psat > P_input) {
        P = psat;
    } else {
        P = P_input;
    }
    double epsRelWater = relEpsilon(T, P, 0);
    double epsilon = epsilon_0 * epsRelWater;
    double dw = density_IAPWS(T, P);
    double tmp = sqrt(2.0 * Avogadro * dw / 1000.);
    double tmp2 = ElectronCharge * ElectronCharge * Avogadro /
                      (epsilon * GasConstant * T);
    double tmp3 = tmp2 * sqrt(tmp2);
    double A_Debye = tmp * tmp3 / (8.0 * Pi);

    // dAdT = - 3/2 Ad/T + 1/2 Ad/dw d(dw)/dT - 3/2 Ad/eps d(eps)/dT
    // dAdT = - 3/2 Ad/T - 1/2 Ad/Vw d(Vw)/dT - 3/2 Ad/eps d(eps)/dT
    if (ifunc == 1 || ifunc == 2) {
        double dAdT = - 1.5 * A_Debye / T;

        double depsRelWaterdT = relEpsilon(T, P, 1);
        dAdT -= A_Debye * (1.5 * depsRelWaterdT / epsRelWater);

        // calculate d(lnV)/dT _constantP, i.e., the cte
        double cte = coeffThermalExp_IAPWS(T, P);
        double contrib2 = - A_Debye * (0.5 * cte);
        dAdT += contrib2;

        if (ifunc == 1) {
            return dAdT;
        }

        if (ifunc == 2) {
            // Get the second derivative of the dielectric constant wrt T
            // -> we will take each of the terms in dAdT and differentiate
            //    it again.
            double d2AdT2 = 1.5 / T * (A_Debye/T - dAdT);
            double d2epsRelWaterdT2 = relEpsilon(T, P, 2);
            d2AdT2 += 1.5 * (- dAdT * depsRelWaterdT / epsRelWater
                             - A_Debye / epsRelWater *
                             (d2epsRelWaterdT2 - depsRelWaterdT * depsRelWaterdT / epsRelWater));
            double deltaT = -0.1;
            double Tdel = T + deltaT;
            double cte_del = coeffThermalExp_IAPWS(Tdel, P);
            double dctedT = (cte_del - cte) / Tdel;
            double contrib3 = 0.5 * (-(dAdT * cte) -(A_Debye * dctedT));
            d2AdT2 += contrib3;
            return d2AdT2;
        }
    }

    // A_Debye = (1/(8 Pi)) sqrt(2 Na dw / 1000)
    //                          (e e/(epsilon R T))^3/2
    //
    // dAdP =  + 1/2 Ad/dw d(dw)/dP - 3/2 Ad/eps d(eps)/dP
    // dAdP =  - 1/2 Ad/Vw d(Vw)/dP - 3/2 Ad/eps d(eps)/dP
    // dAdP =  + 1/2 Ad * kappa  - 3/2 Ad/eps d(eps)/dP
    //
    // where kappa = - 1/Vw d(Vw)/dP_T (isothermal compressibility)
    if (ifunc == 3) {
        double dAdP = 0.0;
        double depsRelWaterdP = relEpsilon(T, P, 3);
        dAdP -= A_Debye * (1.5 * depsRelWaterdP / epsRelWater);
        double kappa = isothermalCompressibility_IAPWS(T,P);
        dAdP += A_Debye * (0.5 * kappa);
        return dAdP;
    }
    return A_Debye;
}

double WaterProps::satPressure(double T)
{
    return m_waterIAPWS->psat(T);
}

double WaterProps::density_IAPWS(double temp, double press)
{
    return m_waterIAPWS->density(temp, press, WATER_LIQUID);
}

double WaterProps::density_IAPWS() const
{
    return m_waterIAPWS->density();
}

double WaterProps::coeffThermalExp_IAPWS(double temp, double press)
{
    double dens = m_waterIAPWS->density(temp, press, WATER_LIQUID);
    if (dens < 0.0) {
        throw CanteraError("WaterProps::coeffThermalExp_IAPWS",
            "Unable to solve for density at T = {} and P = {}", temp, press);
    }
    return m_waterIAPWS->coeffThermExp();
}

double WaterProps::isothermalCompressibility_IAPWS(double temp, double press)
{
    double dens = m_waterIAPWS->density(temp, press, WATER_LIQUID);
    if (dens < 0.0) {
        throw CanteraError("WaterProps::isothermalCompressibility_IAPWS",
            "Unable to solve for density at T = {} and P = {}", temp, press);
    }
    return m_waterIAPWS->isothermalCompressibility();
}

static const double H[4] = {1., 0.978197, 0.579829, -0.202354};

static const double Hij[6][7] = {
    { 0.5132047, 0.2151778, -0.2818107,  0.1778064, -0.04176610,          0.,           0.},
    { 0.3205656, 0.7317883, -1.070786 ,  0.4605040,          0., -0.01578386,           0.},
    { 0.,        1.241044 , -1.263184 ,  0.2340379,          0.,          0.,           0.},
    { 0.,        1.476783 ,         0., -0.4924179,   0.1600435,          0., -0.003629481},
    {-0.7782567,      0.0 ,         0.,  0.       ,          0.,          0.,           0.},
    { 0.1885447,      0.0 ,         0.,  0.       ,          0.,          0.,           0.},
};

static const double rhoStar = 317.763; // kg / m3
static const double presStar = 22.115E6; // Pa

double WaterProps::viscosityWater() const
{
    static const double TStar = 647.27; // Kelvin
    static const double muStar = 55.071E-6; //Pa s
    double temp = m_waterIAPWS->temperature();
    double dens = m_waterIAPWS->density();

    double rhobar = dens/rhoStar;
    double tbar = temp / TStar;
    double tbar2 = tbar * tbar;
    double tbar3 = tbar2 * tbar;
    double mu0bar = std::sqrt(tbar) / (H[0] + H[1]/tbar + H[2]/tbar2 + H[3]/tbar3);

    double tfac1 = 1.0 / tbar - 1.0;
    double tfac2 = tfac1 * tfac1;
    double tfac3 = tfac2 * tfac1;
    double tfac4 = tfac3 * tfac1;
    double tfac5 = tfac4 * tfac1;

    double rfac1 = rhobar - 1.0;
    double rfac2 = rfac1 * rfac1;
    double rfac3 = rfac2 * rfac1;
    double rfac4 = rfac3 * rfac1;
    double rfac5 = rfac4 * rfac1;
    double rfac6 = rfac5 * rfac1;

    double sum = (Hij[0][0]       + Hij[1][0]*tfac1       + Hij[4][0]*tfac4       + Hij[5][0]*tfac5 +
                      Hij[0][1]*rfac1 + Hij[1][1]*tfac1*rfac1 + Hij[2][1]*tfac2*rfac1 + Hij[3][1]*tfac3*rfac1 +
                      Hij[0][2]*rfac2 + Hij[1][2]*tfac1*rfac2 + Hij[2][2]*tfac2*rfac2 +
                      Hij[0][3]*rfac3 + Hij[1][3]*tfac1*rfac3 + Hij[2][3]*tfac2*rfac3 + Hij[3][3]*tfac3*rfac3 +
                      Hij[0][4]*rfac4 + Hij[3][4]*tfac3*rfac4 +
                      Hij[1][5]*tfac1*rfac5 + Hij[3][6]*tfac3*rfac6
                     );
    double mu1bar = std::exp(rhobar * sum);

    // Apply the near-critical point corrections if necessary
    double mu2bar = 1.0;
    if (tbar >= 0.9970 && tbar <= 1.0082 && rhobar >= 0.755 && rhobar <= 1.290) {
        double drhodp = 1.0 / m_waterIAPWS->dpdrho();
        drhodp *= presStar / rhoStar;
        double xsi = rhobar * drhodp;
        if (xsi >= 21.93) {
            mu2bar = 0.922 * std::pow(xsi, 0.0263);
        }
    }

    double mubar = mu0bar * mu1bar * mu2bar;
    return mubar * muStar;
}

double WaterProps::thermalConductivityWater() const
{
    static const double Tstar = 647.27;
    static const double rhostar = 317.763;
    static const double lambdastar = 0.4945;
    static const double presstar = 22.115E6;
    static const double L[4] = {
        1.0000,
        6.978267,
        2.599096,
        -0.998254
    };
    static const double Lji[6][5] = {
        { 1.3293046,    1.7018363,   5.2246158,   8.7127675, -1.8525999},
        {-0.40452437,  -2.2156845, -10.124111,   -9.5000611,  0.93404690},
        { 0.24409490,   1.6511057,   4.9874687,   4.3786606,  0.0},
        { 0.018660751, -0.76736002, -0.27297694, -0.91783782, 0.0},
        {-0.12961068,   0.37283344, -0.43083393,  0.0,        0.0},
        { 0.044809953, -0.11203160,  0.13333849,  0.0,        0.0},
    };

    double temp = m_waterIAPWS->temperature();
    double dens = m_waterIAPWS->density();

    double rhobar = dens/rhostar;
    double tbar = temp / Tstar;
    double tbar2 = tbar * tbar;
    double tbar3 = tbar2 * tbar;
    double lambda0bar = sqrt(tbar) / (L[0] + L[1]/tbar + L[2]/tbar2 + L[3]/tbar3);

    double tfac1 = 1.0 / tbar - 1.0;
    double tfac2 = tfac1 * tfac1;
    double tfac3 = tfac2 * tfac1;
    double tfac4 = tfac3 * tfac1;

    double rfac1 = rhobar - 1.0;
    double rfac2 = rfac1 * rfac1;
    double rfac3 = rfac2 * rfac1;
    double rfac4 = rfac3 * rfac1;
    double rfac5 = rfac4 * rfac1;

    double sum = (Lji[0][0]       + Lji[0][1]*tfac1        + Lji[0][2]*tfac2       + Lji[0][3]*tfac3       + Lji[0][4]*tfac4       +
                      Lji[1][0]*rfac1 + Lji[1][1]*tfac1*rfac1  + Lji[1][2]*tfac2*rfac1 + Lji[1][3]*tfac3*rfac1 + Lji[1][4]*tfac4*rfac1 +
                      Lji[2][0]*rfac2 + Lji[2][1]*tfac1*rfac2  + Lji[2][2]*tfac2*rfac2 + Lji[2][3]*tfac3*rfac2 +
                      Lji[3][0]*rfac3 + Lji[3][1]*tfac1*rfac3  + Lji[3][2]*tfac2*rfac3 + Lji[3][3]*tfac3*rfac3 +
                      Lji[4][0]*rfac4 + Lji[4][1]*tfac1*rfac4  + Lji[4][2]*tfac2*rfac4 +
                      Lji[5][0]*rfac5 + Lji[5][1]*tfac1*rfac5  + Lji[5][2]*tfac2*rfac5
                     );
    double lambda1bar = exp(rhobar * sum);
    double mu0bar = std::sqrt(tbar) / (H[0] + H[1]/tbar + H[2]/tbar2 + H[3]/tbar3);
    double tfac5 = tfac4 * tfac1;
    double rfac6 = rfac5 * rfac1;

    sum = (Hij[0][0]       + Hij[1][0]*tfac1       + Hij[4][0]*tfac4       + Hij[5][0]*tfac5 +
           Hij[0][1]*rfac1 + Hij[1][1]*tfac1*rfac1 + Hij[2][1]*tfac2*rfac1 + Hij[3][1]*tfac3*rfac1 +
           Hij[0][2]*rfac2 + Hij[1][2]*tfac1*rfac2 + Hij[2][2]*tfac2*rfac2 +
           Hij[0][3]*rfac3 + Hij[1][3]*tfac1*rfac3 + Hij[2][3]*tfac2*rfac3 + Hij[3][3]*tfac3*rfac3 +
           Hij[0][4]*rfac4 + Hij[3][4]*tfac3*rfac4 +
           Hij[1][5]*tfac1*rfac5 + Hij[3][6]*tfac3*rfac6
          );
    double mu1bar = std::exp(rhobar * sum);
    double t2r2 = tbar * tbar / (rhobar * rhobar);
    double drhodp = 1.0 / m_waterIAPWS->dpdrho();
    drhodp *= presStar / rhoStar;
    double xsi = rhobar * drhodp;
    double xsipow = std::pow(xsi, 0.4678);
    double rho1 = rhobar - 1.;
    double rho2 = rho1 * rho1;
    double rho4 = rho2 * rho2;
    double temp2 = (tbar - 1.0) * (tbar - 1.0);

    //     beta = M / (rho * Rgas) (d (pressure) / dT) at constant rho
    //
    // Note for ideal gases this is equal to one.
    //
    //     beta = delta (phi0_d() + phiR_d())
    //             - tau delta (phi0_dt() + phiR_dt())
    double beta = m_waterIAPWS->coeffPresExp();
    double dpdT_const_rho = beta * GasConstant * dens / 18.015268;
    dpdT_const_rho *= Tstar / presstar;
    double lambda2bar = 0.0013848 / (mu0bar * mu1bar) * t2r2 * dpdT_const_rho * dpdT_const_rho *
                             xsipow * sqrt(rhobar) * exp(-18.66*temp2 - rho4);
    return (lambda0bar * lambda1bar + lambda2bar) * lambdastar;
}

}
