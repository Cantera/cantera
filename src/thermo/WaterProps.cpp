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

WaterProps::WaterProps(const WaterProps& b)  :
    m_waterIAPWS(0),
    m_own_sub(false)
{
    *this = b;
}

WaterProps::~WaterProps()
{
    if (m_own_sub) {
        delete m_waterIAPWS;
    }
}

WaterProps& WaterProps::operator=(const WaterProps& b)
{
    if (&b == this) {
        return *this;
    }

    if (m_own_sub) {
        delete m_waterIAPWS;
    }
    if (b.m_own_sub) {
        m_waterIAPWS = new WaterPropsIAPWS();
        m_own_sub = true;
    } else {
        m_waterIAPWS = b.m_waterIAPWS;
        m_own_sub = false;
    }

    return *this;
}

doublereal WaterProps::density_T(doublereal T, doublereal P, int ifunc)
{
    static const doublereal Tc = T - 273.15;
    static const doublereal U1 = 288.9414;
    static const doublereal U2 = 508929.2;
    static const doublereal U3 = 68.12963;
    static const doublereal U4 = -3.9863;

    doublereal tmp1 = Tc + U1;
    doublereal tmp4 = Tc + U4;
    doublereal t4t4 = tmp4 * tmp4;
    doublereal tmp3 = Tc + U3;
    doublereal rho = 1000. * (1.0 - tmp1*t4t4/(U2 * tmp3));

    // Impose an ideal gas lower bound on rho. We need this to ensure positivity
    // of rho, even though it is grossly unrepresentative.
    doublereal rhomin = P / (GasConstant * T);
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
        doublereal drhodT = 1000./U2 * (
                                - tmp4 * tmp4 / (tmp3)
                                - tmp1 * 2 * tmp4 / (tmp3)
                                + tmp1 * t4t4 / (tmp3*tmp3)
                            );
        return drhodT;
    } else if (ifunc == 3) {
        return 0.0;
    } else if (ifunc == 2) {
        doublereal t3t3 = tmp3 * tmp3;
        doublereal d2rhodT2 = 1000./U2 *
                               ((-4.0*tmp4-2.0*tmp1)/tmp3 +
                                (2.0*t4t4 + 4.0*tmp1*tmp4)/t3t3
                                - 2.0*tmp1 * t4t4/(t3t3*tmp3));
        return d2rhodT2;
    }
    return rho;
}

doublereal WaterProps::relEpsilon(doublereal T, doublereal P_pascal,
                                  int ifunc)
{
    static const doublereal U1 = 3.4279E2;
    static const doublereal U2 = -5.0866E-3;
    static const doublereal U3 = 9.4690E-7;
    static const doublereal U4 = -2.0525;
    static const doublereal U5 = 3.1159E3;
    static const doublereal U6 = -1.8289E2;
    static const doublereal U7 = -8.0325E3;
    static const doublereal U8 = 4.2142E6;
    static const doublereal U9 = 2.1417;
    doublereal T2 = T * T;

    doublereal eps1000 = U1 * exp(U2 * T + U3 * T2);
    doublereal C = U4 + U5/(U6 + T);
    doublereal B = U7 + U8/T + U9 * T;
    doublereal Pbar = P_pascal * 1.0E-5;
    doublereal tmpBpar = B + Pbar;
    doublereal tmpB1000 = B + 1000.0;
    doublereal ltmp = log(tmpBpar/tmpB1000);
    doublereal epsRel = eps1000 + C * ltmp;

    if (ifunc == 1 || ifunc == 2) {
        doublereal tmpC = U6 + T;
        doublereal dCdT = - U5/(tmpC * tmpC);
        doublereal dBdT = - U8/(T * T) + U9;
        doublereal deps1000dT = eps1000 * (U2 + 2.0 * U3 * T);
        doublereal dltmpdT = (dBdT/tmpBpar - dBdT/tmpB1000);
        if (ifunc == 1) {
            return deps1000dT + dCdT * ltmp + C * dltmpdT;
        }
        doublereal T3 = T2 * T;
        doublereal d2CdT2 = - 2.0 * dCdT / tmpC;
        doublereal d2BdT2 = 2.0 * U8 / (T3);
        doublereal d2ltmpdT2 = (d2BdT2*(1.0/tmpBpar - 1.0/tmpB1000) +
                                dBdT*dBdT*(1.0/(tmpB1000*tmpB1000) - 1.0/(tmpBpar*tmpBpar)));
        doublereal d2eps1000dT2 = (deps1000dT * (U2 + 2.0 * U3 * T) + eps1000 * (2.0 * U3));

        if (ifunc == 2) {
            doublereal d2epsReldT2 = (d2eps1000dT2 + d2CdT2 * ltmp + 2.0 * dCdT * dltmpdT
                                      + C * d2ltmpdT2);
            return d2epsReldT2;
        }
    }
    if (ifunc == 3) {
        doublereal dltmpdP = 1.0E-5 / tmpBpar;
        return C * dltmpdP;
    }
    return epsRel;
}

doublereal WaterProps::ADebye(doublereal T, doublereal P_input, int ifunc)
{
    doublereal psat = satPressure(T);
    doublereal P;
    if (psat > P_input) {
        P = psat;
    } else {
        P = P_input;
    }
    doublereal epsRelWater = relEpsilon(T, P, 0);
    doublereal epsilon = epsilon_0 * epsRelWater;
    doublereal dw = density_IAPWS(T, P);
    doublereal tmp = sqrt(2.0 * Avogadro * dw / 1000.);
    doublereal tmp2 = ElectronCharge * ElectronCharge * Avogadro /
                      (epsilon * GasConstant * T);
    doublereal tmp3 = tmp2 * sqrt(tmp2);
    doublereal A_Debye = tmp * tmp3 / (8.0 * Pi);

    // dAdT = - 3/2 Ad/T + 1/2 Ad/dw d(dw)/dT - 3/2 Ad/eps d(eps)/dT
    // dAdT = - 3/2 Ad/T - 1/2 Ad/Vw d(Vw)/dT - 3/2 Ad/eps d(eps)/dT
    if (ifunc == 1 || ifunc == 2) {
        doublereal dAdT = - 1.5 * A_Debye / T;

        doublereal depsRelWaterdT = relEpsilon(T, P, 1);
        dAdT -= A_Debye * (1.5 * depsRelWaterdT / epsRelWater);

        // calculate d(lnV)/dT _constantP, i.e., the cte
        doublereal cte = coeffThermalExp_IAPWS(T, P);
        doublereal contrib2 = - A_Debye * (0.5 * cte);
        dAdT += contrib2;

        if (ifunc == 1) {
            return dAdT;
        }

        if (ifunc == 2) {
            // Get the second derivative of the dielectric constant wrt T
            // -> we will take each of the terms in dAdT and differentiate
            //    it again.
            doublereal d2AdT2 = 1.5 / T * (A_Debye/T - dAdT);
            doublereal d2epsRelWaterdT2 = relEpsilon(T, P, 2);
            d2AdT2 += 1.5 * (- dAdT * depsRelWaterdT / epsRelWater
                             - A_Debye / epsRelWater *
                             (d2epsRelWaterdT2 - depsRelWaterdT * depsRelWaterdT / epsRelWater));
            doublereal deltaT = -0.1;
            doublereal Tdel = T + deltaT;
            doublereal cte_del = coeffThermalExp_IAPWS(Tdel, P);
            doublereal dctedT = (cte_del - cte) / Tdel;
            doublereal contrib3 = 0.5 * (-(dAdT * cte) -(A_Debye * dctedT));
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
        doublereal dAdP = 0.0;
        doublereal depsRelWaterdP = relEpsilon(T, P, 3);
        dAdP -= A_Debye * (1.5 * depsRelWaterdP / epsRelWater);
        doublereal kappa = isothermalCompressibility_IAPWS(T,P);
        dAdP += A_Debye * (0.5 * kappa);
        return dAdP;
    }
    return A_Debye;
}

doublereal WaterProps::satPressure(doublereal T)
{
    return m_waterIAPWS->psat(T);
}

doublereal WaterProps::density_IAPWS(doublereal temp, doublereal press)
{
    return m_waterIAPWS->density(temp, press, WATER_LIQUID);
}

doublereal WaterProps::density_IAPWS() const
{
    return m_waterIAPWS->density();
}

doublereal WaterProps::coeffThermalExp_IAPWS(doublereal temp, doublereal press)
{
    doublereal dens = m_waterIAPWS->density(temp, press, WATER_LIQUID);
    if (dens < 0.0) {
        throw CanteraError("WaterProps::coeffThermalExp_IAPWS",
            "Unable to solve for density at T = {} and P = {}", temp, press);
    }
    return m_waterIAPWS->coeffThermExp();
}

doublereal WaterProps::isothermalCompressibility_IAPWS(doublereal temp, doublereal press)
{
    doublereal dens = m_waterIAPWS->density(temp, press, WATER_LIQUID);
    if (dens < 0.0) {
        throw CanteraError("WaterProps::isothermalCompressibility_IAPWS",
            "Unable to solve for density at T = {} and P = {}", temp, press);
    }
    return m_waterIAPWS->isothermalCompressibility();
}

static const doublereal H[4] = {1., 0.978197, 0.579829, -0.202354};

static const doublereal Hij[6][7] = {
    { 0.5132047, 0.2151778, -0.2818107,  0.1778064, -0.04176610,          0.,           0.},
    { 0.3205656, 0.7317883, -1.070786 ,  0.4605040,          0., -0.01578386,           0.},
    { 0.,        1.241044 , -1.263184 ,  0.2340379,          0.,          0.,           0.},
    { 0.,        1.476783 ,         0., -0.4924179,   0.1600435,          0., -0.003629481},
    {-0.7782567,      0.0 ,         0.,  0.       ,          0.,          0.,           0.},
    { 0.1885447,      0.0 ,         0.,  0.       ,          0.,          0.,           0.},
};

static const doublereal rhoStar = 317.763; // kg / m3
static const doublereal presStar = 22.115E6; // Pa

doublereal WaterProps::viscosityWater() const
{
    static const doublereal TStar = 647.27; // Kelvin
    static const doublereal muStar = 55.071E-6; //Pa s
    doublereal temp = m_waterIAPWS->temperature();
    doublereal dens = m_waterIAPWS->density();

    doublereal rhobar = dens/rhoStar;
    doublereal tbar = temp / TStar;
    doublereal tbar2 = tbar * tbar;
    doublereal tbar3 = tbar2 * tbar;
    doublereal mu0bar = std::sqrt(tbar) / (H[0] + H[1]/tbar + H[2]/tbar2 + H[3]/tbar3);

    doublereal tfac1 = 1.0 / tbar - 1.0;
    doublereal tfac2 = tfac1 * tfac1;
    doublereal tfac3 = tfac2 * tfac1;
    doublereal tfac4 = tfac3 * tfac1;
    doublereal tfac5 = tfac4 * tfac1;

    doublereal rfac1 = rhobar - 1.0;
    doublereal rfac2 = rfac1 * rfac1;
    doublereal rfac3 = rfac2 * rfac1;
    doublereal rfac4 = rfac3 * rfac1;
    doublereal rfac5 = rfac4 * rfac1;
    doublereal rfac6 = rfac5 * rfac1;

    doublereal sum = (Hij[0][0]       + Hij[1][0]*tfac1       + Hij[4][0]*tfac4       + Hij[5][0]*tfac5 +
                      Hij[0][1]*rfac1 + Hij[1][1]*tfac1*rfac1 + Hij[2][1]*tfac2*rfac1 + Hij[3][1]*tfac3*rfac1 +
                      Hij[0][2]*rfac2 + Hij[1][2]*tfac1*rfac2 + Hij[2][2]*tfac2*rfac2 +
                      Hij[0][3]*rfac3 + Hij[1][3]*tfac1*rfac3 + Hij[2][3]*tfac2*rfac3 + Hij[3][3]*tfac3*rfac3 +
                      Hij[0][4]*rfac4 + Hij[3][4]*tfac3*rfac4 +
                      Hij[1][5]*tfac1*rfac5 + Hij[3][6]*tfac3*rfac6
                     );
    doublereal mu1bar = std::exp(rhobar * sum);

    // Apply the near-critical point corrections if necessary
    doublereal mu2bar = 1.0;
    if (tbar >= 0.9970 && tbar <= 1.0082 && rhobar >= 0.755 && rhobar <= 1.290) {
        doublereal drhodp = 1.0 / m_waterIAPWS->dpdrho();
        drhodp *= presStar / rhoStar;
        doublereal xsi = rhobar * drhodp;
        if (xsi >= 21.93) {
            mu2bar = 0.922 * std::pow(xsi, 0.0263);
        }
    }

    doublereal mubar = mu0bar * mu1bar * mu2bar;
    return mubar * muStar;
}

doublereal WaterProps::thermalConductivityWater() const
{
    static const doublereal Tstar = 647.27;
    static const doublereal rhostar = 317.763;
    static const doublereal lambdastar = 0.4945;
    static const doublereal presstar = 22.115E6;
    static const doublereal L[4] = {
        1.0000,
        6.978267,
        2.599096,
        -0.998254
    };
    static const doublereal Lji[6][5] = {
        { 1.3293046,    1.7018363,   5.2246158,   8.7127675, -1.8525999},
        {-0.40452437,  -2.2156845, -10.124111,   -9.5000611,  0.93404690},
        { 0.24409490,   1.6511057,   4.9874687,   4.3786606,  0.0},
        { 0.018660751, -0.76736002, -0.27297694, -0.91783782, 0.0},
        {-0.12961068,   0.37283344, -0.43083393,  0.0,        0.0},
        { 0.044809953, -0.11203160,  0.13333849,  0.0,        0.0},
    };

    doublereal temp = m_waterIAPWS->temperature();
    doublereal dens = m_waterIAPWS->density();

    doublereal rhobar = dens/rhostar;
    doublereal tbar = temp / Tstar;
    doublereal tbar2 = tbar * tbar;
    doublereal tbar3 = tbar2 * tbar;
    doublereal lambda0bar = sqrt(tbar) / (L[0] + L[1]/tbar + L[2]/tbar2 + L[3]/tbar3);

    doublereal tfac1 = 1.0 / tbar - 1.0;
    doublereal tfac2 = tfac1 * tfac1;
    doublereal tfac3 = tfac2 * tfac1;
    doublereal tfac4 = tfac3 * tfac1;

    doublereal rfac1 = rhobar - 1.0;
    doublereal rfac2 = rfac1 * rfac1;
    doublereal rfac3 = rfac2 * rfac1;
    doublereal rfac4 = rfac3 * rfac1;
    doublereal rfac5 = rfac4 * rfac1;

    doublereal sum = (Lji[0][0]       + Lji[0][1]*tfac1        + Lji[0][2]*tfac2       + Lji[0][3]*tfac3       + Lji[0][4]*tfac4       +
                      Lji[1][0]*rfac1 + Lji[1][1]*tfac1*rfac1  + Lji[1][2]*tfac2*rfac1 + Lji[1][3]*tfac3*rfac1 + Lji[1][4]*tfac4*rfac1 +
                      Lji[2][0]*rfac2 + Lji[2][1]*tfac1*rfac2  + Lji[2][2]*tfac2*rfac2 + Lji[2][3]*tfac3*rfac2 +
                      Lji[3][0]*rfac3 + Lji[3][1]*tfac1*rfac3  + Lji[3][2]*tfac2*rfac3 + Lji[3][3]*tfac3*rfac3 +
                      Lji[4][0]*rfac4 + Lji[4][1]*tfac1*rfac4  + Lji[4][2]*tfac2*rfac4 +
                      Lji[5][0]*rfac5 + Lji[5][1]*tfac1*rfac5  + Lji[5][2]*tfac2*rfac5
                     );
    doublereal lambda1bar = exp(rhobar * sum);
    doublereal mu0bar = std::sqrt(tbar) / (H[0] + H[1]/tbar + H[2]/tbar2 + H[3]/tbar3);
    doublereal tfac5 = tfac4 * tfac1;
    doublereal rfac6 = rfac5 * rfac1;

    sum = (Hij[0][0]       + Hij[1][0]*tfac1       + Hij[4][0]*tfac4       + Hij[5][0]*tfac5 +
           Hij[0][1]*rfac1 + Hij[1][1]*tfac1*rfac1 + Hij[2][1]*tfac2*rfac1 + Hij[3][1]*tfac3*rfac1 +
           Hij[0][2]*rfac2 + Hij[1][2]*tfac1*rfac2 + Hij[2][2]*tfac2*rfac2 +
           Hij[0][3]*rfac3 + Hij[1][3]*tfac1*rfac3 + Hij[2][3]*tfac2*rfac3 + Hij[3][3]*tfac3*rfac3 +
           Hij[0][4]*rfac4 + Hij[3][4]*tfac3*rfac4 +
           Hij[1][5]*tfac1*rfac5 + Hij[3][6]*tfac3*rfac6
          );
    doublereal mu1bar = std::exp(rhobar * sum);
    doublereal t2r2 = tbar * tbar / (rhobar * rhobar);
    doublereal drhodp = 1.0 / m_waterIAPWS->dpdrho();
    drhodp *= presStar / rhoStar;
    doublereal xsi = rhobar * drhodp;
    doublereal xsipow = std::pow(xsi, 0.4678);
    doublereal rho1 = rhobar - 1.;
    doublereal rho2 = rho1 * rho1;
    doublereal rho4 = rho2 * rho2;
    doublereal temp2 = (tbar - 1.0) * (tbar - 1.0);

    //     beta = M / (rho * Rgas) (d (pressure) / dT) at constant rho
    //
    // Note for ideal gases this is equal to one.
    //
    //     beta = delta (phi0_d() + phiR_d())
    //             - tau delta (phi0_dt() + phiR_dt())
    doublereal beta = m_waterIAPWS->coeffPresExp();
    doublereal dpdT_const_rho = beta * GasConstant * dens / 18.015268;
    dpdT_const_rho *= Tstar / presstar;
    doublereal lambda2bar = 0.0013848 / (mu0bar * mu1bar) * t2r2 * dpdT_const_rho * dpdT_const_rho *
                             xsipow * sqrt(rhobar) * exp(-18.66*temp2 - rho4);
    return (lambda0bar * lambda1bar + lambda2bar) * lambdastar;
}

}
