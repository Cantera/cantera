//! @file WaterTransport.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/transport/WaterTransport.h"
#include "cantera/thermo/VPStandardStateTP.h"
#include "cantera/thermo/PDSS_Water.h"
#include "cantera/thermo/WaterSSTP.h"

using namespace std;

namespace {

const double Tstar = 647.27;
const double rhoStar = 317.763; // kg / m3
const double presStar = 22.115E6; // Pa
const double muStar = 55.071E-6; //Pa s

const double H[4] = {1., 0.978197, 0.579829, -0.202354};
const double Hij[6][7] = {
    { 0.5132047, 0.2151778, -0.2818107,  0.1778064, -0.04176610,          0.,           0.},
    { 0.3205656, 0.7317883, -1.070786 ,  0.4605040,          0., -0.01578386,           0.},
    { 0.,        1.241044 , -1.263184 ,  0.2340379,          0.,          0.,           0.},
    { 0.,        1.476783 ,         0., -0.4924179,   0.1600435,          0., -0.003629481},
    {-0.7782567,      0.0 ,         0.,  0.       ,          0.,          0.,           0.},
    { 0.1885447,      0.0 ,         0.,  0.       ,          0.,          0.,           0.},
};

}

namespace Cantera
{

WaterTransport::WaterTransport(thermo_t* thermo, int ndim) :
    Transport(thermo, ndim)
{
    if (thermo) {
        init(thermo);
    }
}

void WaterTransport::init(thermo_t* thermo, int mode, int log_level)
{
    m_thermo = thermo;
}

double WaterTransport::viscosity()
{
    static const double TStar = 647.27; // Kelvin
    double temp = m_thermo->temperature();
    double dens = m_thermo->density();

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

    double sum = Hij[0][0]       + Hij[1][0]*tfac1       + Hij[4][0]*tfac4       + Hij[5][0]*tfac5 +
                 Hij[0][1]*rfac1 + Hij[1][1]*tfac1*rfac1 + Hij[2][1]*tfac2*rfac1 + Hij[3][1]*tfac3*rfac1 +
                 Hij[0][2]*rfac2 + Hij[1][2]*tfac1*rfac2 + Hij[2][2]*tfac2*rfac2 +
                 Hij[0][3]*rfac3 + Hij[1][3]*tfac1*rfac3 + Hij[2][3]*tfac2*rfac3 + Hij[3][3]*tfac3*rfac3 +
                 Hij[0][4]*rfac4 + Hij[3][4]*tfac3*rfac4 +
                 Hij[1][5]*tfac1*rfac5 + Hij[3][6]*tfac3*rfac6;
    double mu1bar = std::exp(rhobar * sum);

    // Apply the near-critical point corrections if necessary
    double mu2bar = 1.0;
    if (tbar >= 0.9970 && tbar <= 1.0082 && rhobar >= 0.755 && rhobar <= 1.290) {
        double drhodp = m_thermo->isothermalCompressibility() * dens;
        drhodp *= presStar / rhoStar;
        double xsi = rhobar * drhodp;
        if (xsi >= 21.93) {
            mu2bar = 0.922 * std::pow(xsi, 0.0263);
        }
    }

    double mubar = mu0bar * mu1bar * mu2bar;
    return mubar * muStar;
}

double WaterTransport::thermalConductivity()
{
    static const double lambdastar = 0.4945;
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

    double temp = m_thermo->temperature();
    double dens = m_thermo->density();

    double rhobar = dens / rhoStar;
    double tbar = temp / Tstar;
    double tbar2 = tbar * tbar;
    double tbar3 = tbar2 * tbar;
    double lambda0bar = sqrt(tbar) / (L[0] + L[1]/tbar + L[2]/tbar2 + L[3]/tbar3);

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

    double sum = (Lji[0][0]       + Lji[0][1]*tfac1        + Lji[0][2]*tfac2       + Lji[0][3]*tfac3       + Lji[0][4]*tfac4       +
                      Lji[1][0]*rfac1 + Lji[1][1]*tfac1*rfac1  + Lji[1][2]*tfac2*rfac1 + Lji[1][3]*tfac3*rfac1 + Lji[1][4]*tfac4*rfac1 +
                      Lji[2][0]*rfac2 + Lji[2][1]*tfac1*rfac2  + Lji[2][2]*tfac2*rfac2 + Lji[2][3]*tfac3*rfac2 +
                      Lji[3][0]*rfac3 + Lji[3][1]*tfac1*rfac3  + Lji[3][2]*tfac2*rfac3 + Lji[3][3]*tfac3*rfac3 +
                      Lji[4][0]*rfac4 + Lji[4][1]*tfac1*rfac4  + Lji[4][2]*tfac2*rfac4 +
                      Lji[5][0]*rfac5 + Lji[5][1]*tfac1*rfac5  + Lji[5][2]*tfac2*rfac5
                     );
    double lambda1bar = exp(rhobar * sum);
    double mu0bar = std::sqrt(tbar) / (H[0] + H[1]/tbar + H[2]/tbar2 + H[3]/tbar3);

    sum = (Hij[0][0]       + Hij[1][0]*tfac1       + Hij[4][0]*tfac4       + Hij[5][0]*tfac5 +
           Hij[0][1]*rfac1 + Hij[1][1]*tfac1*rfac1 + Hij[2][1]*tfac2*rfac1 + Hij[3][1]*tfac3*rfac1 +
           Hij[0][2]*rfac2 + Hij[1][2]*tfac1*rfac2 + Hij[2][2]*tfac2*rfac2 +
           Hij[0][3]*rfac3 + Hij[1][3]*tfac1*rfac3 + Hij[2][3]*tfac2*rfac3 + Hij[3][3]*tfac3*rfac3 +
           Hij[0][4]*rfac4 + Hij[3][4]*tfac3*rfac4 +
           Hij[1][5]*tfac1*rfac5 + Hij[3][6]*tfac3*rfac6
          );
    double mu1bar = std::exp(rhobar * sum);
    double t2r2 = tbar2 / (rhobar * rhobar);
    double kappa = m_thermo->isothermalCompressibility();
    double xsi = rhobar * rhobar * kappa * presStar;
    double xsipow = std::pow(xsi, 0.4678);
    double temp2 = (tbar - 1.0) * (tbar - 1.0);
    double dpdT_const_rho = m_thermo->thermalExpansionCoeff() / kappa;
    dpdT_const_rho *= Tstar / presStar;
    double lambda2bar = 0.0013848 / (mu0bar * mu1bar) * t2r2 * dpdT_const_rho * dpdT_const_rho *
                             xsipow * sqrt(rhobar) * exp(-18.66*temp2 - rfac4);
    return (lambda0bar * lambda1bar + lambda2bar) * lambdastar;
}

}
