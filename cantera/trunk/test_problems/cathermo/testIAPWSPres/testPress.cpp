
#include "cantera/thermo/WaterPropsIAPWS.h"

#include <cstdio>
#include <cmath>

using namespace std;
using namespace Cantera;

double numdpdt(WaterPropsIAPWS* water, double T, double pres)
{
    double rho = water->density(T, pres);
    water->setState_TR(T, rho);
    double presB = water->pressure();
    double Td = T + 0.001;
    water->setState_TR(Td, rho);
    double presd = water->pressure();
    return (presd - presB) / 0.001;
}

int main()
{
#ifdef _MSC_VER
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    double dens, u, s, h;
    WaterPropsIAPWS* water = new WaterPropsIAPWS();

    double T = 273.15 + 100.;
    double rho = 10125. * 18.01 / (8.314472E3 * T);

    water->setState_TR(T, rho);
    double pres = water->pressure();
    printf("pres = %g\n", pres);

    /*
     * Print out the triple point conditions
     */
    T = 273.16;
    pres = water->psat(T);
    printf("psat(%g) = %g\n", T, pres);

    dens = water->density(T, pres, WATER_LIQUID);
    printf("dens (liquid) = %g kg m-3\n", dens);

    u = water->intEnergy();
    if (fabs(u) < 5.0E-7) {
        printf("intEng (liquid) ~= 0.0 J/kmol (less than fabs(5.0E-7))\n");
    } else {
        printf("intEng (liquid) = %g J/kmol\n", u);
    }

    s = water->entropy();
    if (fabs(s) < 1.0E-9) {
        printf("S (liquid) ~= 0.0 J/kmolK (less than fabs(1.0E-9))\n");
    } else {
        printf("S (liquid) = %g J/kmolK\n", s);
    }

    h = water->enthalpy();
    printf("h (liquid) = %g J/kmol\n", h);
    printf("h (liquid) = %g J/kg\n", (h)/18.015268);


    dens = water->density(T, pres, WATER_GAS);
    printf("dens (gas)    = %g kg m-3\n", dens);



    /*
     * Print out the normal boiling point conditions
     */
    T = 373.124;
    pres = water->psat(T);
    printf("psat(%g) = %g\n", T, pres);

    dens = water->density(T, pres, WATER_LIQUID);
    printf("dens (liquid) = %g kg m-3\n", dens);

    double kappa = water->isothermalCompressibility();
    printf("kappa (liquid) = %20.13g kg m-3 \n", kappa);

    double pres2 = pres * 1.001;
    dens = water->density(T, pres2, WATER_LIQUID);
    kappa = water->isothermalCompressibility();
    printf("kappa (liquid) = %20.13g kg m-3 \n", kappa);

    dens = water->density(T, pres, WATER_GAS);
    printf("dens (gas)    = %g kg m-3\n", dens);

    kappa = water->isothermalCompressibility();
    printf("kappa (gas) = %20.13g kg m-3 \n", kappa);

    pres2 = pres * (0.999);
    dens = water->density(T, pres2, WATER_GAS);
    printf("dens (gas)    = %g kg m-3\n", dens);

    kappa = water->isothermalCompressibility();
    printf("kappa (gas) = %20.13g kg m-3 \n", kappa);

    /*
     * Calculate a few test points for the estimated
     * saturation pressure function
     */
    T = 273.15 + 0.;
    pres = water->psat_est(T);
    printf("psat_est(%g) = %g\n", T, pres);

    T = 313.9999;
    pres = water->psat_est(T);
    printf("psat_est(%g) = %g\n", T, pres);

    T = 314.0001;
    pres = water->psat_est(T);
    printf("psat_est(%g) = %g\n", T, pres);

    T = 273.15 + 100.;
    pres = water->psat_est(T);
    printf("psat_est(%g) = %g\n", T, pres);

    T = 647.25;
    pres = water->psat_est(T);
    printf("psat_est(%g) = %g\n", T, pres);


    T = 300;
    pres = 10.;
    rho = water->density(T, pres, WATER_GAS);
    water->setState_TR(T, rho);
    double beta = water->coeffPresExp();

    printf("beta    = %20.13g\n", beta);


    double dpdt = numdpdt(water, T, pres);

    rho = water->density(T,pres);
    double betaNum = dpdt * 18.015268 /(8.314371E3 * rho);
    printf("betaNum = %20.11g\n", betaNum);

    double alpha = water->coeffThermExp();
    printf("alpha = %20.13g\n", alpha);


    T = 300;
    pres = 10.E6;
    rho = water->density(T, pres, WATER_GAS);
    water->setState_TR(T, rho);
    beta = water->coeffPresExp();

    printf("beta    = %20.13g\n", beta);


    dpdt = numdpdt(water, T, pres);

    rho = water->density(T,pres);
    betaNum = dpdt * 18.015268 /(8.314371E3 * rho);

    printf("betaNum = %20.12g\n", betaNum);

    alpha = water->coeffThermExp();
    printf("alpha = %20.13g\n", alpha);

    T = 700;
    pres = 10.E6;
    rho = water->density(T, pres, WATER_GAS);
    water->setState_TR(T, rho);
    beta = water->coeffPresExp();

    printf("beta    = %20.13g\n", beta);


    dpdt = numdpdt(water, T, pres);

    rho = water->density(T,pres);
    betaNum = dpdt * 18.015268 /(8.314371E3 * rho);

    printf("betaNum = %20.10g\n", betaNum);

    alpha = water->coeffThermExp();
    printf("alpha = %20.13g\n", alpha);

    delete water;
    return 0;
}
