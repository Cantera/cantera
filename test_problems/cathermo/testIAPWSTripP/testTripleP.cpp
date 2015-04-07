
#include "cantera/thermo/WaterPropsIAPWS.h"

#include <cstdio>
#include <cmath>

using namespace std;
using namespace Cantera;

int main()
{
#ifdef _MSC_VER
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    double dens, u, s, h, cv, cp, pres;
    WaterPropsIAPWS* water = new WaterPropsIAPWS();

    double T;

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
        printf("intEng (liquid) = ~0.0 J/kmol (fabs(u) < 5.0E-7)\n");
    } else {
        printf("intEng (liquid) = %g J/kmol\n", u);
    }

    s = water->entropy();
    if (fabs(s) < 5.0E-9) {
        printf("S (liquid) = ~0.0 J/kmolK (fabs(s) < 5.0E-9)\n");
    } else {
        printf("S (liquid) = %g J/kmolK\n", s);
    }

    h = water->enthalpy();
    printf("h (liquid) = %g J/kmol\n", h);
    printf("h (liquid) = %g J/kg\n", (h)/18.015268);

    double g_liq = water->Gibbs();
    printf("g (liquid) = %g J/kmol\n", g_liq);


    cv = water->cv();
    printf("cv (liquid) = %g J/kmolK\n", cv);

    cp = water->cp();
    printf("cp (liquid) = %g J/kmolK\n", cp);


    dens = water->density(T, pres, WATER_GAS);
    printf("dens (gas)    = %g kg m-3\n", dens);


    u = water->intEnergy();
    printf("intEng (gas) = %g J/kmol\n", u);

    s = water->entropy();
    printf("S (gas) = %g J/kmolK\n", s);

    h = water->enthalpy();
    printf("h (gas) = %g J/kmol\n", h);
    printf("h (gas) = %g J/kg\n", (h)/18.015268);

    double g_gas = water->Gibbs();
    printf("g (gas) = %g J/kmol\n", g_gas);


    cv = water->cv();
    printf("cv (gas) = %g J/kmolK\n", cv);

    cp = water->cp();
    printf("cp (gas) = %g J/kmolK\n", cp);

    printf("\n");
    double deltaG = g_liq - g_gas;
    if (fabs(deltaG) < 1.0E-5) {
        printf("Delta g = ~0.0 J/kmol ( < 1.0E-5)\n");
    } else {
        printf("Delta g = %g J/kmol\n", g_liq - g_gas);
    }
    delete water;
    return 0;
}
