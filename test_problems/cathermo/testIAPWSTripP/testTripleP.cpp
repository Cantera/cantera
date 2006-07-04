
#include "stdio.h"
#include "WaterPropsIAPWS.h"
#include <new>
using namespace std;

int main () {

    double dens, u, s, h, cv, cp, pres;
    WaterPropsIAPWS *water = new WaterPropsIAPWS();

    double T;
    double rho;


    /* 
     * Print out the triple point conditions
     */
    T = 273.16;
    pres = water->psat(T);
    printf("psat(%g) = %g\n", T, pres);

    dens = water->density(T, pres, WATER_LIQUID);
    printf("dens (liquid) = %g kg m-3\n", dens);

    u = water->intEnergy(T, dens);
    printf("intEng (liquid) = %g J/kmol\n", u);

    s = water->entropy(T, dens);
    printf("S (liquid) = %g J/kmolK\n", s);

    h = water->enthalpy(T, dens);
    printf("h (liquid) = %g J/kmol\n", h);
    printf("h (liquid) = %g J/kg\n", (h)/18.015268);

    double g_liq = water->Gibbs(T, dens);
    printf("g (liquid) = %g J/kmol\n", g_liq);


    cv = water->cv(T, dens);
    printf("cv (liquid) = %g J/kmolK\n", cv);

    cp = water->cp(T, dens);
    printf("cp (liquid) = %g J/kmolK\n", cp);


    dens = water->density(T, pres, WATER_GAS);
    printf("dens (gas)    = %g kg m-3\n", dens);

  
    u = water->intEnergy(T, dens);
    printf("intEng (gas) = %g J/kmol\n", u);

    s = water->entropy(T, dens);
    printf("S (gas) = %g J/kmolK\n", s);

    h = water->enthalpy(T, dens);
    printf("h (gas) = %g J/kmol\n", h);
    printf("h (gas) = %g J/kg\n", (h)/18.015268);

    double g_gas = water->Gibbs(T, dens);
    printf("g (gas) = %g J/kmol\n", g_gas);


    cv = water->cv(T, dens);
    printf("cv (gas) = %g J/kmolK\n", cv);

    cp = water->cp(T, dens);
    printf("cp (gas) = %g J/kmolK\n", cp);

    printf("\n");
    printf("Delta g = %g J/kmol\n", g_liq - g_gas);
 

    delete water;
    return 0;
}
