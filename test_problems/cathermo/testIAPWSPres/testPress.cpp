
#include "stdio.h"
#include "math.h"
#include "WaterPropsIAPWS.h"
#include <new>
using namespace std;

int main () {

    double dens, u, s, h;
    WaterPropsIAPWS *water = new WaterPropsIAPWS();

    double T = 273.15 + 100.;
    double rho = 10125. * 18.01 / (8.314472E3 * T);

    double pres = water->pressure(T, rho);
    printf("pres = %g\n", pres);

    /* 
     * Print out the triple point conditions
     */
    T = 273.16;
    pres = water->psat(T);
    printf("psat(%g) = %g\n", T, pres);

    dens = water->density(T, pres, WATER_LIQUID);
    printf("dens (liquid) = %g kg m-3\n", dens);

    u = water->intEnergy(T, dens);
    if (fabs(u) < 5.0E-7) {
      printf("intEng (liquid) ~= 0.0 J/kmol (less than fabs(5.0E-7))\n");
    } else {
      printf("intEng (liquid) = %g J/kmol\n", u);
    }

    s = water->entropy(T, dens);
    if (fabs(s) < 1.0E-9) {
      printf("S (liquid) ~= 0.0 J/kmolK (less than fabs(1.0E-9))\n");
    } else {
      printf("S (liquid) = %g J/kmolK\n", s);
    }

    h = water->enthalpy(T, dens);
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

    dens = water->density(T, pres, WATER_GAS);
    printf("dens (gas)    = %g kg m-3\n", dens);

    /*
     * Calculate a few test points for the estimated
     * saturation pressure function
     */
    T =	273.15 + 0.;
    pres = water->psat_est(T);
    printf("psat_est(%g) = %g\n", T, pres);
    
    T =	 313.9999;
    pres = water->psat_est(T);
    printf("psat_est(%g) = %g\n", T, pres);
    
    T =	314.0001;
    pres = water->psat_est(T);
    printf("psat_est(%g) = %g\n", T, pres);
    
    T =	273.15 + 100.;
    pres = water->psat_est(T);
    printf("psat_est(%g) = %g\n", T, pres);
    
    T =	647.25;
    pres = water->psat_est(T);
    printf("psat_est(%g) = %g\n", T, pres);







    delete water;
    return 0;
}
