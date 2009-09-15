/*
 * $Id: testWaterPDSS.cpp,v 1.6 2008/12/17 17:31:13 hkmoffa Exp $
 */
#include "ct_defs.h"
#include "ctexceptions.h"
#include "global.h"
#include "xml.h"
#include "ctml.h"
#include "ThermoPhase.h"
#include "VPStandardStateTP.h"

#include "WaterSSTP.h"
#include "WaterTransport.h"
#include <new>

#include <cstdio>

using namespace std;
using namespace Cantera;

double tvalue(double val, double atol = 1.0E-9) {
    double rval = val;
    if (fabs(val) < atol) {
      rval = 0.0; 
    }
    return rval;
}

int main () {

  try {
   
    WaterSSTP * w = new WaterSSTP("waterTPphase.xml", "");

    WaterTransport *wtTran = new WaterTransport(w, 3);
    printf("-------------------------------------------------------------------------\n");
    printf("   T(C)    MPa     Phase         Visc     Visc(paper) \n");
    printf("                                   10-6 kg/m/s    \n");
    printf("-------------------------------------------------------------------------\n");
   
    double T = 273.15 + 25.0;
    double pres = 1.0E5;
    w->setState_TP(T, pres);
    double visc = wtTran->viscosity();
    printf("%8g %10.3g  L  %13.6g      890.5\n",
	   T - 273.15, pres * 1.0E-6, visc * 1.0E6);


    T = 273.15 + 100.0;
    pres = 1.0E5;
    w->setState_TP(T, pres);
    visc = wtTran->viscosity();
    printf("%8g %10.3g  L  %13.6g      281.9\n",
	   T - 273.15, pres * 1.0E-6, visc * 1.0E6);



    T = 273.15 + 100.0;
    pres = 1.0E7;
    w->setState_TP(T, pres);
    visc = wtTran->viscosity();
    printf("%8g %10.3g  L  %13.6g      284.5\n",
	   T - 273.15, pres * 1.0E-6, visc * 1.0E6);

    T = 273.15 + 250.0;
    pres = 5.0E6;
    w->setState_TP(T, pres);
    visc = wtTran->viscosity();
    printf("%8g %10.3g  L  %13.6g      106.4\n",
	   T - 273.15, pres * 1.0E-6, visc * 1.0E6);

    T = 273.15 + 250.0;
    pres = 5.0E7;
    w->setState_TP(T, pres);
    visc = wtTran->viscosity();
    printf("%8g %10.3g  L  %13.6g      117.5\n",
	   T - 273.15, pres * 1.0E-6, visc * 1.0E6);


    T = 273.15 + 350.0;
    pres = 1.75E7;
    w->setState_TP(T, pres);
    visc = wtTran->viscosity();
    printf("%8g %10.3g  L  %13.6g      67.0\n",
	   T - 273.15, pres * 1.0E-6, visc * 1.0E6);


    T = 273.15 + 400.0;
    pres = 1.50E7;
    w->setState_TP(T, pres);
    visc = wtTran->viscosity();
    printf("%8g %10.3g  SC %13.6g      24.93\n",
	   T - 273.15, pres * 1.0E-6, visc * 1.0E6);

    printf("-------------------------------------------------------------------------\n");

    delete w;  
  } catch (CanteraError) {

    showErrors();
    Cantera::appdelete();
    return -1;
  }


  return 0;
}
