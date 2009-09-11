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

    double pres;
    try {
   
      WaterSSTP * w = new WaterSSTP("waterTPphase.xml", "");

      WaterTransport *wtTran = new WaterTransport(w, 3);

      //double T = 273.15 + 150.0;
      double T = 273.15 + 100.0;

      w->setState_TP(T, 1.0E7);
    
      double visc = wtTran->viscosity();

      printf("visc = %g\n", visc);

      delete w;  
    } catch (CanteraError) {

      showErrors();
      Cantera::appdelete();
      return -1;
    }


    return 0;
}
