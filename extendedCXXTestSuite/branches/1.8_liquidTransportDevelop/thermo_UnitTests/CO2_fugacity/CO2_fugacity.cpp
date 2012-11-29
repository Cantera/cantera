/*
 * $Id: testPureWater.cpp 255 2009-11-09 23:36:49Z hkmoffa $
 */

#include "cantera/Cantera.h"
#include "cantera/PureFluid.h"

#include <new>
#include <cstdio>
#include <cmath>

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

      XML_Node *xc = get_XML_File("liquidvapor.xml"); 
      XML_Node * const xs = xc->findNameID("phase", "carbondioxide");
      ThermoPhase *water_tp = newPhase(*xs);
      PureFluidPhase *w = dynamic_cast <PureFluidPhase *>(water_tp);

      /* 
       * Print out the triple point conditions
       */
      double temp = 216.54;
      pres = w->satPressure(temp);
      printf("psat(%g) = %g\n", temp, pres);
      double presLow = 1.0E-2;
      temp = 298.15;
      double oneBar = 1.0E5;
      //double oneBar = 1.01325E5;
 

      // w->setStdState(-393.52E6, 213.79E3, 298.15, 1.0E5);

      printf("Comparisons to NIST: (see http://webbook.nist.gov):\n\n");

      w->setDensity(1.0E-8);
      w->setState_TP(temp, presLow);
      double h = w->enthalpy_mole();
      printf("H0(298.15) = %g J/kmol\n", h);

      double s = w->entropy_mole();
      s -= GasConstant * log(oneBar/presLow);
      printf("S0(298.15) = %g J/kmolK\n", s);


      printf("\n\n");


      double P[30];
      int n = 0;
      temp = 303.;

      double psat = w->satPressure(temp);

      P[n++] = 100.;
      P[n++] = 1.0E5;
      P[n++] = OneAtm;
      
      P[n++] = psat * 0.1;
      P[n++] = psat * 0.2;
      P[n++] = psat * 0.3;
      P[n++] = psat * 0.4;
      P[n++] = psat * 0.5;
      P[n++] = psat * 0.6;
      P[n++] = psat * 0.7;
      P[n++] = psat * 0.8;
      P[n++] = psat * 0.9;
      P[n++] = psat * 1.0;

      printf(" T    = %g\n", temp);
      printf(" psat = %g\n", psat);

      double plow = 100.;
      w->setState_TP(temp, plow);
      double glow =  w->gibbs_mole();
      double RT = GasConstant * temp;

      printf ("           press   chemPot        FugCoeff       Fugacity     \n");
      printf ("           Pascal  J/kmol/S                      (Pascal)      \n");
      for (int i = 0; i < n; i++) {
	
	w->setState_TP(temp, P[i]);
	double g = w->gibbs_mole();
	
	double lnFugRatio = (g - glow) / RT;
	double FugRatio =  exp(lnFugRatio)*plow/P[i];
	double fug = P[i] * FugRatio;

	printf("   %12.4g %12.4g  %12.4g   %12.4g \n", P[i], g, FugRatio, fug);

      }


   
      delete w;  
    } catch (CanteraError) {

      showErrors();
      Cantera::appdelete();
      return -1;
    }


    return 0;
}
