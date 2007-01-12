/*
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 * 
 */

#ifdef SRCDIRTREE
#include "ct_defs.h"
#include "ThermoPhase.h"
#include "IdealGasMix.h"
#include "equil.h"
#else
#include "Cantera.h"
#include "IdealGasMix.h"
#include "equilibrium.h"
#endif


using namespace std;
using namespace Cantera;
int main(int argc, char **argv) {
#ifdef DEBUG_HKM
  ChemEquil_print_lvl = 0;
#endif
  try {

    IdealGasPhase* gas = new IdealGasMix("air_below6000K.xml","air_below6000K");

    vector_fp IndVar2(6, 0.0);
    IndVar2[0] = 1.5E5;
    IndVar2[1] = 3.0E5;
    IndVar2[2] = 9.0E5;
    IndVar2[3] = 2.7E6;
    IndVar2[4] = 6.7E6;
    IndVar2[5] = 1.0E7;

    vector_fp IndVar1(7, 0.0);
    IndVar1[0] = 1.0E-8;
    IndVar1[1] = 1.0E-7;
    IndVar1[2] = 1.0E-6;
    IndVar1[3] = 1.0E-5;
    IndVar1[4] = 1.0E-4;
    IndVar1[5] = 1.0E-3;
    IndVar1[6] = 1.0E-2;
    int nj = 6;
    int ni = 7;

    for (int j=0; j<nj; j++) {
      for (int i=0; i<ni; i++) {
	double offset = -301471.39;
	gas->setState_UV(IndVar2[j]+offset,1.0/IndVar1[i]);
	double tkelvin = gas->temperature();
	double pres = gas->pressure();
	printf("Initial T = %g, pres = %g atm \n", tkelvin, pres/OneAtm);
	beginLogGroup("topEquil", -1);
	equilibrate(*gas,"UV", -1);
	endLogGroup("topEquil");
	cout << report(*gas) << endl;

	tkelvin = gas->temperature();
	pres = gas->pressure();
	printf("Final T = %g, pres = %g atm\n", tkelvin, pres/OneAtm);
   
      }
    }
    delete gas;
  }

  catch (CanteraError) {
    showErrors();
  }
  return 0;
}
