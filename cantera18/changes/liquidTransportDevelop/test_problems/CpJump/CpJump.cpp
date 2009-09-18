/*
 *  $Author: hkmoffa $
 *  $Date: 2008/02/16 21:46:15 $
 *  $Revision: 1.2 $
 *
 *  Copyright 2002 California Institute of Technology
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
using namespace Cantera_CXX;

int main(int argc, char **argv) {
  try {
    IdealGasMix g("bad_air.xml", "air");
    double pres = 1.0E5;
    g.setState_TPX(1000.1, pres, "O2:0.4, N2:0.6");
    equilibrate(g, "TP", -1);
    cout << g;
    double enth = g.enthalpy_mass();
    printf(" enth = %g\n", enth);
    enth -= 2.0E2;
    printf("attempted equil at (H,P) = %10.5g, %10.5g\n", enth, pres);
    g.setState_HP(enth, pres);
    equilibrate(g, "HP", -1);
    cout << g;


    return 0;
  }
  catch (CanteraError) {
    showErrors(cerr);
    cerr << "program terminating." << endl;
    return -1;
  }
}
