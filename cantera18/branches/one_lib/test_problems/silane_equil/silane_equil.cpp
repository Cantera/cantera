
/*
 *  $Author: hkmoffa $
 *  $Date: 2008/08/23 00:55:58 $
 *  $Revision: 1.10 $
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
    IdealGasMix g("silane.xml", "silane");
    g.setState_TPX(1500.0, 100.0, "SIH4:0.01, H2:0.99");
    //g.setState_TPX(1500.0, 1.0132E5, "SIH4:0.01, H2:0.99");
    equilibrate(g, "TP");
    cout << g;
    return 0;
  }
  catch (CanteraError) {
    showErrors(cerr);
    cerr << "program terminating." << endl;
    return -1;
  }
}
