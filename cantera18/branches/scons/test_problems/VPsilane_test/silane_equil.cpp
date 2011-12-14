
/*
 *  $Author: hkmoffa $
 *  $Date: 2009/03/14 00:51:16 $
 *  $Revision: 1.2 $
 *
 *  Copyright 2002 California Institute of Technology
 *
 */

#include "Cantera.h"
#include "IdealGasMix.h"
#include "equilibrium.h"
#include "kernel/IdealSolnGasVPSS.h"
#include "kernel/ThermoFactory.h"
#include "kernel/PrintCtrl.h"

using namespace std;
using namespace Cantera;

int main(int argc, char **argv) {
  try {
    PrintCtrl::GlobalCrop = PrintCtrl::GCT_CROP;
    Cantera::IdealSolnGasVPSS gg("silane.xml", "silane");
     ThermoPhase *g = &gg; 
    //ThermoPhase *g = newPhase("silane.xml", "silane");
    g->setState_TPX(1500.0, 100.0, "SIH4:0.01, H2:0.99");
    //g.setState_TPX(1500.0, 1.0132E5, "SIH4:0.01, H2:0.99");
    Cantera::ChemEquil_print_lvl = 40;
    equilibrate(*g, "TP");
    std::string r = Cantera::report(*g, true);
    cout << r;
    cout << endl;
    return 0;
  }
  catch (CanteraError) {
    showErrors(cerr);
    cerr << "program terminating." << endl;
    return -1;
  }
}
