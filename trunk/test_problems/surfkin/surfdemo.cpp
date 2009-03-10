/*
 * $Id: surfdemo.cpp,v 1.7 2006/06/20 22:36:12 hkmoffa Exp $
 *
 *  Sample program that solves an implicit problem for surface
 *  site fractions.
 */
#ifdef SRCDIRTREE
#include "ct_defs.h"
#include "../silane_equil/IdealGasMix.h"
#include "Interface.h"
#else
#include "Cantera.h"
#include "IdealGasMix.h"
#include "Interface.h"
#endif
using namespace Cantera;

int main() {

    try {
        IdealGasMix gas("gri30.xml", "gri30");
        gas.setState_TPX(1200.0, OneAtm, 
            "H2:2, O2:1, OH:0.01, H:0.01, O:0.01");

        vector<ThermoPhase*> phases;
        phases.push_back(&gas);
        Interface surf("surface.xml", "surface", phases);
        vector_fp cov;
        cov.push_back(0.8);
        cov.push_back(0.2);
        surf.setCoverages(DATA_PTR(cov));
        vector_fp wdot(gas.nSpecies() + surf.nSpecies());
        surf.getNetProductionRates(DATA_PTR(wdot));
        int k;
        for (k = 0; k < gas.nSpecies(); k++)
            cout << gas.speciesName(k) << "   " << wdot[k] << endl;
        for (k = 0; k < surf.nSpecies(); k++)
            cout << surf.speciesName(k) << "   " 
                 << wdot[k+gas.nSpecies()] << endl;
        
    }
    catch (CanteraError) {
        showErrors(cout);
    }
    return 0;
}
      
