//
//     Replace this sample main program with your program
//
//

#include "Cantera.h"
#include "IdealGasMix.h"
#include "Interface.h"

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
        surf.setCoverages(cov.begin());
        vector_fp wdot(gas.nSpecies() + surf.nSpecies());
        surf.getNetProductionRates(wdot.begin());
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

}
      
