// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

/*
 *  Sample program that solves an implicit problem for surface
 *  site fractions.
 */

#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/InterfaceKinetics.h"

#include <iostream>

using namespace Cantera;
using namespace std;

int main()
{
    try {
        auto sol = newSolution("gri30.yaml", "gri30");
        auto gas = sol->thermo();
        gas->setState_TPX(1200.0, OneAtm,
                          "H2:2, O2:1, OH:0.01, H:0.01, O:0.01");

        auto surf = newSolution("surface.xml", "surface", "None", {sol});
        auto surf_ph = std::dynamic_pointer_cast<SurfPhase>(surf->thermo());
        auto surf_kin = std::dynamic_pointer_cast<InterfaceKinetics>(surf->kinetics());

        vector_fp cov { 0.8, 0.2 };
        cout.precision(4);
        surf_ph->setCoverages(cov.data());
        vector_fp wdot(gas->nSpecies() + surf_ph->nSpecies());
        surf_kin->getNetProductionRates(wdot.data());
        for (size_t k = 0; k < gas->nSpecies(); k++) {
            cout << gas->speciesName(k) << "  " << wdot[k] << endl;
        }

        for (size_t k = 0; k < surf_ph->nSpecies(); k++) {
            cout << surf_ph->speciesName(k) << "   "
                 << wdot[k+gas->nSpecies()] << endl;
        }
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
    }
    return 0;
}
