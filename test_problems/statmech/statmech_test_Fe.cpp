/**
 *  @file statmech
 *       test problem for statistical mechanics in cantera
 */

//  Example
//
// Test case to check error thrown if using Fe (not supported species)
//

#include "cantera/transport.h"
#include "cantera/IdealGasMix.h"
#include "cantera/equil/equil.h"

using namespace std;
using namespace Cantera;

int main(int argc, char** argv)
{

    try {
        int k;
        IdealGasMix g("test_stat_Fe.xml");
        int nsp = g.nSpecies();
        double pres = 1.0E5;

        vector_fp Xset(nsp, 0.0);
        Xset[0] =  0.5 ;
        Xset[1] =  0.5;

        g.setState_TPX(1500.0, pres, DATA_PTR(Xset));
        equilibrate(g, "TP", -1);
        g.report();

        vector_fp cp_R(nsp, 0.0);
        g.getCp_R(DATA_PTR(cp_R));

        for (int i=0; i<nsp; i++) {
            std::cout << cp_R[i] << std::endl;
        }

        // error check
        if (cp_R[0] != 0) {
            std::cout << "Error for monotomic Species!\n";
            return 1;
        }
    } catch (CanteraError) {
        // need to get error here because of loading Fe in input file, when
        // no Fe exists in the species information table, in statmech.cpp
        showErrors(cout);
        return 0;
    }

    // Mark it zero!
    return 1;

}
