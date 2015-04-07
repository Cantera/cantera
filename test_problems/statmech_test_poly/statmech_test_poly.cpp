/**
 *  @file statmech
 *       test problem for statistical mechanics in cantera
 */

//  Example
//
// Test case for the statistical mechanics in cantera
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
        IdealGasMix g("test_stat_err.xml");
        int nsp = g.nSpecies();
        double pres = 1.0E5;

        vector_fp Xset(nsp, 0.0);
        Xset[0] =  0.5 ;
        Xset[1] =  0.5;

        g.setState_TPX(1500.0, pres, DATA_PTR(Xset));
        equilibrate(g, "TP", -1);

        vector_fp cp_R(nsp, 0.0);
        g.getCp_R(DATA_PTR(cp_R));

    } catch (CanteraError) {
        // we wanted to catch an error here for incorrectly trying to use poly methods
        // for the statmech species data, so exit successfully
        showErrors(cout);

        // Mark it zero!
        return 0;
    }

    cout << "ERROR" << std::endl;
    // something is wrong: we were suppose to catch an error here: paradox!
    return 1;

}
