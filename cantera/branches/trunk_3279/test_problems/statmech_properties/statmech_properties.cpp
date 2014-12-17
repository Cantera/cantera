/**
 *  @file statmech
 *       test problem for statistical mechanics in cantera
 */

//  Example
//
// Test case for the statistical mechanics in cantera
//

#include <iostream>
#include <string>
#include <vector>
#include <string>
#include <iomanip>

using namespace std;

/*****************************************************************/
/*****************************************************************/

#include "cantera/transport.h"
#include "cantera/IdealGasMix.h"
#include "cantera/equil/equil.h"

#include "cantera/transport/TransportFactory.h"

using namespace Cantera;

int main(int argc, char** argv)
{

    try {
        int k;
        IdealGasMix g("test_stat.xml");
        int nsp = g.nSpecies();
        double pres = 1.0E5;

        vector_fp Xset(nsp, 0.0);
        Xset[0] =  0.5 ;
        Xset[1] =  0.5;

        g.setState_TPX(1500.0, pres, DATA_PTR(Xset));
        equilibrate(g, "TP", -1);

        vector_fp cp_R(nsp, 0.0);
        g.getCp_R(DATA_PTR(cp_R));

        for (size_t i = 0; i < nsp; i++) {
          cout << cp_R[i] << std::endl;
        }

    } catch (CanteraError) {
        showErrors(cout);
        return 1;
    }

    // Mark it zero!
    return 0;

}
