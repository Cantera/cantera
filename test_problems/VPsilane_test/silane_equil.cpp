// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/IdealGasMix.h"
#include "cantera/thermo/IdealSolnGasVPSS.h"
#include "cantera/thermo/ThermoFactory.h"

using namespace std;
using namespace Cantera;

int main(int argc, char** argv)
{
#if defined(_MSC_VER) && _MSC_VER < 1900
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    try {
        IdealSolnGasVPSS gg("silane.xml", "silane");
        ThermoPhase* g = &gg;
        cout.precision(4);
        g->setState_TPX(1500.0, 100.0, "SIH4:0.01, H2:0.99");
        g->equilibrate("TP");
        cout << g->report(true, 0.0) << endl;
        return 0;
    } catch (CanteraError& err) {
        std::cerr << err.what() << std::endl;
        cerr << "program terminating." << endl;
        return -1;
    }
}
