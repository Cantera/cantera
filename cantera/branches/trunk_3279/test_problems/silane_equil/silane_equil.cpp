/*
 *  Copyright 2002 California Institute of Technology
 */

#include "cantera/IdealGasMix.h"
#include "cantera/equilibrium.h"

using namespace std;
using namespace Cantera;

int main(int argc, char** argv)
{
#ifdef _MSC_VER
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    try {
        IdealGasMix g("silane.xml", "silane");
        g.setState_TPX(1500.0, 100.0, "SIH4:0.01, H2:0.99");
        equilibrate(g, "TP");
        //cout << g;
        return 0;
    } catch (CanteraError& err) {
        std::cerr << err.what() << std::endl;
        cerr << "program terminating." << endl;
        return -1;
    }
}
