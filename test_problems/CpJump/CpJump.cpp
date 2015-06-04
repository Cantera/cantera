/*
 *  Copyright 2002 California Institute of Technology
 */

#include "cantera/IdealGasMix.h"

using namespace std;
using namespace Cantera;

int main(int argc, char** argv)
{
#ifdef _MSC_VER
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    try {
        IdealGasMix g("bad_air.xml", "air");
        double pres = 1.0E5;
        g.setState_TPX(1000.1, pres, "O2:0.4, N2:0.6");
        g.equilibrate("TP", "auto");
        //cout << g;
        double enth = g.enthalpy_mass();
        printf(" enth = %g\n", enth);
        enth -= 2.0E2;
        printf("attempted equil at (H,P) = %10.5g, %10.5g\n", enth, pres);
        g.setState_HP(enth, pres);
        g.equilibrate("HP", "auto");
        enth = g.enthalpy_mass();
        printf(" enth = %g\n", enth);
        //cout << g;


        return 0;
    } catch (CanteraError& err) {
        std::cerr << err.what() << std::endl;
        cerr << "program terminating." << endl;
        return -1;
    }
}
