// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/IdealGasPhase.h"

#include <iostream>

using namespace std;
using namespace Cantera;

int main(int argc, char** argv)
{
#if defined(_MSC_VER) && _MSC_VER < 1900
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    try {
        auto sol = newSolution("bad_air.xml", "air");
        auto gas = sol->thermo();
        double pres = 1.0E5;
        gas->setState_TPX(1000.1, pres, "O2:0.4, N2:0.6");
        gas->equilibrate("TP", "auto");
        double enth = gas->enthalpy_mass();
        printf(" enth = %g\n", enth);
        enth -= 2.0E2;
        printf("attempted equil at (H,P) = %10.5g, %10.5g\n", enth, pres);
        gas->setState_HP(enth, pres);
        gas->equilibrate("HP", "auto");
        enth = gas->enthalpy_mass();
        printf(" enth = %g\n", enth);


        return 0;
    } catch (CanteraError& err) {
        std::cerr << err.what() << std::endl;
        cerr << "program terminating." << endl;
        return -1;
    }
}
