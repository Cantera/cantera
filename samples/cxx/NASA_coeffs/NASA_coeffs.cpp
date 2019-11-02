#include "cantera/thermo/IdealGasPhase.h" // defines class IdealGasPhase

#include <iostream>

using namespace Cantera;

// The program is put into a function so that error handling code can
// be conveniently put around the whole thing. See main() below.

void demoprog()
{
    writelog("\n**** Testing modifying NASA polynomial coefficients ****\n");

    auto sol = newSolution("h2o2.yaml", "ohmech");
    auto gas = sol->thermo();
    int nsp = gas->nSpecies();

    int type;
    double c[15];
    double minTemp, maxTemp, refPressure;

    // get a reference to the species thermo property manager
    MultiSpeciesThermo& sp = gas->speciesThermo();

    int n, j;

    // these constants define the location of coefficient "a6" in the
    // cofficient array c. The c array contains Tmid in the first
    // location, followed by the 7 low-temperature coefficients, then
    // the seven high-temperature ones.
    for (n = 0; n < nsp; n++) {
        writelog("\n\n {} (original):", gas->speciesName(n));

        // get the NASA coefficients in array c
        sp.reportParams(n, type, c, minTemp, maxTemp, refPressure);

        // print the unmodified NASA coefficients
        writelog("\n      ");
        for (j = 1; j < 8; j++) {
            writelog("     A{}     ", j);
        }
        writelog("\n  low:");
        for (j = 1; j < 8; j++) {
            writelog(" {:10.4E} ", c[j]);
        }

        writelog("\n high:");
        for (j = 8; j < 15; j++) {
            writelog(" {:10.4E} ", c[j]);
        }
        writelog("\n      ");

        // print the modified NASA coefficients
        writelog("\n\n {} (modified):", gas->speciesName(n).c_str());
        writelog("\n      ");
        for (j = 1; j < 8; j++) {
            writelog("     A{}     ", j);
        }
        writelog("\n  low:");
        for (j = 1; j < 8; j++) {
            writelog(" {:10.4E} ", c[j]);
        }

        writelog("\n high:");
        for (j = 8; j < 15; j++) {
            writelog(" {:10.4E} ", c[j]);
        }
        writelog("\n");
    }
}



int main()
{
    try {
        demoprog();
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
    }
}
