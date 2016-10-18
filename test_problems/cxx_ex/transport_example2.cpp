/////////////////////////////////////////////////////////////
//
//  mixture-averaged transport properties
//
/////////////////////////////////////////////////////////////

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/transport.h"
#include "example_utils.h"
#include "cantera/IdealGasMix.h"

using namespace Cantera;
using std::cout;
using std::endl;

int transport_example2(int job)
{
    try {
        cout << "Multicomponent transport properties." << endl;
        if (job > 0) {
            cout << "Viscosity, thermal conductivity, and thermal diffusion\n"
                 " coefficients at 2 atm for a "
                 << "range of temperatures" << endl;
        }
        if (job <= 1) {
            return 0;
        }

        // create a gas mixture, and set its state

        IdealGasMix gas("gri30.xml", "gri30");
        doublereal temp = 2000.0;
        doublereal pres = 2.0*OneAtm;
        gas.setState_TPX(temp, pres, "H2:1.0, O2:0.5, CH4:0.1, N2:0.2");
        gas.equilibrate("TP");

        // create a transport manager that implements
        // multicomponent transport properties

        Transport* tr = newTransportMgr("Multi",  &gas);
        size_t nsp = gas.nSpecies();

        // create a 2D array to hold the outputs
        int ntemps = 20;
        Array2D output(nsp+3, ntemps);

        // main loop
        for (int i = 0; i < ntemps; i++) {
            temp = 500.0 + 100.0*i;
            gas.setState_TP(temp, pres);
            output(0,i) = temp;
            output(1,i) = tr->viscosity();
            output(2,i) = tr->thermalConductivity();
            tr->getThermalDiffCoeffs(&output(3,i));
        }

        // make a Tecplot data file and an Excel spreadsheet
        std::string plotTitle = "transport example 2: "
                                "multicomponent transport properties";
        plotTransportSoln("tr2.dat", "TEC", plotTitle, gas, output);
        plotTransportSoln("tr2.csv", "XL", plotTitle, gas, output);

        // print final temperature
        cout << "Output files:" << endl
             << "  tr2.csv    (Excel CSV file)" << endl
             << "  tr2.dat    (Tecplot data file)" << endl;
    } catch (CanteraError& err) {
        // handle exceptions thrown by Cantera
        std::cout << err.what() << std::endl;
        cout << " terminating... " << endl;
        appdelete();
        return -1;
    }
    return 0;
}
