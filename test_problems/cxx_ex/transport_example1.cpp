/////////////////////////////////////////////////////////////
//
//  mixture-averaged transport properties
//
/////////////////////////////////////////////////////////////

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/transport.h"
#include "example_utils.h"
#include "cantera/thermo/IdealGasPhase.h"

using namespace Cantera;
using std::cout;
using std::endl;

int transport_example1(int job)
{
    try {
        cout << "Mixture-averaged transport properties." << endl;
        if (job > 0) {
            cout << "Viscosity, thermal conductivity, and mixture-averaged\n"
                 << "diffusion coefficients at 2 atm for a "
                 << "range of temperatures" << endl;
        }
        if (job <= 1) {
            return 0;
        }

        // create a gas mixture, and set its state

        auto sol = newSolution("gri30.yaml", "gri30", "Mix");
        auto gas = sol->thermo();
        double temp = 500.0;
        double pres = 2.0*OneAtm;
        gas->setState_TPX(temp, pres, "H2:1.0, CH4:0.1");

        // create a transport manager that implements
        // mixture-averaged transport properties

        auto tr = sol->transport();

        size_t nsp = gas->nSpecies();


        // create a 2D array to hold the outputs
        int ntemps = 20;
        Array2D output(nsp+3, ntemps);

        // main loop
        for (int i = 0; i < ntemps; i++) {
            temp = 500.0 + 100.0*i;
            gas->setState_TP(temp, pres);
            output(0,i) = temp;
            output(1,i) = tr->viscosity();
            output(2,i) = tr->thermalConductivity();
            tr->getMixDiffCoeffs(&output(3,i));
        }

        // make a Tecplot data file and an Excel spreadsheet
        std::string plotTitle = "transport example 1: "
                                "mixture-averaged transport properties";
        plotTransportSoln("tr1.dat", "TEC", plotTitle, *(sol->thermo()), output);
        plotTransportSoln("tr1.csv", "XL", plotTitle, *(sol->thermo()), output);

        // print final temperature
        cout << "Output files:" << endl
             << "  tr1.csv    (Excel CSV file)" << endl
             << "  tr1.dat    (Tecplot data file)" << endl;
    } catch (CanteraError& err) {
        // handle exceptions thrown by Cantera
        std::cout << err.what() << std::endl;
        cout << " terminating... " << endl;
        appdelete();
        return -1;
    }
    return 0;
}
