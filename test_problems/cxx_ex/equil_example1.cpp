/////////////////////////////////////////////////////////////
//
//  chemical equilibrium
//
/////////////////////////////////////////////////////////////

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "example_utils.h"

#include "cantera/base/Solution.h"
#include "cantera/thermo/IdealGasPhase.h"

#include <iostream>

using namespace Cantera;
//-------------------------------------------------------------------

// utility functions for plotting

template<class G, class V>
void makeEquilDataLabels(const G& gas, V& names)
{
    size_t nsp = gas->nSpecies();
    names.resize(nsp + 2);
    names[0]  = "Temperature (K)";
    names[1]  = "Pressure (Pa)";
    for (size_t k = 0; k < nsp; k++) {
        names[2+k] = gas->speciesName(k);
    }
}

template<class G, class A>
void plotEquilSoln(const std::string& fname, const std::string& fmt,
                   const std::string& title, const G& gas, const A& soln)
{
    std::vector<std::string> names;
    makeEquilDataLabels(gas, names);
    writePlotFile(fname, fmt, title, names, soln);
}
//-----------------------------------------------------------------


// Equilibrium example. This is written as a function so that one
// driver program can run multiple examples.
// The action taken depends on input parameter job:
//     job = 0:   print a one-line description of the example.
//     job = 1:   print a longer description
//     job = 2:   print description, then run the example.


int equil_example1(int job)
{
    std::cout << "Chemical equilibrium." << std::endl;
    if (job > 0) {
        std::cout << "Equilibrium composition and pressure for a "
                  << "range of temperatures at constant density."
                  << std::endl << std::endl;
    }
    if (job <= 1) {
        return 0;
    }

    // create a gas mixture, and set its state
    auto sol = newSolution("silane.xml", "silane");
    auto gas = sol->thermo();
    size_t nsp = gas->nSpecies();

    int ntemps = 50;   // number of temperatures
    Array2D output(nsp+2, ntemps);

    // main loop
    double temp;
    double thigh = gas->maxTemp();
    double tlow = 500.0;
    double dt = (thigh - tlow)/(ntemps);
    double pres = 0.01*OneAtm;
    for (int i = 0; i < ntemps; i++) {
        temp = tlow + dt*i;
        if (temp > gas->maxTemp()) {
            break;
        }
        gas->setState_TPX(temp, pres, "SIH4:0.01, H2:0.99");

        gas->equilibrate("TP");
        output(0,i) = temp;
        output(1,i) = gas->pressure();
        gas->getMoleFractions(&output(2,i));

    }

    // make a Tecplot data file and an Excel spreadsheet
    std::string plotTitle = "equilibrium example 1: "
                            "chemical equilibrium";
    plotEquilSoln("eq1.dat", "TEC", plotTitle, gas, output);
    plotEquilSoln("eq1.csv", "XL", plotTitle, gas, output);

    std::cout << "Output files:" << std::endl
              << "  eq1.csv    (Excel CSV file)" << std::endl
              << "  eq1.dat    (Tecplot data file)" << std::endl;

    return 0;

}
