/////////////////////////////////////////////////////////////
//
//  chemical equilibrium
//
//  $Author: hkmoffa $
//  $Revision: 1.15 $
//  $Date: 2008/02/16 21:33:37 $
//
//  copyright California Institute of Technology 2002
//
/////////////////////////////////////////////////////////////

// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include <cantera/Cantera.h>
#include <time.h>
#include "example_utils.h"
#include <cantera/equilibrium.h>

#include <cantera/IdealGasMix.h>

using namespace Cantera;
using namespace Cantera_CXX;
//-------------------------------------------------------------------

// utility functions for plotting

template<class G, class V>
void makeEquilDataLabels(const G& gas, V& names) {
    int nsp = gas.nSpecies();
    names.resize(nsp + 2);
    names[0]  = "Temperature (K)";
    names[1]  = "Pressure (Pa)";
    int k;
    for (k = 0; k < nsp; k++) names[2+k] = gas.speciesName(k);
}

template<class G, class A>
void plotEquilSoln(string fname, string fmt, string title, const G& gas, 
    const A& soln) {
    vector<string> names;
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


int equil_example1(int job) {

    cout << "Chemical equilibrium." << endl;
    if (job > 0) {
        cout << "Equilibrium composition and pressure for a "
             << "range of temperatures at constant density." << endl;
    }
    if (job <= 1) return 0; 

    // header
    writeCanteraHeader(cout);

    // create a gas mixture, and set its state

    //IdealGasMix gas("silane.cti", "silane");
    IdealGasMix gas("silane.xml", "silane");
    int nsp = gas.nSpecies();

    int ntemps = 50;   // number of temperatures
    Array2D output(nsp+2, ntemps);
        
    // main loop
    doublereal temp;
    doublereal thigh = gas.maxTemp();
    doublereal tlow = 500.0;
    doublereal dt = (thigh - tlow)/(ntemps);
    doublereal pres = 0.01*OneAtm;
    clock_t t0 = clock();
    for (int i = 0; i < ntemps; i++) {
        temp = tlow + dt*i;
        if (temp > gas.maxTemp()) break;
        gas.setState_TPX(temp, pres, "SIH4:0.01, H2:0.99");

        //        equilibrate(gas,"TP",1,1.0e-9,1000,100,15);
        equilibrate(gas,"TP");
        output(0,i) = temp;
        output(1,i) = gas.pressure();
        gas.getMoleFractions(&output(2,i));

    }
    clock_t t1 = clock();

    // make a Tecplot data file and an Excel spreadsheet
    string plotTitle = "equilibrium example 1: "
                       "chemical equilibrium";
    plotEquilSoln("eq1.dat", "TEC", plotTitle, gas, output);
    plotEquilSoln("eq1.csv", "XL", plotTitle, gas, output);

    // print timing data
    doublereal tmm = 1.0*(t1 - t0)/CLOCKS_PER_SEC;
    cout << " time = " << tmm << endl << endl;

    cout << "Output files:" << endl
         << "  eq1.csv    (Excel CSV file)" << endl
         << "  eq1.dat    (Tecplot data file)" << endl;

    return 0;

}
