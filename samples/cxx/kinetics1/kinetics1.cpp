/////////////////////////////////////////////////////////////
//
// zero-dimensional kinetics example program
//
/////////////////////////////////////////////////////////////

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zerodim.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "example_utils.h"

using namespace Cantera;
using std::cout;
using std::endl;

int kinetics1(int np, void* p)
{
    cout << "Constant-pressure ignition of a "
         << "hydrogen/oxygen/nitrogen"
         " mixture \nbeginning at T = 1001 K and P = 1 atm." << endl;

    // create an ideal gas mixture that corresponds to GRI-Mech 3.0
    auto sol = newSolution("gri30.yaml", "gri30", "None");
    auto gas = sol->thermo();

    // set the state
    gas->setState_TPX(1001.0, OneAtm, "H2:2.0, O2:1.0, N2:4.0");
    int nsp = gas->nSpecies();

    // create a reactor
    IdealGasConstPressureReactor r;

    // 'insert' the gas into the reactor and environment.  Note
    // that it is ok to insert the same gas object into multiple
    // reactors or reservoirs. All this means is that this object
    // will be used to evaluate thermodynamic or kinetic
    // quantities needed.
    r.insert(sol);

    double dt = 1.e-5; // interval at which output is written
    int nsteps = 100; // number of intervals

    // create a 2D array to hold the output variables,
    // and store the values for the initial state
    Array2D soln(nsp+4, 1);
    saveSoln(0, 0.0, *(sol->thermo()), soln);

    // create a container object to run the simulation
    // and add the reactor to it
    ReactorNet sim;
    sim.addReactor(r);

    // main loop
    clock_t t0 = clock(); // save start time
    for (int i = 1; i <= nsteps; i++) {
        double tm = i*dt;
        sim.advance(tm);
        cout << "time = " << tm << " s" << endl;
        saveSoln(tm, *(sol->thermo()), soln);
    }
    clock_t t1 = clock(); // save end time


    // make a Tecplot data file and an Excel spreadsheet
    std::string plotTitle = "kinetics example 1: constant-pressure ignition";
    plotSoln("kin1.dat", "TEC", plotTitle, *(sol->thermo()), soln);
    plotSoln("kin1.csv", "XL", plotTitle, *(sol->thermo()), soln);


    // print final temperature and timing data
    double tmm = 1.0*(t1 - t0)/CLOCKS_PER_SEC;
    cout << " Tfinal = " << r.temperature() << endl;
    cout << " time = " << tmm << endl;
    cout << " number of residual function evaluations = "
         << sim.integrator().nEvals() << endl;
    cout << " time per evaluation = " << tmm/sim.integrator().nEvals()
         << endl << endl;
    cout << "Output files:" << endl
         << "  kin1.csv    (Excel CSV file)" << endl
         << "  kin1.dat    (Tecplot data file)" << endl;

    return 0;
}


int main()
{
    try {
        int retn = kinetics1(0, 0);
        appdelete();
        return retn;
    } catch (CanteraError& err) {
        // handle exceptions thrown by Cantera
        std::cout << err.what() << std::endl;
        cout << " terminating... " << endl;
        appdelete();
        return -1;
    }
}
