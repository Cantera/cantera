/////////////////////////////////////////////////////////////
//
//  zero-dimensional kinetics example program
//
/////////////////////////////////////////////////////////////

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zerodim.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "example_utils.h"
using namespace Cantera;
using namespace std;

// Kinetics example. This is written as a function so that one
// driver program can run multiple examples.
// The action taken depends on input parameter job:
//     job = 0:   print a one-line description of the example.
//     job = 1:   print a longer description
//     job = 2:   print description, then run the example.

// Note: although this simulation can be done in C++, as shown here,
// it is much easier in Python or Matlab!

int kinetics_example1(int job)
{
    try {

        cout << "Ignition simulation using class IdealGasPhase "
             << "with file gri30.yaml."
             << endl;

        if (job >= 1) {
            cout << "Constant-pressure ignition of a "
                 << "hydrogen/oxygen/nitrogen"
                 " mixture \nbeginning at T = 1001 K and P = 1 atm." << endl;
        }
        if (job < 2) {
            return 0;
        }

        // create an ideal gas mixture that corresponds to GRI-Mech
        // 3.0
        auto sol = newSolution("gri30.yaml", "gri30", "None");
        auto gas = sol->thermo();

        // set the state
        gas->setState_TPX(1001.0, OneAtm, "H2:2.0, O2:1.0, N2:4.0");
        size_t kk = gas->nSpecies();

        // create a reactor
        Reactor r;
        r.insert(sol);

        // create a reservoir to represent the environment
        Reservoir env;
        env.insert(sol);

        // create a flexible, insulating wall between the reactor and the
        // environment
        Wall w;
        w.install(r,env);

        // set the "Vdot coefficient" to a large value, in order to
        // approach the constant-pressure limit; see the documentation
        // for class Reactor
        w.setExpansionRateCoeff(1.e9);
        w.setArea(1.0);

        // create a container object to run the simulation
        // and add the reactor to it
        ReactorNet sim;
        sim.setVerbose(false);
        sim.addReactor(r);

        double tm;
        double dt = 1.e-5;    // interval at which output is written
        int nsteps = 100;     // number of intervals

        // create a 2D array to hold the output variables,
        // and store the values for the initial state
        Array2D soln(kk+4, 1);
        saveSoln(0, 0.0, *(sol->thermo()), soln);

        // main loop
        for (int i = 1; i <= nsteps; i++) {
            tm = i*dt;
            sim.advance(tm);
            saveSoln(tm, *(sol->thermo()), soln);
        }

        // make a Tecplot data file and an Excel spreadsheet
        string plotTitle = "kinetics example 1: constant-pressure ignition";
        plotSoln("kin1.dat", "TEC", plotTitle, *(sol->thermo()), soln);
        plotSoln("kin1.csv", "XL", plotTitle, *(sol->thermo()), soln);

        // print final temperature
        cout << " Tfinal = " << r.temperature() << endl;
        cout << "Output files:" << endl
             << "  kin1.csv    (Excel CSV file)" << endl
             << "  kin1.dat    (Tecplot data file)" << endl;
    } catch (CanteraError& err) {
        // handle exceptions thrown by Cantera
        std::cout << err.what() << std::endl;
        cout << " terminating... " << endl;
        appdelete();
        return -1;
    }
    return 0;
}
