/////////////////////////////////////////////////////////////
//
//  zero-dimensional kinetics example program
//
//  copyright California Institute of Technology 2006
//
/////////////////////////////////////////////////////////////

#include "cantera/zerodim.h"
#include "cantera/IdealGasMix.h"
#include "example_utils.h"

using namespace Cantera;
using std::cout;
using std::endl;

// Kinetics example. This is written as a function so that one
// driver program can run multiple examples.
// The action taken depends on input parameter job:
//     job = 0:   print a one-line description of the example.
//     job = 1:   print a longer description
//     job = 2:   print description, then run the example.
//

// Note: although this simulation can be done in C++, as shown here,
// it is much easier in Python or Matlab!

int kinetics_example3(int job)
{
    try {

        cout << "Ignition simulation using class IdealGasMix "
             << "with file gri30.cti."
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
        IdealGasMix* gg = new IdealGasMix("gri30.xml", "gri30");
        IdealGasMix& gas = *gg;

        // set the state
        gas.setState_TPX(1001.0, OneAtm, "H2:2.0, O2:1.0, N2:4.0");
        size_t kk = gas.nSpecies();

        // create a reactor
        ConstPressureReactor r;

        // create a reservoir to represent the environment
        Reservoir env;

        // specify the thermodynamic property and kinetics managers
        r.setThermoMgr(gas);
        r.setKineticsMgr(gas);
        env.setThermoMgr(gas);

        // create a flexible, insulating wall between the reactor and the
        // environment
        Wall w;
        w.install(r,env);

        w.setArea(1.0);

        // create a container object to run the simulation
        // and add the reactor to it
        ReactorNet& sim = *(new ReactorNet());
        sim.setVerbose(false);
        sim.addReactor(r);

        double tm;
        double dt = 1.e-5;    // interval at which output is written
        int nsteps = 100;     // number of intervals

        // create a 2D array to hold the output variables,
        // and store the values for the initial state
        Array2D soln(kk+4, 1);
        saveSoln(0, 0.0, gas, soln);

        // main loop
        for (int i = 1; i <= nsteps; i++) {
            tm = i*dt;
            sim.advance(tm);
            saveSoln(tm, gas, soln);
        }

        // make a Tecplot data file and an Excel spreadsheet
        std::string plotTitle = "kinetics example 3: constant-pressure ignition";
        plotSoln("kin3.dat", "TEC", plotTitle, gas, soln);
        plotSoln("kin3.csv", "XL", plotTitle, gas, soln);


        // print final temperature
        cout << " Tfinal = " << r.temperature() << endl;
        cout << "Output files:" << endl
             << "  kin3.csv    (Excel CSV file)" << endl
             << "  kin3.dat    (Tecplot data file)" << endl;

        delete gg;
    } catch (CanteraError& err) {
        // handle exceptions thrown by Cantera
        std::cout << err.what() << std::endl;
        cout << " terminating... " << endl;
        appdelete();
        return -1;
    }
    return 0;
}
