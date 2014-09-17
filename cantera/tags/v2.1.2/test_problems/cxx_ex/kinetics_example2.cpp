/////////////////////////////////////////////////////////////
//
//  zero-dimensional kinetics example program
//
//  copyright California Institute of Technology 2002
//
/////////////////////////////////////////////////////////////

#include "cantera/GRI30.h"
#include "cantera/zerodim.h"
#include "example_utils.h"

using namespace Cantera;
/**
 * Same as kinetics_example1, except that it uses class GRI30 instead
 * of class IdealGasMix.
 */

// Note: although this simulation can be done in C++, as shown here,
// it is much easier in Python or Matlab!

int kinetics_example2(int job)
{
    suppress_deprecation_warnings();
    try {
        std::cout << "Ignition simulation using class GRI30." << std::endl;

        if (job >= 1) {
            std::cout <<
                      "Constant-pressure ignition of a hydrogen/oxygen/nitrogen"
                      " mixture \nbeginning at T = 1001 K and P = 1 atm." << std::endl;
        }
        if (job < 2) {
            return 0;
        }

        // create a GRI30 object
        GRI30 gas;
        gas.setState_TPX(1001.0, OneAtm, "H2:2.0, O2:1.0, N2:4.0");
        size_t kk = gas.nSpecies();

        // create a reactor
        Reactor r;

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

        // set the "Vdot coefficient" to a large value, in order to
        // approach the constant-pressure limit; see the documentation
        // for class Reactor
        w.setExpansionRateCoeff(1.e9);
        w.setArea(1.0);

        // create a container object to run the simulation
        // and add the reactor to it
        ReactorNet* sim_ptr = new ReactorNet();
        ReactorNet& sim = *sim_ptr;
        sim.setVerbose(false);
        sim.addReactor(&r);

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
        std::string plotTitle = "kinetics example 2: constant-pressure ignition";
        plotSoln("kin2.dat", "TEC", plotTitle, gas, soln);
        plotSoln("kin2.csv", "XL", plotTitle, gas, soln);


        // print final temperature
        std::cout << " Tfinal = " << r.temperature() << std::endl;
        std::cout << "Output files:" << std::endl
                  << "  kin2.csv    (Excel CSV file)" << std::endl
                  << "  kin2.dat    (Tecplot data file)" << std::endl;
    } catch (CanteraError& err) {
        // handle exceptions thrown by Cantera
        std::cout << err.what() << std::endl;
        std::cout << " terminating... " << std::endl;
        appdelete();
        return -1;
    }
    return 0;
}
