/////////////////////////////////////////////////////////////
//
//  zero-dimensional kinetics example program
//
//  copyright California Institute of Technology 2002
//
/////////////////////////////////////////////////////////////

#include "cantera/zerodim.h"
#include "cantera/IdealGasMix.h"
#include "cantera/numerics.h"

#include <time.h>
#include "example_utils.h"

using namespace Cantera;
using std::cout;
using std::endl;

int kinetics1(int np, void* p)
{

    cout << "Constant-pressure ignition of a "
         << "hydrogen/oxygen/nitrogen"
         " mixture \nbeginning at T = 1001 K and P = 1 atm." << endl;

    // create an ideal gas mixture that corresponds to GRI-Mech
    // 3.0
    IdealGasMix gas("gri30.cti", "gri30");

    // set the state
    gas.setState_TPX(1001.0, OneAtm, "H2:2.0, O2:1.0, N2:4.0");
    int nsp = gas.nSpecies();

    // create a reactor
    Reactor r;

    // create a reservoir to represent the environment
    Reservoir env;

    // 'insert' the gas into the reactor and environment.  Note
    // that it is ok to insert the same gas object into multiple
    // reactors or reservoirs. All this means is that this object
    // will be used to evaluate thermodynamic or kinetic
    // quantities needed.
    r.insert(gas);
    env.insert(gas);

    //r.addHomogenRxnSens(0);

    // create a wall between the reactor and the environment
    Wall w;
    w.install(r,env);

    // The wall "expansion rate coefficient" controls how fast it
    // moves in response to a pressure difference. Set it to a
    // large value to approach the constant-pressure limit, so
    // that the wall moves to counteract even small pressure
    // differences
    w.setExpansionRateCoeff(1.e9);

    // set the wall to have unit area (arbitrary)
    w.setArea(1.0);

    double tm;
    double dt = 1.e-5;    // interval at which output is written
    int nsteps = 100;     // number of intervals

    // create a 2D array to hold the output variables,
    // and store the values for the initial state
    Array2D soln(nsp+4, 1);
    saveSoln(0, 0.0, gas, soln);

    // create a container object to run the simulation
    // and add the reactor to it
    ReactorNet sim;
    sim.addReactor(&r, false);

    // main loop
    clock_t t0 = clock();        // save start time
    for (int i = 1; i <= nsteps; i++) {
        tm = i*dt;
        sim.advance(tm);
        cout << "time = " << tm << " s" << endl;
        saveSoln(tm, gas, soln);
    }
    clock_t t1 = clock();        // save end time


    // make a Tecplot data file and an Excel spreadsheet
    std::string plotTitle = "kinetics example 1: constant-pressure ignition";
    plotSoln("kin1.dat", "TEC", plotTitle, gas, soln);
    plotSoln("kin1.csv", "XL", plotTitle, gas, soln);


    // print final temperature and timing data
    doublereal tmm = 1.0*(t1 - t0)/CLOCKS_PER_SEC;
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


#ifndef CXX_DEMO

int main()
{

    try {
        int retn = kinetics1(0, 0);
        appdelete();
        return retn;
    }
    // handle exceptions thrown by Cantera
    catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        cout << " terminating... " << endl;
        appdelete();
        return -1;
    }
}

#endif
