/////////////////////////////////////////////////////////////
//
//  reaction path diagrams
//
/////////////////////////////////////////////////////////////

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zerodim.h"
#include "example_utils.h"
#include "cantera/kinetics/ReactionPath.h"
#include "cantera/thermo/IdealGasPhase.h"

using namespace Cantera;
using std::cout;
using std::endl;

void writeRxnPathDiagram(double time, ReactionPathBuilder& b,
                         Kinetics& kin, std::ostream& logfile, std::ostream& outfile)
{
    // create a new empty diagram
    ReactionPathDiagram d;

    // show the details of which reactions contribute to the flux
    d.show_details = false;

    // set the threshold for the minimum flux relative value that will
    // be plotted
    d.threshold = 0.001;

    // color for bold lines
    d.bold_color = "orange";

    // color for normal-weight lines
    d.normal_color = "steelblue";

    // color for dashed lines
    d.dashed_color = "gray";

    // options for the 'dot' program
    d.dot_options = "center=1;size=\"6,9\";ratio=auto";

    // minimum relative flux for bold lines
    d.bold_min = 0.0;

    // maximum relative flux for dashed lines
    d.dashed_max = 0.01;

    // minimum relative flux for labels
    d.label_min = 0.01;

    // autoscale
    d.scale = -1;

    // set to either NetFlow or OneWayFlow
    d.flow_type = NetFlow; //OneWayFlow;

    // arrow width. If < 0, then scale with flux value
    d.arrow_width = -2.0;

    // title
    d.title = fmt::format("time = {} (s)", time);

    // build the diagram following elemental nitrogen
    b.build(kin, "N", logfile, d);

    // write an input file for 'dot'
    d.exportToDot(outfile);

}


int rxnpath_example1(int job)
{
    try {

        cout << "Reaction path diagram movies with file gri30.cti." << endl;
        if (job >= 1) {
            cout << "Generate reaction path diagrams following nitrogen\n"
                 << "as a function of time for constant-pressure ignition of a\n"
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
        gas->setState_TPX(1001.0, OneAtm, "H2:2.0, O2:1.0, N2:4.0");

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

        double tm;
        double dt = 1.e-5;    // interval at which output is written
        int nsteps = 100;     // number of intervals

        // create a container object to run the simulation
        // and add the reactor to it
        ReactorNet& sim = *(new ReactorNet());
        sim.addReactor(r);

        // create a reaction path diagram builder
        ReactionPathBuilder b;
        std::ofstream rplog("rp1.log");   // log file
        std::ofstream rplot("rp1.dot");   // output file
        b.init(rplog, *(sol->kinetics()));   // initialize

        // main loop
        for (int i = 1; i <= nsteps; i++) {
            tm = i*dt;
            sim.advance(tm);
            writeRxnPathDiagram(tm, b, *(sol->kinetics()), rplog, rplot);
        }

        // print final temperature
        cout << "Output files:" << endl
             << "  rp1.log    (log file)" << endl
             << "  rp1.dot    (input file for dot)" << endl;
        cout << "To generate the diagrams in Postscript, execute the command" << endl << endl
             << "dot -Tps rp1.dot > rp1.ps" << endl << endl
             << "Get dot for Windows here: http://blue.caltech.edu/dot.exe" << endl;
    } catch (CanteraError& err) {
        // handle exceptions thrown by Cantera
        std::cout << err.what() << std::endl;
        cout << " terminating... " << endl;
        appdelete();
        return -1;
    }
    return 0;
}
