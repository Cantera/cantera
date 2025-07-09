/*
 * Freely-propagating flame
 * ========================
 *
 * C++ demo program to compute flame speeds using GRI-Mech.
 *
 * Usage:
 *
 *     flamespeed [equivalence_ratio] [refine_grid] [loglevel]
 *
 * Requires: cantera >= 3.1
 *
 * .. tags:: C++, combustion, 1D flow, premixed flame, flame speed, saving output
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/onedim.h"
#include "cantera/oneD/DomainFactory.h"
#include "cantera/base/stringUtils.h"
#include <fstream>

using namespace Cantera;
using fmt::print;

int flamespeed(double phi, std::string refine_grid, int loglevel)
{
    try {
        auto sol = newSolution("gri30.yaml", "gri30", "mixture-averaged");
        auto gas = sol->thermo();
        double temp = 300.0; // K
        double pressure = 1.0*OneAtm; //atm
        double uin = 0.3; //m/sec

        size_t nsp = gas->nSpecies();
        vector<double> x(nsp, 0.0);

        gas->setEquivalenceRatio(phi, "CH4", "O2:0.21,N2:0.79");
        gas->setState_TP(temp, pressure);
        gas->getMoleFractions(x.data());

        double rho_in = gas->density();

        vector<double> yin(nsp);
        gas->getMassFractions(&yin[0]);

        gas->equilibrate("HP");
        vector<double> yout(nsp);
        gas->getMassFractions(&yout[0]);
        double rho_out = gas->density();
        double Tad = gas->temperature();
        print("phi = {}, Tad = {}\n", phi, Tad);

        //=============  build each domain ========================

        //-------- step 1: create the flow -------------

        auto flow = newDomain<Flow1D>("gas-flow", sol, "flow");
        flow->setFreeFlow();

        // create an initial grid
        int nz = 6;
        double lz = 0.1;
        vector<double> z(nz);
        double dz = lz/((double)(nz-1));
        for (int iz = 0; iz < nz; iz++) {
            z[iz] = ((double)iz)*dz;
        }

        flow->setupGrid(nz, &z[0]);

        //------- step 2: create the inlet  -----------------------

        auto inlet = newDomain<Inlet1D>("inlet", sol);

        inlet->setMoleFractions(x.data());
        double mdot = uin * rho_in;
        inlet->setMdot(mdot);
        inlet->setTemperature(temp);

        //------- step 3: create the outlet  ---------------------

        auto outlet = newDomain<Outlet1D>("outlet", sol);

        //=================== create the container and insert the domains =====

        vector<shared_ptr<Domain1D>> domains { inlet, flow, outlet };
        Sim1D flame(domains);

        //----------- Supply initial guess----------------------

        vector<double> locs{0.0, 0.3, 0.7, 1.0};
        vector<double> value;

        double uout = inlet->mdot()/rho_out;
        value = {uin, uin, uout, uout};
        flame.setInitialGuess("velocity",locs,value);
        value = {temp, temp, Tad, Tad};
        flame.setInitialGuess("T",locs,value);

        for (size_t i=0; i<nsp; i++) {
            value = {yin[i], yin[i], yout[i], yout[i]};
            flame.setInitialGuess(gas->speciesName(i),locs,value);
        }

        inlet->setMoleFractions(x.data());
        inlet->setMdot(mdot);
        inlet->setTemperature(temp);

        flame.show();

        int flowdomain = 1;
        double ratio = 10.0;
        double slope = 0.08;
        double curve = 0.1;

        flame.setRefineCriteria(flowdomain,ratio,slope,curve);

        // Save initial guess to container file

        // Solution is saved in HDF5 or YAML file format
        string fileName;
        if (usesHDF5()) {
            // Cantera is compiled with native HDF5 support
            fileName = "flamespeed.h5";
        } else {
            fileName = "flamespeed.yaml";
        }
        flame.save(fileName, "initial-guess", "Initial guess", true);

        // Solve freely propagating flame

        // Linearly interpolate to find location where this temperature would
        // exist. The temperature at this location will then be fixed for
        // remainder of calculation.
        flame.setFixedTemperature(0.5 * (temp + Tad));
        flow->solveEnergyEqn();

        flame.solve(loglevel,refine_grid);
        double flameSpeed_mix = flame.value(flowdomain,
                                            flow->componentIndex("velocity"),0);
        print("Flame speed with mixture-averaged transport: {} m/s\n", flameSpeed_mix);
        flame.save(fileName, "mix", "Solution with mixture-averaged transport", true);

        // now switch to multicomponent transport
        flow->setTransportModel("multicomponent");
        flame.solve(loglevel, refine_grid);
        double flameSpeed_multi = flame.value(flowdomain,
                                              flow->componentIndex("velocity"),0);
        print("Flame speed with multicomponent transport: {} m/s\n", flameSpeed_multi);
        flame.save(fileName, "multi", "Solution with multicomponent transport", true);

        // now enable Soret diffusion
        flow->enableSoret(true);
        flame.solve(loglevel, refine_grid);
        double flameSpeed_full = flame.value(flowdomain,
                                             flow->componentIndex("velocity"),0);
        print("Flame speed with multicomponent transport + Soret: {} m/s\n",
              flameSpeed_full);
        flame.save(fileName, "soret",
                   "Solution with mixture-averaged transport and Soret", true);

        vector<double> zvec,Tvec,COvec,CO2vec,Uvec;

        print("\n{:9s}\t{:8s}\t{:5s}\t{:7s}\n",
              "z (m)", "T (K)", "U (m/s)", "Y(CO)");
        for (size_t n = 0; n < flow->nPoints(); n++) {
            Tvec.push_back(flame.value(flowdomain,flow->componentIndex("T"),n));
            COvec.push_back(flame.value(flowdomain,
                                        flow->componentIndex("CO"),n));
            CO2vec.push_back(flame.value(flowdomain,
                                         flow->componentIndex("CO2"),n));
            Uvec.push_back(flame.value(flowdomain,
                                       flow->componentIndex("velocity"),n));
            zvec.push_back(flow->z(n));
            print("{:9.6f}\t{:8.3f}\t{:5.3f}\t{:7.5f}\n",
                  flow->z(n), Tvec[n], Uvec[n], COvec[n]);
        }

        print("\nAdiabatic flame temperature from equilibrium is: {}\n", Tad);
        print("Flame speed for phi={} is {} m/s.\n", phi, Uvec[0]);

        std::ofstream outfile("flamespeed.csv", std::ios::trunc);
        outfile << "  Grid,   Temperature,   Uvec,   CO,    CO2\n";
        for (size_t n = 0; n < flow->nPoints(); n++) {
            print(outfile, " {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}\n",
                  flow->z(n), Tvec[n], Uvec[n], COvec[n], CO2vec[n]);
        }
    } catch (CanteraError& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << "program terminating." << std::endl;
        return -1;
    }
    return 0;
}

int main(int argc, char** argv)
{
    double phi;
    int loglevel = 1;
    std::string refine_grid = "refine";
    if (argc >= 2) {
        phi = fpValue(argv[1]);
    } else {
        print("Enter phi: ");
        std::cin >> phi;
    }
    if (argc >= 3) {
        refine_grid = std::string(argv[2]);
    }
    if (argc >= 4) {
        loglevel = std::stoi(argv[3]);
    }
    return flamespeed(phi, refine_grid, loglevel);
}
