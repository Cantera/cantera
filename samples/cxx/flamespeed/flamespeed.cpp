/*!
 * @file flamespeed.cpp
 * C++ demo program to compute flame speeds using GRI-Mech.
 */

#include "cantera/oneD/Sim1D.h"
#include "cantera/oneD/Boundary1D.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/transport.h"
#include <fstream>

using namespace Cantera;
using fmt::print;

int flamespeed(double phi)
{
    try {
        auto sol = newSolution("gri30.yaml", "gri30", "None");
        auto gas = sol->thermo();
        double temp = 300.0; // K
        double pressure = 1.0*OneAtm; //atm
        double uin = 0.3; //m/sec

        size_t nsp = gas->nSpecies();
        vector_fp x(nsp, 0.0);

        double C_atoms = 1.0;
        double H_atoms = 4.0;
        double ax = C_atoms + H_atoms / 4.0;
        double fa_stoic = 1.0 / (4.76 * ax);
        x[gas->speciesIndex("CH4")] = 1.0;
        x[gas->speciesIndex("O2")] = 0.21 / phi / fa_stoic;
        x[gas->speciesIndex("N2")] = 0.79 / phi/ fa_stoic;

        gas->setState_TPX(temp, pressure, x.data());
        double rho_in = gas->density();

        vector_fp yin(nsp);
        gas->getMassFractions(&yin[0]);

        gas->equilibrate("HP");
        vector_fp yout(nsp);
        gas->getMassFractions(&yout[0]);
        double rho_out = gas->density();
        double Tad = gas->temperature();
        print("phi = {}, Tad = {}\n", phi, Tad);

        //=============  build each domain ========================


        //-------- step 1: create the flow -------------

        StFlow flow(gas);
        flow.setFreeFlow();

        // create an initial grid
        int nz = 6;
        double lz = 0.1;
        vector_fp z(nz);
        double dz = lz/((double)(nz-1));
        for (int iz = 0; iz < nz; iz++) {
            z[iz] = ((double)iz)*dz;
        }

        flow.setupGrid(nz, &z[0]);

        // specify the objects to use to compute kinetic rates and
        // transport properties

        std::unique_ptr<Transport> trmix(newTransportMgr("Mix", sol->thermo().get()));
        std::unique_ptr<Transport> trmulti(newTransportMgr("Multi", sol->thermo().get()));

        flow.setTransport(*trmix);
        flow.setKinetics(*sol->kinetics());
        flow.setPressure(pressure);

        //------- step 2: create the inlet  -----------------------

        Inlet1D inlet;

        inlet.setMoleFractions(x.data());
        double mdot=uin*rho_in;
        inlet.setMdot(mdot);
        inlet.setTemperature(temp);


        //------- step 3: create the outlet  ---------------------

        Outlet1D outlet;

        //=================== create the container and insert the domains =====

        std::vector<Domain1D*> domains { &inlet, &flow, &outlet };
        Sim1D flame(domains);

        //----------- Supply initial guess----------------------

        vector_fp locs{0.0, 0.3, 0.7, 1.0};
        vector_fp value;

        double uout = inlet.mdot()/rho_out;
        value = {uin, uin, uout, uout};
        flame.setInitialGuess("velocity",locs,value);
        value = {temp, temp, Tad, Tad};
        flame.setInitialGuess("T",locs,value);

        for (size_t i=0; i<nsp; i++) {
            value = {yin[i], yin[i], yout[i], yout[i]};
            flame.setInitialGuess(gas->speciesName(i),locs,value);
        }

        inlet.setMoleFractions(x.data());
        inlet.setMdot(mdot);
        inlet.setTemperature(temp);

        flame.showSolution();

        int flowdomain = 1;
        double ratio = 10.0;
        double slope = 0.08;
        double curve = 0.1;

        flame.setRefineCriteria(flowdomain,ratio,slope,curve);

        int loglevel=1;

        // Solve freely propagating flame

        // Linearly interpolate to find location where this temperature would
        // exist. The temperature at this location will then be fixed for
        // remainder of calculation.
        flame.setFixedTemperature(0.5 * (temp + Tad));
        flow.solveEnergyEqn();
        bool refine_grid = true;

        flame.solve(loglevel,refine_grid);
        double flameSpeed_mix = flame.value(flowdomain,
                                            flow.componentIndex("velocity"),0);
        print("Flame speed with mixture-averaged transport: {} m/s\n",
              flameSpeed_mix);

        // now switch to multicomponent transport
        flow.setTransport(*trmulti);
        flame.solve(loglevel, refine_grid);
        double flameSpeed_multi = flame.value(flowdomain,
                                              flow.componentIndex("velocity"),0);
        print("Flame speed with multicomponent transport: {} m/s\n",
              flameSpeed_multi);

        // now enable Soret diffusion
        flow.enableSoret(true);
        flame.solve(loglevel, refine_grid);
        double flameSpeed_full = flame.value(flowdomain,
                                             flow.componentIndex("velocity"),0);
        print("Flame speed with multicomponent transport + Soret: {} m/s\n",
              flameSpeed_full);

        vector_fp zvec,Tvec,COvec,CO2vec,Uvec;

        print("\n{:9s}\t{:8s}\t{:5s}\t{:7s}\n",
              "z (m)", "T (K)", "U (m/s)", "Y(CO)");
        for (size_t n = 0; n < flow.nPoints(); n++) {
            Tvec.push_back(flame.value(flowdomain,flow.componentIndex("T"),n));
            COvec.push_back(flame.value(flowdomain,
                                        flow.componentIndex("CO"),n));
            CO2vec.push_back(flame.value(flowdomain,
                                         flow.componentIndex("CO2"),n));
            Uvec.push_back(flame.value(flowdomain,
                                       flow.componentIndex("velocity"),n));
            zvec.push_back(flow.grid(n));
            print("{:9.6f}\t{:8.3f}\t{:5.3f}\t{:7.5f}\n",
                  flow.grid(n), Tvec[n], Uvec[n], COvec[n]);
        }

        print("\nAdiabatic flame temperature from equilibrium is: {}\n", Tad);
        print("Flame speed for phi={} is {} m/s.\n", phi, Uvec[0]);

        std::ofstream outfile("flamespeed.csv", std::ios::trunc);
        outfile << "  Grid,   Temperature,   Uvec,   CO,    CO2\n";
        for (size_t n = 0; n < flow.nPoints(); n++) {
            print(outfile, " {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}\n",
                  flow.grid(n), Tvec[n], Uvec[n], COvec[n], CO2vec[n]);
        }
    } catch (CanteraError& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << "program terminating." << std::endl;
        return -1;
    }
    return 0;
}

int main()
{
    double phi;
    print("Enter phi: ");
    std::cin >> phi;
    return flamespeed(phi);
}
