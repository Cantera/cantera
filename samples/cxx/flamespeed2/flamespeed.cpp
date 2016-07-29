/*!
 * @file flamespeed.cpp
 * C++ demo program to compute flame speeds using GRI-Mech.
 */

#include "cantera/oneD/Sim1D.h"
#include "cantera/oneD/Inlet1D.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/IdealGasMix.h"
#include "cantera/transport.h"
#include <fstream>

using namespace Cantera;
using fmt::print;

int flamespeed(double phi)
{
    try {
        IdealGasMix gas("gri30_ion.xml","gri30_ion_mix");

        doublereal temp = 300.0; // K
        doublereal pressure = 1.0*OneAtm; //atm
        doublereal uin = 0.3; //m/sec

        size_t nsp = gas.nSpecies();
        vector_fp x(nsp, 0.0);

        doublereal C_atoms = 1.0;
        doublereal H_atoms = 4.0;
        doublereal ax = C_atoms + H_atoms / 4.0;
        doublereal fa_stoic = 1.0 / (4.76 * ax);
        x[gas.speciesIndex("CH4")] = 1.0;
        x[gas.speciesIndex("O2")] = 0.21 / phi / fa_stoic;
        x[gas.speciesIndex("N2")] = 0.79 / phi/ fa_stoic;

        gas.setState_TPX(temp, pressure, x.data());
        doublereal rho_in = gas.density();

        vector_fp yin(nsp);
        gas.getMassFractions(&yin[0]);

        gas.equilibrate("HP");
        vector_fp yout(nsp);
        gas.getMassFractions(&yout[0]);
        doublereal rho_out = gas.density();
        doublereal Tad = gas.temperature();
        print("phi = {}, Tad = {}\n", phi, Tad);

        //=============  build each domain ========================


        //-------- step 1: create the flow -------------

        FreeFlame flow(&gas);

        // create an initial grid
        int nz = 6;
        doublereal lz = 0.1;
        vector_fp z(nz);
        doublereal dz = lz/((doublereal)(nz-1));
        for (int iz = 0; iz < nz; iz++) {
            z[iz] = ((doublereal)iz)*dz;
        }

        flow.setupGrid(nz, &z[0]);

        // specify the objects to use to compute kinetic rates and
        // transport properties

        std::unique_ptr<Transport> trmix(newTransportMgr("Mix", &gas));
        std::unique_ptr<Transport> trmulti(newTransportMgr("Multi", &gas));

        flow.setTransport(*trmix);
        flow.setKinetics(gas);
        flow.setPressure(pressure);

        //------- step 2: create the inlet  -----------------------

        Inlet1D inlet;

        inlet.setMoleFractions(x.data());
        doublereal mdot=uin*rho_in;
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
        flame.setInitialGuess("u",locs,value);
        value = {temp, temp, Tad, Tad};
        flame.setInitialGuess("T",locs,value);

        for (size_t i=0; i<nsp; i++) {
            value = {yin[i], yin[i], yout[i], yout[i]};
            flame.setInitialGuess(gas.speciesName(i),locs,value);
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
        double flameSpeed_mix = flame.value(flowdomain,flow.componentIndex("u"),0);
        print("Flame speed with mixture-averaged transport: {} m/s\n",
              flameSpeed_mix);

        // now enable ambipolar diffusion
        flow.enableAmbipolar(true);
        flame.solve(loglevel, refine_grid);
        double flameSpeed_full = flame.value(flowdomain,flow.componentIndex("u"),0);
        print("Flame speed with ambipolar diffusion: {} m/s\n",
              flameSpeed_full);

        vector_fp zvec,Tvec,HCOxvec,H3Oxvec,Evec,Uvec;

        print("\n{:9s}\t{:8s}\t{:5s}\t{:7s}\n",
              "z (m)", "T (K)", "U (m/s)", "Y(CO)");
        for (size_t n = 0; n < flow.nPoints(); n++) {
            Tvec.push_back(flame.value(flowdomain,flow.componentIndex("T"),n));
            HCOxvec.push_back(flame.value(flowdomain,flow.componentIndex("HCO+"),n));
            H3Oxvec.push_back(flame.value(flowdomain,flow.componentIndex("H3O+"),n));
            Evec.push_back(flame.value(flowdomain,flow.componentIndex("E"),n));
            Uvec.push_back(flame.value(flowdomain,flow.componentIndex("u"),n));
            zvec.push_back(flow.grid(n));
            print("{:9.6f}\t{:8.3f}\t{:5.3f}\t{:7.5f}\n",
                  flow.grid(n), Tvec[n], Uvec[n], HCOxvec[n]);
        }

        print("\nAdiabatic flame temperature from equilibrium is: {}\n", Tad);
        print("Flame speed for phi={} is {} m/s.\n", phi, Uvec[0]);

        std::ofstream outfile("flamespeed.csv", std::ios::trunc);
        outfile << "  Grid,   Temperature,   Uvec,   HCO+,    H3O+,   E\n";
        for (size_t n = 0; n < flow.nPoints(); n++) {
            print(outfile, " {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}, {:11.3e}\n",
                  flow.grid(n), Tvec[n], Uvec[n], HCOxvec[n], H3Oxvec[n], Evec[n]);
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
