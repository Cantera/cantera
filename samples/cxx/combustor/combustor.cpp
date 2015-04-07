
// A combustor. Two separate stream - one pure methane and the other
// air, both at 300 K and 1 atm flow into an adiabatic combustor where
// they mix. We are interested in the steady-state burning
// solution. Since at 300 K no reaction will occur between methane and
// air, we need to use an 'igniter' to initiate the chemistry. A simple
// igniter is a pulsed flow of atomic hydrogen. After the igniter is
// turned off, the system approaches the steady burning solution."""

#include "cantera/zerodim.h"
#include "cantera/IdealGasMix.h"

#include <fstream>

using namespace Cantera;

void runexample()
{

    // use reaction mechanism GRI-Mech 3.0
    IdealGasMix gas("gri30.cti", "gri30");
    int nsp = gas.nSpecies();

    // create a reservoir for the fuel inlet, and set to pure methane.
    Reservoir fuel_in;
    gas.setState_TPX(300.0, OneAtm, "CH4:1.0");
    fuel_in.insert(gas);
    double fuel_mw = gas.meanMolecularWeight();

    // create a reservoir for the air inlet
    Reservoir air_in;
    IdealGasMix air("air.cti");
    gas.setState_TPX(300.0, OneAtm, "N2:0.78, O2:0.21, AR:0.01");
    double air_mw = air.meanMolecularWeight();
    air_in.insert(gas);

    // to ignite the fuel/air mixture, we'll introduce a pulse of radicals.
    // The steady-state behavior is independent of how we do this, so we'll
    // just use a stream of pure atomic hydrogen.
    gas.setState_TPX(300.0, OneAtm, "H:1.0");
    Reservoir igniter;
    igniter.insert(gas);


    // create the combustor, and fill it in initially with N2
    gas.setState_TPX(300.0, OneAtm, "N2:1.0");
    Reactor combustor;
    combustor.insert(gas);
    combustor.setInitialVolume(1.0);


    // create a reservoir for the exhaust. The initial composition
    // doesn't matter.
    Reservoir exhaust;
    exhaust.insert(gas);


    // lean combustion, phi = 0.5
    double equiv_ratio = 0.5;

    // compute fuel and air mass flow rates
    double factor = 0.1;
    double air_mdot = factor*9.52*air_mw;
    double fuel_mdot = factor*equiv_ratio*fuel_mw;

    // create and install the mass flow controllers. Controllers
    // m1 and m2 provide constant mass flow rates, and m3 provides
    // a short Gaussian pulse only to ignite the mixture
    MassFlowController m1;
    m1.install(fuel_in, combustor);
    m1.setMassFlowRate(fuel_mdot);

    // Now create the air mass flow controller.  Note that this connects
    // two reactors with different reaction mechanisms and different
    // numbers of species. Downstream and upstream species are matched by
    // name.
    MassFlowController m2;
    m2.install(air_in, combustor);
    m2.setMassFlowRate(air_mdot);


    // The igniter will use a Gaussian 'functor' object to specify the
    // time-dependent igniter mass flow rate.
    double A = 0.1;
    double FWHM = 0.2;
    double t0 = 0.5;
    Gaussian igniter_mdot(A, t0, FWHM);

    MassFlowController m3;
    m3.install(igniter, combustor);
    m3.setFunction(&igniter_mdot);

    // put a valve on the exhaust line to regulate the pressure
    Valve v;
    v.install(combustor, exhaust);
    double Kv = 1.0;
    v.setParameters(1, &Kv);

    // the simulation only contains one reactor
    ReactorNet sim;
    sim.addReactor(combustor);

    // take single steps to 6 s, writing the results to a CSV file
    // for later plotting.
    double tfinal = 1.0;
    double tnow = 0.0;
    double tres;
    int k;

    std::ofstream f("combustor_cxx.csv");
    std::vector<size_t> k_out;
    k_out.push_back(gas.speciesIndex("CH4"));
    k_out.push_back(gas.speciesIndex("O2"));
    k_out.push_back(gas.speciesIndex("CO2"));
    k_out.push_back(gas.speciesIndex("H2O"));
    k_out.push_back(gas.speciesIndex("CO"));
    k_out.push_back(gas.speciesIndex("OH"));
    k_out.push_back(gas.speciesIndex("H"));
    k_out.push_back(gas.speciesIndex("C2H6"));

    while (tnow < tfinal) {
        tnow += 0.005;
        sim.advance(tnow);
        tres = combustor.mass()/v.massFlowRate();
        f << tnow << ", "
          << combustor.temperature() << ", "
          << tres << ", ";
        ThermoPhase& c = combustor.contents();
        for (size_t i = 0; i < k_out.size(); i++) {
            f << c.moleFraction(k_out[i]) << ", ";
        }
        f << std::endl;
    }
    f.close();
}

int main()
{

    try {
        runexample();
        return 0;
    }
    // handle exceptions thrown by Cantera
    catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        std::cout << " terminating... " << std::endl;
        appdelete();
        return 1;
    }
}
