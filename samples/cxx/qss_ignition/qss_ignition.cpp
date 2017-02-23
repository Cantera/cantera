
// A zero-D homogeneous constant-pressure auto-ignition case to demonstrate
// the usage of GasQSSKinetics class. We are interested in the ignition-delay
// time."""

#include "cantera/kinetics/GasQSSKinetics.h"
#include "cantera/kinetics/importKinetics.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/zerodim.h"

#include <fstream>

using namespace Cantera;

void runexample()
{
    // use reaction mechanism DME_SK39_LuoLu with 30 bulk and 9 qss species
    IdealGasPhase gas_bulk("DME_RED30_LuoLu.cti", "DME_SK39");
    IdealGasPhase gas_qss("DME_RED30_LuoLu.cti", "QSS");

    // decale shared-pointer for GasQSSKinetics
    std::shared_ptr<GasQSSKinetics> kin_qss = make_shared<GasQSSKinetics>();
    // create vector of thermo
    std::vector<ThermoPhase *> th_qss{&gas_bulk, &gas_qss};
    // import kinetics
    importKinetics(gas_qss.xml(), th_qss, kin_qss.get());

    // set initial conditions
    gas_bulk.setState_TPX(1500.0, OneAtm, "CH3OCH3:0.33, O2:1.00, N2:3.76");

    // create a reactor
    IdealGasConstPressureReactor r;
    r.setThermoMgr(gas_bulk);
    r.setKineticsMgr(*kin_qss);

    // create a container object to run the simulation
    // and add the reactor to it
    ReactorNet sim;
    sim.addReactor(r);

    // advance reactor till T >= T0 + 500 K
    const double dt = 1.0e-7;
    double tm = 0.0;
    while (true) {
        sim.advance(tm);
        tm += dt;
        if (gas_bulk.temperature() >= 2000.0) break;
    }

    std::ofstream f("qss_ignition_cxx.csv");
    std::vector<size_t> k_out {
        gas_bulk.speciesIndex("CH3OCH3"),
        gas_bulk.speciesIndex("O2"),
        gas_bulk.speciesIndex("CO2"),
        gas_bulk.speciesIndex("H2O"),
        gas_bulk.speciesIndex("CO"),
        gas_bulk.speciesIndex("OH"),
        gas_bulk.speciesIndex("H")
    };
    f << tm << ", "
      << gas_bulk.temperature() << ", ";
    for (size_t i = 0; i < k_out.size(); i++) {
        f << gas_bulk.moleFraction(k_out[i]) << ", ";
    }
    f << std::endl;
}

int main()
{
    try {
        runexample();
        return 0;
    } catch (CanteraError& err) {
        // handle exceptions thrown by Cantera
        std::cout << err.what() << std::endl;
        std::cout << " terminating... " << std::endl;
        appdelete();
        return 1;
    }
}
