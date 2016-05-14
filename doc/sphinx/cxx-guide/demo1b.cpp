#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/transport.h"

using namespace Cantera;

// The actual code is put into a function that can be called from the main
// program.
void simple_demo2()
{
    // Create a new phase
    std::unique_ptr<ThermoPhase> gas(newPhase("gri30.cti", "gri30_mix"));

    // List of phases participating in reactions (just one for homogeneous
    // kinetics)
    std::vector<ThermoPhase*> phases{gas.get()};

    // Create the Kinetics object. Based on the phase definition used, this will
    // be a GasKinetics object.
    std::unique_ptr<Kinetics> kin(newKineticsMgr(gas->xml(), phases));

    // Set an "interesting" mixture state where we will observe non-zero reacton
    // rates.
    gas->setState_TPX(500.0, 2.0*OneAtm, "CH4:1.0, O2:1.0, N2:3.76");
    gas->equilibrate("HP");
    gas->setState_TP(gas->temperature() - 100, gas->pressure());

    // Get the net reaction rates
    vector_fp wdot(kin->nReactions());
    kin->getNetRatesOfProgress(wdot.data());

    writelog("Net reaction rates for reactions involving CO2\n");
    size_t kCO2 = gas->speciesIndex("CO2");
    for (size_t i = 0; i < kin->nReactions(); i++) {
        if (kin->reactantStoichCoeff(kCO2, i)
            || kin->productStoichCoeff(kCO2, i)) {
            writelog("{:3d}  {:30s}  {: .8e}\n",
                i, kin->reactionString(i), wdot[i]);
        }
    }
    writelog("\n");

    // Create a Transport object. Based on the transport model specified in the
    // "gri30_mix" phase, this will be a MixGasTransport object.
    std::unique_ptr<Transport> trans(newDefaultTransportMgr(gas.get()));
    writelog("T        viscosity     thermal conductivity\n");
    writelog("------   -----------   --------------------\n");
    for (size_t n = 0; n < 5; n++) {
        double T = 300 + 100 * n;
        gas->setState_TP(T, gas->pressure());
        writelog("{:.1f}    {:.4e}    {:.4e}\n",
            T, trans->viscosity(), trans->thermalConductivity());
    }
}

// the main program just calls function simple_demo2 within a 'try' block, and
// catches exceptions that might be thrown
int main()
{
    try {
        simple_demo2();
    } catch (std::exception& err) {
        std::cout << err.what() << std::endl;
    }
}

