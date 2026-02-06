#include "cantera/core.h"
#include "cantera/kinetics/Reaction.h"
#include <iostream>

using namespace Cantera;

void simple_demo2()
{
    auto sol = newSolution("gri30.yaml", "gri30");
    auto gas = sol->thermo();

    // Access the Kinetics object. Based on the phase definition in the input file,
    // this will reference a BulkKinetics object.
    auto kin = sol->kinetics();

    // Set an "interesting" state where we will observe non-zero reaction rates.
    gas->setState_TPX(500.0, 2.0*OneAtm, "CH4:1.0, O2:1.0, N2:3.76");
    gas->equilibrate("HP");
    gas->setState_TP(gas->temperature() - 100, gas->pressure());

    // Get the net reaction rates.
    vector<double> wdot(kin->nReactions());
    kin->getNetRatesOfProgress(wdot);

    writelog("Net reaction rates for reactions involving CO2\n");
    size_t kCO2 = gas->speciesIndex("CO2");
    for (size_t i = 0; i < kin->nReactions(); i++) {
        if (kin->reactantStoichCoeff(kCO2, i) || kin->productStoichCoeff(kCO2, i)) {
            auto rxn = kin->reaction(i);
            writelog("{:3d}  {:30s}  {: .8e}\n", i, rxn->equation(), wdot[i]);
        }
    }
    writelog("\n");

    // Access the Transport object. Based on the phase definition in the input file,
    // this will reference a MixTransport object.
    auto trans = sol->transport();
    writelog("T        viscosity     thermal conductivity\n");
    writelog("------   -----------   --------------------\n");
    for (size_t n = 0; n < 5; n++) {
        double T = 300 + 100 * n;
        gas->setState_TP(T, gas->pressure());
        writelog("{:.1f}    {:.4e}    {:.4e}\n",
                 T, trans->viscosity(), trans->thermalConductivity());
    }
}

int main()
{
    try {
        simple_demo2();
    } catch (std::exception& err) {
        std::cout << err.what() << std::endl;
        return 1;
    }
    return 0;
}
