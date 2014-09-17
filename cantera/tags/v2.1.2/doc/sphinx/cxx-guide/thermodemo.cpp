#include "cantera/thermo.h"

using namespace Cantera;

void thermo_demo(const std::string& file, const std::string& phase)
{
    ThermoPhase* gas = newPhase(file, phase);
    gas->setState_TPX(1500.0, 2.0*OneAtm, "O2:1.0, H2:3.0, AR:1.0");

    // temperature, pressure, and density
    std::cout << gas->temperature() << std::endl;
    std::cout << gas->pressure() << std::endl;
    std::cout << gas->density() << std::endl;

    // molar thermodynamic properties
    std::cout << gas->enthalpy_mole() << std::endl;
    std::cout << gas->entropy_mole() << std::endl;

    // specific (per unit mass) thermodynamic properties
    std::cout << gas->enthalpy_mass() << std::endl;
    std::cout << gas->entropy_mass() << std::endl;

    // chemical potentials of the species
    int numSpecies = gas->nSpecies();
    vector_fp mu(numSpecies);
    gas->getChemPotentials(&mu[0]);
    int n;
    for (n = 0; n < numSpecies; n++) {
        std::cout << gas->speciesName(n) << " " << mu[n] << std::endl;
    }
}

int main(int argc, char** argv)
{
    try {
        thermo_demo("h2o2.cti","ohmech");
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return 1;
    }
    return 0;
}
