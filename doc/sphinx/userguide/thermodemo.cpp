#include "cantera/core.h"
#include <iostream>

using namespace Cantera;

void thermo_demo(const string& file, const string& phase)
{
    // Create a new Solution object
    auto sol = newSolution(file, phase);
    auto gas = sol->thermo();

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
    size_t numSpecies = gas->nSpecies();
    vector<double> mu(numSpecies);
    gas->getChemPotentials(mu.data());
    for (size_t n = 0; n < numSpecies; n++) {
        std::cout << gas->speciesName(n) << " " << mu[n] << std::endl;
    }
}

int main(int argc, char** argv)
{
    try {
        thermo_demo("h2o2.yaml", "ohmech");
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return 1;
    }
    return 0;
}
