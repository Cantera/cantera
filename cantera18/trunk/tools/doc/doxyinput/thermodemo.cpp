
#include <cantera/Cantera.h>

void thermo_demo(string file, string phase) {
    ThermoPhase* gas = newPhase(file, phase);
    gas->setState_TPX(1500.0, 2.0*OneAtm, "O2:1.0, H2:3.0, AR:1.0");

    // temperature, pressure, and density
    cout << gas->temperature() << endl;
    cout << gas->pressure() << endl;
    cout << gas->density() << endl;

    // molar thermodynamic properties
    cout << gas->enthalpy_mole() << endl;
    cout << gas->entropy_mole() << endl;

    // specific (per unit mass) thermodynamic properties
    cout << gas->enthalpy_mass() << endl;
    cout << gas->entropy_mass() << endl;

    // chemical potentials of the species
    int numSpecies = gas->nSpecies();
    vector_fp mu(numSpecies);
    gas->getChemPotentials(mu.begin());
    int n;
    for (n = 0; n < numSpecies; n++) {
        cout << gas->speciesName(n) << " " << mu[n] << endl;
    }
}
 
int main() {

    try {
        thermo_demo("h2o2.cti","ohmech");
    }
    catch (CanteraError) {
        showErrors();
    }
}

