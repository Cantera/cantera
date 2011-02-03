
#include <cantera/Cantera.h>
#include <cantera/IdealGasMix.h>    // defines class IdealGasMix
#include <cantera/equilibrium.h>    // chemical equilibrium
#include <cantera/transport.h>      // transport properties

void demoprog() {

    IdealGasMix gas("h2o2.cti","ohmech");
    double temp = 1200.0;
    double pres = OneAtm;
    gas.setState_TPX(temp, pres, "H2:1, O2:1, AR:2");

    //    Thermodynamic properties

    printf("\n\nInitial state:\n\n");
    printf(
        "Temperature:    %14.5g K\n"
        "Pressure:       %14.5g Pa\n"
        "Density:        %14.5g kg/m3\n"
        "Molar Enthalpy: %14.5g J/kmol\n"
        "Molar Entropy:  %14.5g J/kmol-K\n"
        "Molar cp:       %14.5g J/kmol-K\n",
        gas.temperature(), gas.pressure(), gas.density(),
        gas.enthalpy_mole(), gas.entropy_mole(), gas.cp_mole());
        
    // set the gas to the equilibrium state with the same specific
    // enthalpy and pressure
    equilibrate(gas,"HP");

    printf("\n\nEquilibrium state:\n\n");
    printf(
        "Temperature:    %14.5g K\n"
        "Pressure:       %14.5g Pa\n"
        "Density:        %14.5g kg/m3\n"
        "Molar Enthalpy: %14.5g J/kmol\n"
        "Molar Entropy:  %14.5g J/kmol-K\n"
        "Molar cp:       %14.5g J/kmol-K\n",
        gas.temperature(), gas.pressure(), gas.density(),
        gas.enthalpy_mole(), gas.entropy_mole(), gas.cp_mole());


    //   Reaction information
 
    int irxns = gas.nReactions();
    double* qf = new double[irxns];
    double* qr = new double[irxns];
    double* q  = new double[irxns];

    // since the gas has been set to an equilibrium state, the forward
    // and reverse rates of progress should be equal for all
    // reversible reactions, and the net rates should be zero.
    gas.getFwdRatesOfProgress(qf);
    gas.getRevRatesOfProgress(qr);
    gas.getNetRatesOfProgress(q);

    printf("\n\n");
    for (int i = 0; i < irxns; i++) {
        printf("%30s %14.5g %14.5g %14.5g  kmol/m3/s\n", 
            gas.reactionString(i).c_str(), qf[i], qr[i], q[i]);
    }


    // transport properties

    Transport* tr = newTransportMgr("Mix", &gas, 1);

    printf("\n\nViscosity:            %14.5g Pa-s\n", tr->viscosity());
    printf("Thermal conductivity: %14.5g W/m/K\n", tr->thermalConductivity());

    int nsp = gas.nSpecies();
    double* diff = new double[nsp];
    tr->getMixDiffCoeffs(diff);
    int k;
    printf("\n\n%20s  %26s\n", "Species","Diffusion Coefficient");
    for (k = 0; k < nsp; k++) {
        printf("%20s  %14.5g m2/s \n", gas.speciesName(k).c_str(), diff[k]);
    }

    // clean up
    delete qf;
    delete qr;
    delete q;
    delete diff;
    delete tr;
}
     

 
int main() {

    try {
        demoprog();
    }
    catch (CanteraError) {
        showErrors(cout);
    }
}

