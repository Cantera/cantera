///////////////////////////////////////////////////////////////////////
//
//     This demonstration program builds an object representing a
//     reacting gas mixture, and uses it to compute thermodynamic
//     properties, chemical equilibrium, and transport properties.
//
///////////////////////////////////////////////////////////////////////

// Include cantera header files. They should be included in the form
// "cantera/*.h". These headers are designed for use in C++ programs and
// provide a simplified interface to the Cantera header files. If you need
// to include core headers directly, use the format "cantera/module/*.h".

#include "cantera/IdealGasMix.h"    // defines class IdealGasMix
#include "cantera/transport.h"      // transport properties

// All Cantera kernel names are in namespace Cantera. You can either
// reference everything as Cantera::<name>, or include the following
// 'using namespace' line.
using namespace Cantera;

// The program is put into a function so that error handling code can
// be conveniently put around the whole thing. See main() below.

void demoprog()
{

    printf("\n\n**** C++ Test Program ****\n\n");

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
    gas.equilibrate("HP");

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
    vector_fp qf(irxns);
    vector_fp qr(irxns);
    vector_fp q(irxns);

    // since the gas has been set to an equilibrium state, the forward
    // and reverse rates of progress should be equal for all
    // reversible reactions, and the net rates should be zero.
    gas.getFwdRatesOfProgress(&qf[0]);
    gas.getRevRatesOfProgress(&qr[0]);
    gas.getNetRatesOfProgress(&q[0]);

    printf("\n\n");
    for (int i = 0; i < irxns; i++) {
        printf("%30s %14.5g %14.5g %14.5g  kmol/m3/s\n",
               gas.reactionString(i).c_str(), qf[i], qr[i], q[i]);
    }


    // transport properties

    // create a transport manager for the gas that computes
    // mixture-averaged properties
    std::auto_ptr<Transport> tr(newTransportMgr("Mix", &gas, 1));

    // print the viscosity, thermal conductivity, and diffusion
    // coefficients
    printf("\n\nViscosity:            %14.5g Pa-s\n", tr->viscosity());
    printf("Thermal conductivity: %14.5g W/m/K\n", tr->thermalConductivity());

    int nsp = gas.nSpecies();
    vector_fp diff(nsp);
    tr->getMixDiffCoeffs(&diff[0]);
    int k;
    printf("\n\n%20s  %26s\n", "Species","Diffusion Coefficient");
    for (k = 0; k < nsp; k++) {
        printf("%20s  %14.5g m2/s \n", gas.speciesName(k).c_str(), diff[k]);
    }
}



int main()
{

    try {
        demoprog();
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
    }
}

