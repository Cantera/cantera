//
//     Replace this sample main program with your program
//
//

#include "Cantera.h"
#include "IdealGasMix.h"
#include "equilibrium.h"

using namespace Cantera;

int main() {

    IdealGasMix gas("h2o2.xml");
    double temp = 1200.0;
    double pres = OneAtm;
    gas.setState_TPX(temp, pres, "H2:2, O2:1, OH:0.01, H:0.01, O:0.01");
    equilibrate(gas,"HP");

    printf("\n\n**** C++ Test Program ****\n\n");


    //    Thermodynamic properties

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
    double* q = new double[irxns];

    gas.getFwdRatesOfProgress(qf);
    gas.getRevRatesOfProgress(qr);
    gas.getNetRatesOfProgress(q);

    for (int i = 0; i < irxns; i++) {
        printf("%30s %14.5g %14.5g %14.5g \n", 
            gas.reactionString(i).c_str(), qf[i], qr[i], q[i]);
    }

}
      
