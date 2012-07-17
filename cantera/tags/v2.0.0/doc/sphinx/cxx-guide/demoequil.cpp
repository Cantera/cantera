#include "cantera/equilibrium.h"
#include "cantera/thermo.h"

using namespace Cantera;

void equil_demo()
{
    ThermoPhase* gas = newPhase("h2o2.cti","ohmech");
    gas->setState_TPX(1500.0, 2.0*OneAtm, "O2:1.0, H2:3.0, AR:1.0");
    equilibrate(*gas, "TP");
    std::cout << gas->report() << std::endl;
}

int main()
{

    try {
        equil_demo();
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
    }
}

