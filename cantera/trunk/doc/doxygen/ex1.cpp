#include "cantera/thermo.h"
#include <iostream>

int main()
{
    ThermoPhase* gas = newPhase("h2o2.cti","ohmech");
    cout << gas->temperature() << endl;
    return 0;
}
