#include "cantera/thermo.h"
#include <iostream>

int main(int argc, char** argv)
{
    Cantera::ThermoPhase* gas = Cantera::newPhase("h2o2.cti","ohmech");
    std::cout << gas->temperature() << std::endl;
    return 0;
}
