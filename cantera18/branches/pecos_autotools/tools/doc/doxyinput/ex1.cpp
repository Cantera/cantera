#include <cantera/Cantera.h>
int main() {
    ThermoPhase* gas = newPhase("h2o2.cti","ohmech");
    cout << gas->temperature() << endl;
    return 0;
}
