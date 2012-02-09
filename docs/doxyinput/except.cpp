#include <cantera/Cantera.h>
//
// artifical example of throwing and catching a CanteraError exception.
//
void mycode() {
    ThermoPhase* gas = newPhase("h2o2.cti","ohmech");
    if (gas->temperature() < 3000.0) {
        throw CanteraError("mycode","test of exception throwing and catching");
    }
}

int main() {
    try {
        mycode();
    }
    catch (CanteraError) {
        showErrors();
        error("program terminating.");
    }
}
