#include <cantera/Cantera.h>

void demoprog() {
 // Calls Cantera
}

int main() {
    try {
        demoprog();
    }
    catch (CanteraError) {
        showErrors(cout);
    }
}

