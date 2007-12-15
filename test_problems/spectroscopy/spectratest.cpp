
#include "Cantera.h"
#include "spectra.h"
#include <iostream>
using namespace std;
using namespace Cantera;

int main() {

    Rotor* r = new Rotor(1.0);
    double w = 8065.0;
    cout << "eV: " << wnum_to_eV(w) << endl;
    delete r;
}
