
#include "Cantera.h"
#include "spectra.h"
#include <iostream>
using namespace std;
using namespace Cantera;

int main() {
    double B;
    double T;
    cout << "enter B, T: ";
    cin >> B >> T;
    Rotor* r = new Rotor(B);
    double theta = wnum_to_J(B)/Boltzmann;
    cout << "theta = " << theta << endl;
    double pop, cpop = 0.0;
    for (int j = 0; j < 50; j++) {
        pop = r->population(j, T);
        cpop += pop;
        if (cpop > 0.999) break;
        cout << j << ", " << r->energy_w(j) << ", " 
             << r->frequency(j, j+1) << ", " 
             << pop << ", " << cpop << ", " << r->partitionFunction(T) << ", " << T/theta << endl;
    }
}
