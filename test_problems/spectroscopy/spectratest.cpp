
#include "Cantera.h"
#include "spectra.h"
#include <iostream>
using namespace std;
using namespace CanteraSpectra;

int main() {
//     double B;
//     double T;
//     cout << "enter B, T: ";
//     cin >> B >> T;
//     Rotor* r = new Rotor(B);
//     double theta = wnum_to_J(B)/Boltzmann;
//     cout << "theta = " << theta << endl;
//     double pop, cpop = 0.0;
//     for (int j = 0; j < 50; j++) {
//         pop = r->population(j, T);
//         cpop += pop;
//         if (cpop > 0.999) break;
//         cout << j << ", " << r->energy_w(j) << ", " 
//              << r->frequency(j, j+1) << ", " 
//              << pop << ", " << cpop << ", " << r->partitionFunction(T) << ", " << T/theta << endl;
//     }

    // test line broading classes
    double gam = 2.0e0; 
    double sigma = 5.0;

    LineBroadener* lor = new Lorentzian(gam);
    LineBroadener* gaus = new Gaussian(sigma);
    Voigt* voig = new Voigt(sigma, gam);
    //voig->testv();

    double dnu = 0.1;
    double nu;
    double sum = 0.0, sumg = 0.0, sumlor = 0.0;
    for (int n = -2000; n < 2000; n++) {
        //cout << n << endl;
        nu = n*dnu;
        sumg += gaus->profile(nu)*dnu;
        sum += voig->profile(nu)*dnu;
        sumlor += lor->profile(nu)*dnu;
        cout << nu << ", " << (*lor)(nu) << ", " << (*gaus)(nu) 
             << ", " << (*voig)(nu) << endl;
    }
    cout << "Voigt area = " << sum << endl;
    cout << "Gaussian area = " << sumg << endl;
    cout << "Lorentzian area = " << sumlor << endl;
}
