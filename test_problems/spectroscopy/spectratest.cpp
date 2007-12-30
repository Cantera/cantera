
#include "Cantera.h"
#include "spectra.h"
#include "kernel/Nuclei.h"
#include <iostream>
using namespace std;
using namespace CanteraSpectra;

int main() {

    Nucleus* a = CanteraSpectra::HydrogenNucleus();
    Nucleus* b = HydrogenNucleus();
    if (*a == *b) {
        cout << "a and b and indistinguishable" << endl;
    }

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
