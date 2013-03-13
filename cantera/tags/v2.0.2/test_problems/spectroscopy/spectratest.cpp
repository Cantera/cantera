#include "cantera/spectra.h"
#include "spectra/Nuclei.h"
#include <iostream>
using namespace std;
using namespace Cantera;

int main()
{

    Nucleus* a = HydrogenNucleus();
    Nucleus* b = HydrogenNucleus();
    if (*a == *b) {
        cout << "a and b and indistinguishable" << endl;
    } else {
        cout << "\nwhy are a and b not indistinguishable?\n";
        return 1;
    }

    // test line broading classes
    double gam = 2.0e0;
    double sigma = 5.0;

    LineBroadener* lor = new LorentzianProfile(gam);
    LineBroadener* gaus = new GaussianProfile(sigma);
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
        //cout << nu << ", " << (*lor)(nu) << ", " << (*gaus)(nu)
        //     << ", " << (*voig)(nu) << endl;
    }

    /* old output
    cout << "Voigt area = " << sum << endl;
    cout << "Gaussian area = " << sumg << endl;
    cout << "Lorentzian area = " << sumlor << endl;
    */

    // 'blessed' output:
    // Voigt area = 0.99363
    // Gaussian area = 1
    // Lorentzian area = 0.993634

    // guessing a sane tolerance
    double TOL = .0001;

    if (abs(sum-.99363) > TOL) {
        cout << "\nVOIGT AREA REGRESSION TEST FAILURE\n";
        return 1;
    }

    if (abs(sumg-1.0) > TOL) {
        cout << "\nGAUSSIAN AREA REGRESSION TEST FAILURE\n";
        return 1;
    }
    if (abs(sumlor-.993634) > TOL) {
        cout << "\nLORENTZIAN AREA REGRESSION TEST FAILURE\n";
        return 1;
    }

    // steady as she goes
    return 0;

}
