/**
 *  @file runDiamond.cpp
 */

#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include <iostream>

using namespace std;
using namespace Cantera;

int main(int argc, char** argv)
{
#if defined(_MSC_VER) && _MSC_VER < 1900
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    int i, k;

    try {
        auto gas = newThermo("diamond.yaml", "gas");
        size_t nsp = gas->nSpecies();
        cout.precision(4);
        cout << "Number of species = " << nsp << endl;

        auto diamond = newThermo("diamond.yaml", "diamond");
        size_t nsp_diamond = diamond->nSpecies();
        cout << "Number of species in diamond = " << nsp_diamond << endl;

        auto diamond100 = newThermo("diamond.yaml", "diamond_100");
        size_t nsp_d100 = diamond100->nSpecies();
        cout << "Number of species in diamond_100 = " << nsp_d100 << endl;

        auto kin = newKinetics({gas, diamond, diamond100},
                               "diamond.yaml", "diamond_100");
        InterfaceKinetics& ikin = dynamic_cast<InterfaceKinetics&>(*kin);
        size_t nr = kin->nReactions();
        cout << "Number of reactions = " << nr << endl;

        double x[20];
        for (i = 0; i < 20; i++) {
            x[i] = 0.0;
        }
        x[0] = 0.0010;
        x[1] = 0.9888;
        x[2] = 0.0002;
        x[3] = 0.0100;
        double p = 20.0*OneAtm/760.0;

        gas->setState_TPX(1200., p, x);

        for (i = 0; i < 20; i++) {
            x[i] = 0.0;
        }
        size_t i0 = diamond100->speciesIndex("c6H*");
        x[i0] = 0.1;
        size_t i1 = diamond100->speciesIndex("c6HH");
        x[i1] = 0.9;
        diamond100->setState_TX(1200., x);

        for (i = 0; i < 20; i++) {
            x[i] = 0.0;
        }
        x[0] = 1.0;
        diamond100->setState_TP(1200., p);

        ikin.advanceCoverages(100.);

        // Throw some asserts in here to test that they compile
        AssertTrace(p == p);
        AssertThrow(p == p, "main");
        AssertThrowMsg(i == 20, "main", "are you kidding");

        double src[20];
        for (i = 0; i < 20; i++) {
            src[i] = 0.0;
        }
        kin->getNetProductionRates(src);
        double sum = 0.0;
        double naH = 0.0;
        for (k = 0; k < 13; k++) {
            if (k < 4) {
                naH = gas->nAtoms(k, 0);
            } else if (k == 4) {
                naH = 0;
            } else if (k > 4) {
                int itp = k - 5;
                naH = diamond100->nAtoms(itp, 0);
            }
            writelog("{} {} {}\n", k, naH, src[k]);
            sum += naH * src[k];
        }

        writelog("sum = {}\n", sum);
        double mwd = diamond->molecularWeight(0);
        double dens = diamond->density();
        double gr = src[4] * mwd / dens;
        gr *= 1.0E6 * 3600.;
        cout << "growth rate = " << gr << " microns per hour" << endl;

        diamond100->getMoleFractions(x);
        cout << "Coverages:" << endl;
        for (k = 0; k < 8; k++) {
            cout << k << "   " << diamond100->speciesName(k)
                 << "   "
                 << x[k] << endl;
        }
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
    }

    return 0;
}
