/**
 *  @file runDiamond.cpp
 */

#include "cantera/core.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/kinetics/KineticsFactory.h"
#include <iostream>

using namespace std;
using namespace Cantera;

int main(int argc, char** argv)
{
    int i, k;

    try {
        CanteraError::setStackTraceDepth(20);
        auto iface = newInterface("diamond.yaml", "diamond_100");
        auto gas = iface->adjacent("gas")->thermo();
        size_t nsp = gas->nSpecies();
        cout.precision(4);
        cout << "Number of species = " << nsp << endl;

        auto diamond = iface->adjacent("diamond")->thermo();
        size_t nsp_diamond = diamond->nSpecies();
        cout << "Number of species in diamond = " << nsp_diamond << endl;

        auto diamond100 = iface->thermo();
        size_t nsp_d100 = diamond100->nSpecies();
        cout << "Number of species in diamond_100 = " << nsp_d100 << endl;

        auto ikin = iface->kinetics();
        size_t nr = ikin->nReactions();
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
        diamond100->setMoleFractions(x);
        diamond100->setTemperature(1200.);

        for (i = 0; i < 20; i++) {
            x[i] = 0.0;
        }
        x[0] = 1.0;
        diamond100->setState_TP(1200., p);

        ikin->advanceCoverages(100.);

        // Throw some asserts in here to test that they compile
        AssertTrace(p == p);
        AssertThrow(p == p, "main");
        AssertThrowMsg(i == 20, "main", "are you kidding");

        double src[20];
        for (i = 0; i < 20; i++) {
            src[i] = 0.0;
        }
        ikin->getNetProductionRates(src);
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
