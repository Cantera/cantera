/**
 *  @file runDiamond.cpp
 */

//  Example
//
// Note that this example needs updating. It works fine, but is
// written in a way that is less than transparent or
// user-friendly. This could be rewritten using class Interface to
// make things simpler.

#include "cantera/kinetics.h"
#include <cstdio>

using namespace std;
using namespace Cantera;

void printDbl(double val)
{
    if (fabs(val) < 5.0E-200) {
        cout << " nil";
    } else {
        cout << val;
    }
}

int main(int argc, char** argv)
{
#ifdef _MSC_VER
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    string infile = "frac.xml";
    double x[10], kc[10];
    double cdot[10], ddot[10];
    //double fwd_rop[10];
    try {
        XML_Node* xc = new XML_Node();
        string path = findInputFile(infile);
        ctml::get_CTML_Tree(xc, path);

        XML_Node* const xg = xc->findNameID("phase", "gas");
        ThermoPhase* gasTP = newPhase(*xg);
        size_t nsp = gasTP->nSpecies();
        cout.precision(4);
        cout << "Number of species = " << nsp << endl;

        vector<ThermoPhase*> phaseList;
        phaseList.push_back(gasTP);
        GasKinetics* iKin_ptr = new GasKinetics();
        importKinetics(*xg, phaseList, iKin_ptr);
        size_t nr = iKin_ptr->nReactions();
        cout << "Number of reactions = " << nr << endl;

        size_t iH2 = gasTP->speciesIndex("H2");
        size_t iH = gasTP->speciesIndex("H");
        size_t iO2 = gasTP->speciesIndex("O2");
        size_t iOH = gasTP->speciesIndex("OH");
        size_t iH2O = gasTP->speciesIndex("H2O");

        for (size_t i = 0; i < nsp; i++) {
            x[i] = 0.0;
        }
        x[iH2O] = 1.0/2.0;
        x[iOH] = 0.1/2.0;
        x[iH]  = 0.2/2.0;
        x[iO2] = 0.3/2.0;
        x[iH2] = 0.4/2.0;

        double p = OneAtm;

        gasTP->setState_TPX(2000., p, x);


        double src[20];
        for (size_t i = 0; i < 20; i++) {
            src[i] = 0.0;
        }
        iKin_ptr->getNetProductionRates(src);

        double fwd_rop[10];
        iKin_ptr->getFwdRatesOfProgress(fwd_rop);
        cout << "fwd_rop[0] = " << fwd_rop[0] << endl;
        cout << "fwd_rop[1] = " << fwd_rop[1] << endl;

        iKin_ptr->getCreationRates(cdot);
        iKin_ptr->getDestructionRates(ddot);

        for (size_t k = 0; k < nsp; k++) {
            string sss = gasTP->speciesName(k);
            cout << k << "  " << sss << "  ";
            printDbl(src[k]);
            cout << endl;
        }

        printf("Creation Rates: \n");
        for (size_t k = 0; k < nsp - 1; k++) {
            string sss = gasTP->speciesName(k);
            cout << k << "  " << sss << "  ";
            cout << cdot[k] << "  ";
            cout << cdot[k] / fwd_rop[0] << " ";
            cout << endl;
        }
        string sss = gasTP->speciesName(iH2O);
        cout << iH2O << "  " << sss << "  ";
        cout << cdot[iH2O] << "  ";
        cout << cdot[iH2O] / fwd_rop[1] << " ";
        cout << endl;


        printf("Destruction Rates: \n");
        for (size_t k = 0; k < nsp-1; k++) {
            string sss = gasTP->speciesName(k);
            cout << k << "  " << sss << "  ";
            cout << ddot[k] << "  ";
            cout << ddot[k] / fwd_rop[1] << " ";
            cout << endl;
        }
        sss = gasTP->speciesName(iH2O);
        cout << iH2O << "  " << sss << "  ";
        cout << ddot[iH2O] << "  ";
        cout << ddot[iH2O] / fwd_rop[0] << " ";
        cout << endl;


        double c[10];
        gasTP->getConcentrations(c);

        double order_H2 = 0.8;
        double order_OH = 2.0;
        double order_O2 = 1.0;

        double kf[10];
        iKin_ptr->getFwdRateConstants(kf);
        printf("kf[0] = %g\n", kf[0]);
        printf("kf[1] = %g\n", kf[1]);

        //double cprod0 = c[iH2O];
        double cprod1 = pow(c[iH2], order_H2) * pow(c[iOH], order_OH) * pow(c[iO2], order_O2);

        printf("equal numbers 0: %g %g \n", kf[0] * c[iH2O], fwd_rop[0]);

        printf("equal numbers 1: %g %g\n", kf[1] * cprod1, fwd_rop[1]);

        iKin_ptr->getEquilibriumConstants(kc);

        printf("Equilibrium constants for irreversible fractional rxns:\n");
        printf("Kc[0] = %g\n", kc[0]);
        printf("Kc[1] = %g\n", kc[1]);

        delete iKin_ptr;
        iKin_ptr = 0;
        delete gasTP;
        delete xc;
        appdelete();


    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return -1;
    }

    return 0;
}
/***********************************************************/
