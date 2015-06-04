#include "cantera/kinetics.h"

#include <cstdio>

using namespace Cantera;
using namespace std;

int main(int argc, char** argv)
{
#ifdef _MSC_VER
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    try {
        if (argc != 2) {
            cout << "Error: no input file specified.\n"
                 "Choose either 'noxNeg.cti' or 'noxNeg_blessed.xml" << endl;
            exit(-1);
        }
        std::string infile(argv[1]);

        size_t i;
        double x[20];
        double cdot[20], ddot[20];

        XML_Node* xc = get_XML_File(infile);
        XML_Node* const xg = xc->findNameID("phase", "air");
        ThermoPhase* gasTP = newPhase(*xg);
        size_t nsp = gasTP->nSpecies();
        cout << "Number of species = " << nsp << endl;


        vector<ThermoPhase*> phaseList;
        phaseList.push_back(gasTP);
        GasKinetics* iKin_ptr = new GasKinetics();
        importKinetics(*xg, phaseList, iKin_ptr);
        size_t nr = iKin_ptr->nReactions();
        cout << "Number of reactions = " << nr << endl;

        size_t iH = gasTP->speciesIndex("H");
        size_t iO2 = gasTP->speciesIndex("O2");
        size_t iH2O = gasTP->speciesIndex("H2O");
        size_t iNH = gasTP->speciesIndex("NH");
        size_t iNO = gasTP->speciesIndex("NO");
        size_t iN2O = gasTP->speciesIndex("N2O");

        for (i = 0; i < nsp; i++) {
            x[i] = 0.0;
        }
        x[iH2O] = 1.0 /2.0;
        x[iH]   = 0.2 /2.0;
        x[iO2]  = 0.3 /2.0;
        x[iNH]  = 0.05/2.0;
        x[iNO]  = 0.05/2.0;
        x[iN2O]  = 0.05/2.0;

        double p = OneAtm;

        gasTP->setState_TPX(2000., p, x);


        double src[20];
        for (i = 0; i < 20; i++) {
            src[i] = 0.0;
        }
        iKin_ptr->getNetProductionRates(src);

        for (i = 0; i < nsp; i++) {
            string sSt = gasTP->speciesName(i);
            printf("rop [ %.4d:%s ] = %.5g \n", (int) i, sSt.c_str(), src[i]);
        }

        size_t nReactions = iKin_ptr->nReactions();
        cout << "number of reactions = " << nReactions << endl;

        double fwd_rop[20];
        double rev_rop[20];
        iKin_ptr->getFwdRatesOfProgress(fwd_rop);
        iKin_ptr->getRevRatesOfProgress(rev_rop);
        for (i = 0; i < nReactions; i++) {
            printf("fwd_rop[%3d] = %13g    rev_rop[%3d] = %13g\n", (int) i,
                   fwd_rop[i], (int) i, rev_rop[i]);
        }



        iKin_ptr->getCreationRates(cdot);
        iKin_ptr->getDestructionRates(ddot);


        return 0;
    } catch (CanteraError& err) {
        std::cerr << err.what() << std::endl;
        cerr << "program terminating." << endl;
        return -1;
    }
}
