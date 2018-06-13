#include "cantera/thermo/HMWSoln.h"

#include <cstdio>

using namespace std;
using namespace Cantera;

void pAtable(HMWSoln* HMW)
{
    size_t nsp = HMW->nSpecies();
    double acMol[100];
    double mf[100];
    double activities[100];
    double moll[100];

    HMW->getMolalityActivityCoefficients(acMol);
    HMW->getMoleFractions(mf);
    HMW->getActivities(activities);
    HMW->getMolalities(moll);
    string sName;
    printf("            Name      Activity  ActCoeffMolal "
           "   MoleFract      Molality\n");
    for (size_t k = 0; k < nsp; k++) {
        sName = HMW->speciesName(k);
        printf("%16s %13.4f %13.4f %13.4f %13.4f\n",
               sName.c_str(), activities[k], acMol[k], mf[k], moll[k]);
    }

}

int main(int argc, char** argv)
{
#if defined(_MSC_VER) && _MSC_VER < 1900
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

    int retn = 0;
    try {

        HMWSoln* HMW = new HMWSoln("HMW_NaCl_tc.xml");

        size_t nsp = HMW->nSpecies();

        double mu0[100];
        double moll[100];
        string sName;

        HMW->getMolalities(moll);
        moll[1] = 6.0997;
        moll[2] = 2.1628E-9;
        moll[3] = 6.0997;
        moll[4] =1.3977E-6;
        double Temp = 150 + 273.15;
        /*
         * Equalize charge balance and dump into Cl-
         */
        double sum = -moll[1] + moll[2] + moll[3] - moll[4];
        moll[1] += sum;

        HMW->setState_TPM(Temp, OneAtm, moll);

        printf("       Temperature = %.4f K\n", Temp);
        HMW->printCoeffs();
        pAtable(HMW);

        HMW->setState_TP(Temp, 1.01325E5);
        HMW->getStandardChemPotentials(mu0);
        // translate from J/kmol to kJ/gmol
        for (size_t k = 0; k < nsp; k++) {
            mu0[k] *= 1.0E-6;
        }

        printf("           Species   Standard chemical potentials (kJ/gmol) \n");
        printf("------------------------------------------------------------\n");
        for (size_t k = 0; k < nsp; k++) {
            sName = HMW->speciesName(k);
            printf("%16s %16.9g\n", sName.c_str(), mu0[k]);
        }
        printf("------------------------------------------------------------\n");
        printf(" Some DeltaSS values:               Delta(mu_0)\n");
        double deltaG;
        size_t i1, i2, j1;
        double RT = 8.314472E-3 * 298.15;

        i1 = HMW->speciesIndex("Na+");
        i2 = HMW->speciesIndex("Cl-");
        deltaG = -432.6304 - mu0[i1] - mu0[i2];
        printf(" NaCl(S): Na+ + Cl- -> NaCl(S): %14.5g kJ/gmol \n",
               deltaG);
        printf("                                : %14.5g (dimensionless) \n",
               deltaG/RT);
        printf("                                : %14.5g (dimensionless/ln10) \n",
               deltaG/(RT * log(10.0)));

        i1 = HMW->speciesIndex("H+");
        i2 = HMW->speciesIndex("H2O(L)");
        j1 = HMW->speciesIndex("OH-");
        if (i1 == npos || i2 == npos || j1 == npos) {
            printf("problems\n");
            exit(-1);
        }
        deltaG = mu0[j1] + mu0[i1] - mu0[i2];
        printf(" OH-: H2O(L) - H+ -> OH-: %14.5g kJ/gmol \n",
               deltaG);
        printf("                                : %14.5g (dimensionless) \n",
               deltaG/RT);
        printf("                                : %14.5g (dimensionless/ln10) \n",
               deltaG/(RT * log(10.0)));



        printf("------------------------------------------------------------\n");


        delete HMW;
        HMW = 0;
        Cantera::appdelete();

        return retn;

    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return -1;
    }
}
