/**
 *  @file DH_graph_1
 */

#include "fileLog.h"
#include "cantera/thermo/DebyeHuckel.h"

#include <cstdio>

using namespace std;
using namespace Cantera;

int main(int argc, char** argv)
{
    int retn = 0;
    size_t i;
    string fName = "DH_graph_1.log";
    fileLog* fl = new fileLog(fName);
    try {
        std::string iFile = (argc > 1) ? argv[1] : "DH_NaCl.xml";
        setLogger(fl);

        DebyeHuckel* DH = new DebyeHuckel(iFile, "NaCl_electrolyte");

        size_t nsp = DH->nSpecies();
        double acMol[100];
        double mf[100];
        double moll[100];
        DH->getMoleFractions(mf);
        string sName;

        DH->setState_TP(298.15, 1.01325E5);

        size_t i1 = DH->speciesIndex("Na+");
        size_t i2 = DH->speciesIndex("Cl-");
        size_t i3 = DH->speciesIndex("H2O(L)");
        for (i = 1; i < nsp; i++) {
            moll[i] = 0.0;
        }
        DH->setMolalities(moll);
        double Itop = 10.;
        double Ibot = 0.0;
        double ISQRTtop = sqrt(Itop);
        double ISQRTbot = sqrt(Ibot);
        double ISQRT;
        double Is = 0.0;
        size_t its = 100;
        printf("              Is,     sqrtIs,     meanAc,"
               "  log10(meanAC),     acMol_Na+,"
               ",     acMol_Cl-,   ac_Water\n");
        for (i = 0; i < its; i++) {
            ISQRT = ISQRTtop*((double)i)/(its - 1.0)
                    + ISQRTbot*(1.0 - (double)i/(its - 1.0));
            Is = ISQRT * ISQRT;
            moll[i1] = Is;
            moll[i2] = Is;
            DH->setMolalities(moll);
            DH->getMolalityActivityCoefficients(acMol);
            double meanAC = sqrt(acMol[i1] * acMol[i2]);
            printf("%15g, %15g, %15g, %15g, %15g, %15g, %15g\n",
                   Is, ISQRT, meanAC, log10(meanAC),
                   acMol[i1], acMol[i2], acMol[i3]);
        }


        delete DH;
        DH = 0;
        /*
         * This delete the file logger amongst other things.
         */
        appdelete();

        return retn;

    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        delete fl;
        return -1;
    }
}
