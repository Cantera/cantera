/**
 *  @file HMW_graph_1.cpp
 */

#include "cantera/thermo.h"
#include "cantera/thermo/HMWSoln.h"

#include "TemperatureTable.h"

#include <cstdio>

using namespace std;
using namespace Cantera;

int main(int argc, char** argv)
{
    int retn = 0;
    size_t i;

    try {
        std::string iFile = (argc > 1) ? argv[1] : "HMW_NaCl.xml";
        double Enth0_RT[20], pmEnth[20], molarEnth;

        HMWSoln* HMW = new HMWSoln(iFile, "NaCl_electrolyte");


        /*
         * Load in and initialize the
         */
        Cantera::ThermoPhase* solid = newPhase("NaCl_Solid.xml","NaCl(S)");


        size_t nsp = HMW->nSpecies();
        double mf[100];
        double moll[100];
        for (i = 0; i < 100; i++) {
            mf[i] = 0.0;
            moll[i] = 0.0;
        }
        HMW->getMoleFractions(mf);
        string sName;

        TemperatureTable TTable(15, false, 273.15, 25., 0, 0);


        HMW->setState_TP(298.15, 1.01325E5);

        size_t i1 = HMW->speciesIndex("Na+");
        size_t i2 = HMW->speciesIndex("Cl-");
        for (i = 1; i < nsp; i++) {
            moll[i] = 0.0;
        }
        HMW->setMolalities(moll);

        double Is = 0.0;

        /*
         * Set the Pressure
         */
        double pres = OneAtm;

        /*
         * Fix the molality
         */
        Is = 6.146;
        moll[i1] = Is;
        moll[i2] = Is;
        HMW->setState_TPM(298.15, pres, moll);
        double Xmol[30];
        HMW->getMoleFractions(Xmol);

        /*
         * ThermoUnknowns
         */
        double T;
        double H0_NaCl = 0.0, H0_Naplus = 0.0, H0_Clminus = 0.0, Delta_H0s, H0_H2O = 0.0;
        double H_NaCl = 0.0, H_Naplus = 0.0, H_Clminus = 0.0, H_H2O = 0.0;
        double molarEnth0;

        printf("A_L/RT: Comparison to Pitzer's book, p. 99, can be made.\n");
        printf("        Agreement to 3-4 sig digits for Aphi and A_L/RT\n");
        printf("\n");

        printf("Delta_H0: Heat of Solution per mole of salt (standard states)\n");
        printf("          rxn for the ss heat of soln:     "
               "NaCl(s) -> Na+(aq) + Cl-(aq)\n");

        printf("\n");
        printf("Delta_Hs: Heat of Solution per mole of salt\n");
        printf("          rxn for heat of soln:     "
               " n1 H2O(l,pure) + n2 NaCl(s) -> n2 MX(aq) + n1 H2O(l) \n");
        printf("          Delta_Hs = (n1 h_H2O_bar + n2 h_MX_bar "
               "- n1 h_H2O_0 - n2 h_MX_0)/n2\n");
        printf("\n");
        printf("phiL:     phiL, calculated from the program, is checked\n");
        printf("          against analytical formula in L_standalone program.\n");
        printf("          (comparison against Eq. 12, Silvester and Pitzer)\n");

        /*
         * Create a Table of NaCl Enthalpy Properties as a Function
         * of the Temperature
         */
        printf("\n\n");
        printf("            T,          Pres,          Aphi,        A_L/RT,"
               "      Delta_H0,"
               "      Delta_Hs,             L,          phiL,"
               "   L_rel_molal,"
               "     MolarEnth,   MolarEnth0\n");
        printf("       Kelvin,           bar, sqrt(kg/gmol), sqrt(kg/gmol),"
               "   kJ/gmolSalt,"
               "   kJ/gmolSalt,   kJ/gmolSoln,   kJ/gmolSalt,"
               "   kJ/gmolSalt,       kJ/gmol,    kJ/gmol\n");
        for (i = 0; i < TTable.NPoints + 1; i++) {
            if (i == TTable.NPoints) {
                T = 323.15;
            } else {
                T = TTable.T[i];
            }
            /*
             * RT is in units of J/kmolK
             */
            double RT = GasConstant * T;

            /*
             * Make sure we are at the saturation pressure or above.
             */
            pres = std::max(HMW->satPressure(T), OneAtm);

            HMW->setState_TPM(T, pres, moll);

            solid->setState_TP(T, pres);

            /*
             * Get the Standard State DeltaH
             */
            solid->getEnthalpy_RT(Enth0_RT);
            H0_NaCl = Enth0_RT[0] * RT * 1.0E-6;
            HMW->getEnthalpy_RT(Enth0_RT);
            H0_H2O    = Enth0_RT[0] * RT * 1.0E-6;
            H0_Naplus = Enth0_RT[i1] * RT * 1.0E-6;
            H0_Clminus = Enth0_RT[i2] * RT * 1.0E-6;

            /*
             * Calculate the standard state heat of solution
             * for NaCl(s) -> Na+ + Cl-
             *   units: kJ/gmolSalt
             */
            Delta_H0s = H0_Naplus + H0_Clminus - H0_NaCl;


            solid->getPartialMolarEnthalpies(pmEnth);
            H_NaCl = pmEnth[0] * 1.0E-6;
            HMW->getPartialMolarEnthalpies(pmEnth);
            H_H2O     = pmEnth[0] * 1.0E-6;
            H_Naplus  = pmEnth[i1] * 1.0E-6;
            H_Clminus = pmEnth[i2] * 1.0E-6;
            molarEnth = HMW->enthalpy_mole() * 1.0E-6;

            double Delta_Hs = (Xmol[0] * H_H2O +
                               Xmol[i1] * H_Naplus +
                               Xmol[i2] * H_Clminus
                               - Xmol[0] * H0_H2O
                               - Xmol[i1] * H_NaCl);
            Delta_Hs /= Xmol[i1];

            /*
             * Calculate the relative enthalpy, L, from the
             * partial molar quantities. units kJ/gmolSolution
             */
            double L = (Xmol[0] * (H_H2O - H0_H2O) +
                        Xmol[i1] * (H_Naplus - H0_Naplus) +
                        Xmol[i2] * (H_Clminus - H0_Clminus));

            /*
             * Calculate the apparent relative molal enthalpy, phiL,
             * units of kJ/gmolSaltAdded
             */
            double phiL = L / Xmol[i1];
            double Aphi = HMW->A_Debye_TP() / 3.0;
            double AL = HMW->ADebye_L();
            double LrelMol = HMW->relative_molal_enthalpy() * 1.0E-6;

            for (size_t k = 0; k < nsp; k++) {
                Enth0_RT[k] *= RT * 1.0E-6;
            }

            molarEnth0 = 0.0;
            for (size_t k = 0; k < nsp; k++) {
                molarEnth0 += Xmol[k] * Enth0_RT[k];
            }

            if (i != TTable.NPoints+1) {
                printf("%13.4f, %13.4f, %13.4f, %13.4f, %13.4f, %13.4f, %13.4f, "
                       "%13.4f, %13.4f, %13.4f, %13.4f\n",
                       T, pres*1.0E-5,  Aphi, AL/RT, Delta_H0s, Delta_Hs,
                       L, phiL,
                       LrelMol,
                       molarEnth , molarEnth0);

            }

        }

        printf("Breakdown of Enthalpy Calculation at 323.15 K, 1atm:\n");

        printf(" Species     MoleFrac        Molal          H0      "
               "    partH     (partH - H0)\n");
        printf("  H2O(L)");
        printf("%13.4f %13.4f %13.4f %13.4f %13.4f\n", Xmol[0], moll[0], H0_H2O , H_H2O,  H_H2O-H0_H2O);
        printf("  Na+   ");
        printf("%13.4f %13.4f %13.4f %13.4f %13.4f\n", Xmol[i1], moll[i1],
               H0_Naplus , H_Naplus,  H_Naplus -H0_Naplus);
        printf("  Cl-   ");
        printf("%13.4f %13.4f %13.4f %13.4f %13.4f\n", Xmol[i2], moll[i2],
               H0_Clminus , H_Clminus,  H_Clminus - H0_Clminus);



        delete HMW;
        HMW = 0;
        delete solid;
        solid = 0;
        Cantera::appdelete();
        return retn;

    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        Cantera::appdelete();
        return -1;
    }
}
