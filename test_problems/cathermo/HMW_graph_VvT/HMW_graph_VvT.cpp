/**
 *  @file HMW_graph_VvT
 */

#include "cantera/thermo.h"
#include "TemperatureTable.h"
#include "cantera/thermo/HMWSoln.h"

#include <cstdio>

using namespace std;
using namespace Cantera;

int main(int argc, char** argv)
{
    int retn = 0;
    size_t i;

    try {
        std::string iFile = (argc > 1) ? argv[1] : "HMW_NaCl.xml";
        double V0[20], pmV[20];

        HMWSoln* HMW = new HMWSoln(iFile, "NaCl_electrolyte");


        /*
         * Load in and initialize the
         */
        Cantera::ThermoPhase* solid = newPhase("NaCl_Solid.xml","NaCl(S)");


        size_t nsp = HMW->nSpecies();
        double mf[100];
        double moll[100];
        HMW->getMoleFractions(mf);
        string sName;

        TemperatureTable TTable(15, false, 273.15, 25., 0, 0);


        HMW->setState_TP(298.15, 1.01325E5);

        size_t i1 = HMW->speciesIndex("Na+");
        size_t i2 = HMW->speciesIndex("Cl-");
        for (i = 0; i < nsp; i++) {
            moll[i] = 0.0;
        }
        HMW->setMolalities(moll);

        /*
         * Set the Pressure
         */
        double pres = OneAtm;

        /*
         * Fix the molality
         */
        double Is = 6.146;
        moll[i1] = Is;
        moll[i2] = Is;
        HMW->setState_TPM(298.15, pres, moll);
        double Xmol[30];
        HMW->getMoleFractions(Xmol);
        double meanMW = HMW->meanMolecularWeight();

        /*
         * ThermoUnknowns
         */
        double T;

        double V0_NaCl = 0.0, V0_Naplus = 0.0, V0_Clminus = 0.0, Delta_V0s = 0.0, V0_H2O = 0.0;
        double V_NaCl = 0.0, V_Naplus = 0.0, V_Clminus = 0.0, V_H2O = 0.0;
        double molarV0;
        printf("A_V   : Comparison to Pitzer's book, p. 99, can be made.\n");
        printf("        Agreement to 3  sig digits \n");
        printf("\n");

        printf("Delta_V0: Heat Capacity of Solution per mole of salt (standard states)\n");
        printf("           rxn for the ss heat of soln:     "
               "NaCl(s) -> Na+(aq) + Cl-(aq)\n");

        printf("\n");
        printf("Delta_Vs: Delta volume of Solution per mole of salt\n");
        printf("          rxn for heat of soln:     "
               " n1 H2O(l,pure) + n2 NaCl(s) -> n2 MX(aq) + n1 H2O(l) \n");
        printf("          Delta_Hs = (n1 h_H2O_bar + n2 h_MX_bar "
               "- n1 h_H2O_0 - n2 h_MX_0)/n2\n");
        printf("\n");
        printf("phiV:     phiV, calculated from the program, is checked\n");
        printf("          against analytical formula in V_standalone program.\n");
        printf("          (comparison against Pitzer book, p. 97, eqn. 96)\n");

        /*
         * Create a Table of NaCl Enthalpy Properties as a Function
         * of the Temperature
         */
        printf("\n\n");
        printf("            T,          Pres,         Aphi,            A_V,"
               "      Delta_V0,"
               "      Delta_Vs,           Vex,          phiV,"
               "        MolarV,     MolarV0\n");
        printf("       Kelvin,         bar, sqrt(kg/gmol),sqrt(kg/gmol)cm3/gmol,"
               "cm**3/gmolSalt,"
               "cm**3/gmolSalt,cm**3/gmolSoln,cm**3/gmolSalt,"
               "cm**3/gmol,   cm**3/gmol\n");
        for (i = 0; i < TTable.NPoints + 1; i++) {
            if (i == TTable.NPoints) {
                T = 323.15;
            } else {
                T = TTable.T[i];
            }
            /*
             * Make sure we are at the saturation pressure or above.
             */
            pres = std::max(HMW->satPressure(T), OneAtm);

            HMW->setState_TPM(T, pres, moll);

            solid->setState_TP(T, pres);

            /*
             * Get the Standard State volumes m3/kmol
             */
            solid->getStandardVolumes(V0);
            V0_NaCl = V0[0];
            HMW->getStandardVolumes(V0);
            V0_H2O     = V0[0];
            V0_Naplus  = V0[i1];
            V0_Clminus = V0[i2];

            /*
             * Calculate the standard state volume change of solution
             * for NaCl(s) -> Na+ + Cl-
             *   units: m3 / kmol
             */
            Delta_V0s = V0_Naplus + V0_Clminus - V0_NaCl;

            double dd = solid->density();
            double MW_NaCl = solid->meanMolecularWeight();
            V_NaCl = MW_NaCl / dd;

            /*
             * Get the partial molar volumes
             */
            HMW->getPartialMolarVolumes(pmV);
            V_H2O     = pmV[0];
            V_Naplus  = pmV[i1];
            V_Clminus = pmV[i2];

            /*
             * Calculate the molar volume of solution
             */
            double dsoln = HMW->density();
            meanMW = HMW->meanMolecularWeight();
            double molarV = meanMW / dsoln;

            /*
             * Calculate the delta volume of solution for the reaction
             *                NaCl(s) -> Na+ + Cl-
             */
            double Delta_Vs = (Xmol[0] * V_H2O +
                               Xmol[i1] * V_Naplus +
                               Xmol[i2] * V_Clminus
                               - Xmol[0] * V0_H2O
                               - Xmol[i1] * V_NaCl);
            Delta_Vs /= Xmol[i1];


            /*
             * Calculate the apparent molar volume, J, from the
             * partial molar quantities, units m3/kmol
             */
            double Vex = (Xmol[0] * (V_H2O    - V0_H2O) +
                          Xmol[i1] * (V_Naplus - V0_Naplus) +
                          Xmol[i2] * (V_Clminus - V0_Clminus));

            /*
             * Calculate the apparent relative molal volume, phiV,
             * units of m3/kmol
             */
            double phiV = Vex / Xmol[i1];

            double Aphi = HMW->A_Debye_TP(T, pres) / 3.0;
            double Av = HMW->ADebye_V(T, pres) * 1.0E3;

            molarV0 = 0.0;
            for (size_t k = 0; k < nsp; k++) {
                molarV0 += Xmol[k] * V0[k];
            }

            if (i != TTable.NPoints+1) {
                printf("%13.4g, %13.4g, %13.4g, %13.4g, %13.4g, %13.4g, "
                       "%13.5g, %13.4g, %13.4g, %13.4g\n",
                       T, pres*1.0E-5,  Aphi, Av, Delta_V0s*1.0E3, Delta_Vs*1.0E3,
                       Vex*1.0E3, phiV*1.0E3, molarV*1.0E3 , molarV0*1.0E3);
            }

        }

        printf("Breakdown of Volume Calculation at 323.15 K, 1atm:\n");

        printf(" Species     MoleFrac        Molal          V0      "
               "    partV     (partV - V0)\n");
        printf("  H2O(L)");
        printf("%13.4g %13.4g %13.4g %13.4g %13.4g\n", Xmol[0], moll[0], V0_H2O*1.E3, V_H2O*1.E3,
               (V_H2O-V0_H2O)*1.E3);
        printf("  Na+   ");
        printf("%13.4g %13.4g %13.4g %13.4g %13.4g\n", Xmol[i1], moll[i1],
               V0_Naplus*1.E3 , V_Naplus*1.E3, (V_Naplus -V0_Naplus)*1.E3);
        printf("  Cl-   ");
        printf("%13.4g %13.4g %13.4g %13.4g %13.4g\n", Xmol[i2], moll[i2],
               V0_Clminus*1.E3, V_Clminus*1.E3, (V_Clminus - V0_Clminus)*1.E3);

        printf(" NaCl(s)");
        double dd = V_NaCl*1.E3 - V0_NaCl*1.E3;
        if (fabs(dd) < 1.0E-12) {
            dd = 0.0;
        }
        printf("%13.4g               %13.4g %13.4g %13.4g\n", 1.0,
               V0_NaCl*1.E3 , V_NaCl*1.E3,  dd);


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
