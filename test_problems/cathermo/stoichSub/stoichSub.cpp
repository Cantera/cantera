/**
 *
 *  @file HMW_graph_1.cpp
 */

#include "cantera/thermo/StoichSubstance.h"
#include "TemperatureTable.h"
#include <iostream>

using namespace std;
using namespace Cantera;

int main(int argc, char** argv)
{
    try {
        std::string iFile = (argc > 1) ? argv[1] : "NaCl_Solid.xml";
        std::string file_ID = iFile + "#NaCl(S)";
        XML_Node* xm = get_XML_NameID("phase", file_ID, 0);
        StoichSubstance* solid = new StoichSubstance(*xm);

        size_t nsp = solid->nSpecies();
        if (nsp != 1) {
            throw CanteraError("main","Should just be one species");
        }
        string sName;

        TemperatureTable TTable(8, true, 300, 100., 0, 0);

        /*
         * Set the Pressure
         */
        double pres = OneAtm;
        double T = 298.15;
        solid->setState_TP(T, pres);

        /*
         * ThermoUnknowns
         */
        double mu0_RT[20], mu[20], cp_r[20];
        double enth_RT[20];
        double entrop_RT[20], intE_RT[20];
        double mu_NaCl, enth_NaCl, entrop_NaCl;
        double cp_NaCl;
        /*
         * Create a Table of NaCl  Properties as a Function
         * of the Temperature
         */

        double RT = GasConstant * T;
        solid->getEnthalpy_RT(enth_RT);
        double enth_NaCl_298 = enth_RT[0] * RT * 1.0E-6;

        printf(" Data from http://webbook.nist.gov\n");
        printf("\n");

        printf("           T,    Pres,    molarGibbs0,    Enthalpy,      Entropy,         Cp  ,"
               "  -(G-H298)/T,     H-H298 ");
        printf("\n");

        printf("      Kelvin,    bars,       kJ/gmol,      kJ/gmol,      J/gmolK,     J/gmolK ,"
               "      J/gmolK,     J/gmol");
        printf("\n");

        for (size_t i = 0; i < TTable.NPoints; i++) {
            T = TTable.T[i];

            // GasConstant is in J/kmol
            RT = GasConstant * T;
            pres = OneAtm;

            solid->setState_TP(T, pres);
            /*
            * Get the Standard State DeltaH
            */
            solid->getGibbs_RT(mu0_RT);

            solid->getEnthalpy_RT(enth_RT);
            enth_NaCl = enth_RT[0] * RT * 1.0E-6;

            solid->getChemPotentials(mu);
            mu_NaCl = mu[0] * 1.0E-6;

            solid->getEntropy_R(entrop_RT);
            entrop_NaCl = entrop_RT[0] * GasConstant * 1.0E-3;

            solid->getIntEnergy_RT(intE_RT);

            solid->getCp_R(cp_r);
            cp_NaCl = cp_r[0] * GasConstant * 1.0E-3;

            double pbar = pres * 1.0E-5;

            printf("%10g, %10g, %12g, %12g, %12g, %12g, %12g, %12g",
                   T, pbar, mu_NaCl, enth_NaCl, entrop_NaCl, cp_NaCl, -1.0E3*(mu_NaCl-enth_NaCl_298)/T, enth_NaCl-enth_NaCl_298);
            printf("\n");
        }

        delete solid;
        solid = 0;
        Cantera::appdelete();
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        Cantera::appdelete();
        return -1;
    }
    return 0;
}
