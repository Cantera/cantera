/**
 *  @file ISSPTester.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

//  Example
//
//  Read a mechanism and a thermodynamics file for the
//  class IdealSolidSolnPhase in order to test that it's
//  working correctly
//

#include "cantera/thermo/IdealSolnGasVPSS.h"
#include <iostream>
#include <cstdio>

using namespace std;
using namespace Cantera;

int main(int argc, char** argv)
{
#if defined(_MSC_VER) && _MSC_VER < 1900
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

    try {
        double Tkelvin = 1200.;
        IdealSolnGasVPSS issp("IdealSolidSolnPhaseExample.xml");
        issp.setState_TPX(Tkelvin, OneAtm,
                          "C2H2-graph:0.3, C-graph:0.6, H2-solute:0.1");
        double hm = issp.enthalpy_mole();
        printf("molar enthalpy   = %13.5g J kg-1\n", hm);

        double um = issp.intEnergy_mole();
        printf("molar intEnergy  = %13.5g J kg-1\n", um);


        double sm = issp.entropy_mole();
        printf("molar entropy    = %13.5g J kg-1 K-1\n", sm);

        double gm = issp.gibbs_mole();
        printf("molar gibbs      = %13.5g J kg-1\n", gm);

        double cpm = issp.cp_mole();
        printf("molar Cp         = %13.5g J kg-1 K-1\n", cpm);

        double dens = issp.density();
        printf("mixture density = %13.5g kg m-3\n", dens);

        double mdens = issp.molarDensity();
        printf("molar density = %13.5g kmol m-3\n", mdens);

        double mmw = issp.meanMolecularWeight();
        printf("mean molecular weight = %13.5g kg kmol-1\n", mmw);

        size_t n = issp.nSpecies();

        double HiSS[20], muiSS[20],SiSS[20], CpiSS[20], VoliSS[20];
        double RT = GasConstant * Tkelvin;
        issp.getStandardChemPotentials(muiSS);
        issp.getEnthalpy_RT(HiSS);
        issp.getEntropy_R(SiSS);
        issp.getCp_R(CpiSS);
        issp.getStandardVolumes(VoliSS);


        for (size_t i = 0; i < n; i++) {
            HiSS[i] *= RT;
            SiSS[i] *= RT;
            CpiSS[i] *= GasConstant;
        }

        printf(" Printout of standard state properties\n");
        printf("            Name         mu_i       H_i_SS   "
               "    S_i_SS      Cp_i_SS     Vol_i_SS\n");
        for (size_t i = 0; i < n; i++) {
            string sn = issp.speciesName(i);
            printf(" %15s %12.5g %12.5g %12.5g %12.5g %12.5g\n", sn.c_str(), muiSS[i],
                   HiSS[i], SiSS[i], CpiSS[i], VoliSS[i]);
        }



        double HiPM[20], mui[20],SiPM[20], CpiPM[20], VoliPM[20];

        issp.getChemPotentials(mui);
        issp.getPartialMolarEnthalpies(HiPM);
        issp.getPartialMolarEntropies(SiPM);
        issp.getPartialMolarCp(CpiPM);
        issp.getPartialMolarVolumes(VoliPM);
        printf(" Printout of Partial molar properties\n");
        printf("            Name         mu_i       H_i_PM   "
               "    S_i_PM      Cp_i_PM     Vol_i_PM\n");
        for (size_t i = 0; i < n; i++) {
            string sn = issp.speciesName(i);
            printf(" %15s %12.5g %12.5g %12.5g %12.5g %12.5g\n", sn.c_str(), mui[i],
                   HiPM[i], SiPM[i], CpiPM[i], VoliPM[i]);
        }
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return -1;
    }

    return 0;
}
/***********************************************************/
