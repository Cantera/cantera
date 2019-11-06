/**
 *  @file mixGasTransport.cpp
 *       test problem for mixture transport
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

//  Example
//
// Test case for mixture transport in a gas
// The basic idea is to set up a gradient of some kind.
// Then the resulting transport coefficients out.
// Essentially all of the interface routines should be
// exercised and the results dumped out.
//
// A blessed solution test will make sure that the actual
// solution doesn't change as a function of time or
// further development.

// perhaps, later, an analytical solution could be added

#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/transport/MixTransport.h"

#include <iostream>

using namespace std;
using namespace Cantera;

int main(int argc, char** argv)
{
#if defined(_MSC_VER) && _MSC_VER < 1900
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    string infile = "diamond.xml";

    try {
        auto sol = newSolution("gri30.yaml", "gri30", "Mix");
        auto gas = sol->thermo();
        size_t nsp = gas->nSpecies();
        double pres = 1.0E5;
        vector_fp Xset(nsp, 0.0);
        Xset[0] =  0.269205 ;
        Xset[1] =  0.000107082;
        Xset[2] =  1.36377e-09 ;
        Xset[3] =  4.35475e-10;
        Xset[4] =  4.34036e-06 ;
        Xset[5] =  0.192249;
        Xset[6] =  3.59356e-13;
        Xset[7] =  2.78061e-12 ;
        Xset[8] =  4.7406e-18   ;
        Xset[9] =  4.12955e-17 ;
        Xset[10] = 2.58549e-14  ;
        Xset[11] = 8.96502e-16 ;
        Xset[12] = 6.09056e-11   ;
        Xset[13] = 7.56752e-09  ;
        Xset[14] = 0.192253;
        Xset[15] = 0.0385036;
        Xset[16] = 1.49596e-08   ;
        Xset[17] = 2.22378e-08     ;
        Xset[18] =   1.43096e-13   ;
        Xset[19] =   1.45312e-15 ;
        Xset[20] =  1.96948e-12 ;
        Xset[21] =   8.41937e-19;
        Xset[22] =  3.18852e-13 ;
        Xset[23] =  7.93625e-18 ;
        Xset[24] = 3.20653e-15  ;
        Xset[25] = 1.15149e-19 ;
        Xset[26] = 1.61189e-18  ;
        Xset[27] =   1.4719e-15 ;
        Xset[28] =  5.24728e-13 ;
        Xset[29] = 6.90582e-17  ;
        Xset[30] = 6.37248e-12   ;
        Xset[31] =5.93728e-11   ;
        Xset[32] =   2.71219e-09  ;
        Xset[33] =    2.66645e-06 ;
        Xset[34] =   6.57142e-11 ;
        Xset[35] =   9.52453e-08 ;
        Xset[36] =   1.26006e-14;
        Xset[37] =   3.49802e-12;
        Xset[38] =    1.19232e-11 ;
        Xset[39] =    7.17782e-13    ;
        Xset[40] =    1.85347e-07   ;
        Xset[41] =    8.25325e-14   ;
        Xset[42] =  5.00914e-20 ;
        Xset[43] = 1.54407e-16 ;
        Xset[44] =3.07176e-11 ;
        Xset[45] =4.93198e-08 ;
        Xset[46] =4.84792e-12 ;
        Xset[47] = 0.307675  ;
        Xset[48] =0;
        Xset[49] =6.21649e-29;
        Xset[50] = 8.42393e-28 ;
        Xset[51] = 6.77865e-18;
        Xset[52] = 2.19225e-16;
        double T1 = 1500.;

        double sum = 0.0;
        for (size_t k = 0; k < nsp; k++) {
            sum += Xset[k];
        }
        for (size_t k = 0; k < nsp; k++) {
            Xset[k] /= sum;
        }

        vector_fp X2set(nsp, 0.0);
        X2set[0]  = 0.25 ;
        X2set[5]  = 0.17;
        X2set[14] = 0.15;
        X2set[15] = 0.05;
        X2set[47] =  0.38 ;
        double T2 = 1200.;

        double dist = 0.1;

        vector_fp X3set(nsp, 0.0);
        X3set[0]  = 0.27 ;
        X3set[5]  = 0.15;
        X3set[14] = 0.18;
        X3set[15] = 0.06;
        X3set[47] = 0.36 ;
        double T3 = 1400.;

        vector_fp grad_T(3, 0.0);
        Array2D grad_X(nsp, 2, 0.0);

        for (size_t k = 0; k < nsp; k++) {
            grad_X(k,0) = (X2set[k] - Xset[k])/dist;
            grad_X(k,1) = (X3set[k] - Xset[k])/dist;
        }

        grad_T[0] = (T2 - T1) / dist;
        grad_T[1] = (T3 - T1) / dist;

        auto tran = sol->transport();
        gas->setState_TPX(1500.0, pres, Xset.data());

        vector_fp mixDiffs(nsp, 0.0);

        tran->getMixDiffCoeffs(mixDiffs.data());
        printf(" Dump of the mixture Diffusivities:\n");
        for (size_t k = 0; k < nsp; k++) {
            string sss = gas->speciesName(k);
            printf("    %15s %13.5g\n", sss.c_str(), mixDiffs[k]);
        }

        vector_fp specVisc(nsp, 0.0);

        tran->getSpeciesViscosities(specVisc.data());
        printf(" Dump of the species viscosities:\n");
        for (size_t k = 0; k < nsp; k++) {
            string sss = gas->speciesName(k);
            printf("    %15s %13.5g\n", sss.c_str(), specVisc[k]);
        }

        vector_fp thermDiff(nsp, 0.0);
        tran->getThermalDiffCoeffs(thermDiff.data());
        printf(" Dump of the Thermal Diffusivities :\n");
        for (size_t k = 0; k < nsp; k++) {
            string sss = gas->speciesName(k);
            printf("    %15s %13.5g\n", sss.c_str(), thermDiff[k]);
        }

        printf("Viscosity and thermal Cond vs. T\n");
        for (size_t k = 0; k < 10; k++) {
            T1 = 400. + 100. * k;
            gas->setState_TPX(T1, pres, Xset.data());
            double visc = tran->viscosity();
            double cond = tran->thermalConductivity();
            printf("    %13.4g %13.4g %13.4g\n", T1, visc, cond);
        }

        gas->setState_TPX(T1, pres, Xset.data());

        Array2D Bdiff(nsp, nsp, 0.0);
        printf("Binary Diffusion Coefficients H2 vs species\n");

        tran->getBinaryDiffCoeffs(nsp, Bdiff.ptrColumn(0));
        for (size_t k = 0; k < nsp; k++) {
            string sss = gas->speciesName(k);
            printf(" H2 -   %15s %13.4g %13.4g\n", sss.c_str(), Bdiff(0,k), Bdiff(k,0));
        }


        vector_fp specMob(nsp, 0.0);

        tran->getMobilities(specMob.data());
        printf(" Dump of the species mobilities:\n");
        for (size_t k = 0; k < nsp; k++) {
            string sss = gas->speciesName(k);
            printf("    %15s %13.4g\n", sss.c_str(), specMob[k]);
        }

        Array2D fluxes(nsp, 2, 0.0);

        tran->getSpeciesFluxes(2, grad_T.data(), nsp,
                               grad_X.ptrColumn(0), nsp, fluxes.ptrColumn(0));
        printf(" Dump of the species fluxes:\n");
        double sum1 = 0.0;
        double sum2 = 0.0;
        double max1 = 0.0;
        double max2 = 0.0;
        for (size_t k = 0; k < nsp; k++) {
            string sss = gas->speciesName(k);
            printf("    %15s %13.4g %13.4g\n", sss.c_str(), fluxes(k,0), fluxes(k,1));
            sum1 += fluxes(k,0);
            if (fabs(fluxes(k,0)) > max1) {
                max1 = fabs(fluxes(k,0));
            }
            sum2 += fluxes(k,1);
            if (fabs(fluxes(k,1)) > max2) {
                max2 = fabs(fluxes(k,0));
            }
        }

        // Make sure roundoff error doesn't interfere with the printout.
        // these should be zero.
        if (fabs(sum1) * 1.0E14 > max1) {
            printf("sum in x direction = %13.4g\n", sum1);
        } else {
            printf("sum in x direction = 0\n");
        }
        if (fabs(sum2) * 1.0E14 > max2) {
            printf("sum in y direction = %13.4g\n", sum1);
        } else {
            printf("sum in y direction = 0\n");
        }


    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return -1;
    }

    return 0;
}
/***********************************************************/
