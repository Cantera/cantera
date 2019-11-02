// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/IdealGasPhase.h"

#include <iostream>

using namespace std;
using namespace Cantera;

int main(int argc, char** argv)
{
#if defined(_MSC_VER) && _MSC_VER < 1900
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    try {
        auto sol = newSolution("air_below6000K.cti", "air_below6000K");
        auto gas = sol->thermo();

        vector_fp IndVar2(6, 0.0);
        IndVar2[0] = 1.5E5;
        IndVar2[1] = 3.0E5;
        IndVar2[2] = 9.0E5;
        IndVar2[3] = 2.7E6;
        IndVar2[4] = 6.7E6;
        IndVar2[5] = 1.0E7;

        vector_fp IndVar1(7, 0.0);
        IndVar1[0] = 1.0E-8;
        IndVar1[1] = 1.0E-7;
        IndVar1[2] = 1.0E-6;
        IndVar1[3] = 1.0E-5;
        IndVar1[4] = 1.0E-4;
        IndVar1[5] = 1.0E-3;
        IndVar1[6] = 1.0E-2;
        int nj = 6;
        int ni = 7;
        FILE* FF = fopen("table.csv","w");
        size_t kk = gas->nSpecies();
        vector_fp Xmol(kk, 0.0);
        const std::vector<string> &snames = gas->speciesNames();
        fprintf(FF,"Temperature,  Pressure,");
        for (size_t k = 0; k < kk; k++) {
            fprintf(FF, "%11s", snames[k].c_str());
            if (k < kk-1) {
                fprintf(FF,",");
            }
        }
        fprintf(FF,"\n");

        for (int j=0; j<nj; j++) {
            for (int i=0; i<ni; i++) {
                double offset = -301471.39;
                gas->setState_UV(IndVar2[j]+offset,1.0/IndVar1[i]);
                double tkelvin = gas->temperature();
                double pres = gas->pressure();
                printf("Initial T = %g, pres = %g atm\n", tkelvin, pres/OneAtm);

                gas->equilibrate("UV", "element_potential", 1e-12);

                tkelvin = gas->temperature();
                pres = gas->pressure();
                printf("Final T = %g, pres = %g atm\n", tkelvin, pres/OneAtm);
                cout << "enthalpy = " << gas->enthalpy_mass() << endl;
                cout << "entropy = " << gas->entropy_mass() << endl;
                cout << "Gibbs function = " << gas->gibbs_mass() << endl;
                cout << "heat capacity c_p = " << gas->cp_mass() << endl;
                cout << "heat capacity c_v = " << gas->cv_mass() << endl << endl;

                gas->getMoleFractions(Xmol.data());
                fprintf(FF,"%10.4g, %10.4g,", tkelvin, pres);
                for (size_t k = 0; k < kk; k++) {
                    if (fabs(Xmol[k]) < 1.0E-130) {
                        fprintf(FF," %10.4g", 0.0);
                    } else {
                        fprintf(FF," %10.4g", Xmol[k]);
                    }
                    if (k < kk-1) {
                        fprintf(FF,",");
                    }
                }
                fprintf(FF,"\n");
            }
        }
        fclose(FF);
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return -1;
    }
    return 0;
}
