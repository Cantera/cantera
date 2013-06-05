#include "cantera/thermo.h"
#include "cantera/transport.h"
#include "cantera/thermo/HMWSoln.h"
#include "cantera/transport/SimpleTransport.h"

#include <memory>
#include <cstdio>

using namespace std;
using namespace Cantera;

int main(int argc, char** argv)
{
    try {
        int log_level = 3;

        HMWSoln HMW("HMW_NaCl_pdss.xml", "NaCl_electrolyte");

        auto_ptr<Transport> tran(newDefaultTransportMgr(&HMW, log_level));

        SimpleTransport& tranSimple = dynamic_cast<SimpleTransport&>(*tran.get());
        size_t nsp = HMW.nSpecies();

        HMW.setState_TP(30+273.13, OneAtm);

        double visc = tranSimple.viscosity();
        printf("visc = %g\n", visc);

        vector_fp x(nsp, 0.0);

        tranSimple.getSpeciesViscosities(&x[0]);
        for (size_t k = 0; k < nsp; k++) {
            printf("sp visc (%s) = %g\n", HMW.speciesName(k).c_str(), x[k]);
        }

        double cond = tranSimple.thermalConductivity();
        printf("cond = %g\n", cond);

        tranSimple.getMixDiffCoeffs(&x[0]);
        for (size_t k = 0; k < nsp; k++) {
            printf("sp diff (%s) = %g\n", HMW.speciesName(k).c_str(), x[k]);
        }

        tranSimple.getMobilities(&x[0]);
        for (size_t k = 0; k < nsp; k++) {
            printf("Mobility (%s) = %g\n", HMW.speciesName(k).c_str(), x[k]);
        }

        vector_fp gradX(nsp, 0.0);
        gradX[1] = 1.0;
        double gradT = 0.0;

        tranSimple.getSpeciesFluxes(1, &gradT, 5, &gradX[0], 5, &x[0]);
        for (size_t k = 0; k < nsp; k++) {
            string spName = HMW.speciesName(k);
            printf("SpeciesFlux (%s) = %g\n", spName.c_str(), x[k]);
        }

        gradX[1] = 0.0;
        double gradV = 1.0;

        tranSimple.set_Grad_T(&gradT);
        tranSimple.set_Grad_V(&gradV);
        tranSimple.set_Grad_X(&gradX[0]);

        tranSimple.getSpeciesFluxesExt(5, &x[0]);
        for (size_t k = 0; k < nsp; k++) {
            printf("SpeciesFlux (%s) = %g\n", HMW.speciesName(k).c_str(), x[k]);
        }

        Cantera::appdelete();
        return 0;

    } catch (CanteraError) {

        showErrors();
        return -1;
    }
}
