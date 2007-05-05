// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ThermoPhase.h"
#include <stdio.h>

using namespace std;

namespace Cantera {

    /**
     * Format a summary of the mixture state for output.
     */           
    std::string report(const ThermoPhase& th, bool show_thermo) {

        try {
        char p[200];
        string s = "";

        sprintf(p, "\n       temperature    %12.6g  K\n", th.temperature());
        s += p;
        sprintf(p, "          pressure    %12.6g  Pa\n", th.pressure());
        s += p;
        sprintf(p, "           density    %12.6g  kg/m^3\n", th.density());
        s += p;
        sprintf(p, "  mean mol. weight    %12.6g  amu\n", th.meanMolecularWeight());
        s += p;

        if (show_thermo) {
        sprintf(p, "\n");
        s += p;
        sprintf(p, "                          1 kg            1 kmol\n");
        s += p;
        sprintf(p, "                       -----------      ------------\n");
        s += p;
        sprintf(p, "          enthalpy    %12.6g     %12.4g     J\n", 
            th.enthalpy_mass(), th.enthalpy_mole());
        s += p;
        sprintf(p, "   internal energy    %12.6g     %12.4g     J\n", 
            th.intEnergy_mass(), th.intEnergy_mole());
        s += p;
        sprintf(p, "           entropy    %12.6g     %12.4g     J/K\n", 
            th.entropy_mass(), th.entropy_mole());
        s += p;
        sprintf(p, "    Gibbs function    %12.6g     %12.4g     J\n", 
            th.gibbs_mass(), th.gibbs_mole());
        s += p;
        sprintf(p, " heat capacity c_p    %12.6g     %12.4g     J/K\n", 
            th.cp_mass(), th.cp_mole());
        s += p;
        sprintf(p, " heat capacity c_v    %12.6g     %12.4g     J/K\n", 
            th.cv_mass(), th.cv_mole());
        s += p;
        }

        int kk = th.nSpecies();
        array_fp x(kk);
        array_fp y(kk);
        th.getMoleFractions(DATA_PTR(x));
        th.getMassFractions(DATA_PTR(y));

        int k;

        sprintf(p, "\n                           X                 Y   \n");
        s += p;
        sprintf(p, "                     -------------     ------------\n");
        s += p;
        for (k = 0; k < kk; k++) {
                sprintf(p, "%18s   %12.6e     %12.6e\n", 
                    th.speciesName(k).c_str(), x[k], y[k]);
                s += p;
        }
        return s;
        }
        catch (CanteraError) {
            return std::string("<error>");
        }
    }

    /**
     * Format a composition list for output.
     */           
    std::string formatCompList(const Phase& mix, int xyc) {

        const doublereal Threshold = 1.e-20;

        char p[200];
        std::string s = "";
        int kk = mix.nSpecies();
        array_fp zz(kk);
        switch (xyc) {
        case 0: mix.getMoleFractions(DATA_PTR(zz)); break;
        case 1: mix.getMassFractions(DATA_PTR(zz)); break;
        case 2: mix.getConcentrations(DATA_PTR(zz)); break;
        default: return "error: xyc must be 0, 1, or 2";
        }

        doublereal z;
        int k;
        for (k = 0; k < kk; k++) {
            z = fabs(zz[k]);
            if (z < Threshold) zz[k] = 0.0;
        }

        for (k = 0; k < kk; k++) {
            sprintf(p, "%18s\t %12.6e\n", mix.speciesName(k).c_str(), 
                zz[k]);
            s += p;
        }
        return s;
    }

}

