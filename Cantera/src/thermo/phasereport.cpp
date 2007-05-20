
// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ThermoPhase.h"
#include "PureFluidPhase.h"
#include <stdio.h>
#include "mix_defs.h"

using namespace std;

namespace Cantera {

    /**
     * Format a summary of the mixture state for output.
     */           
    string report(const ThermoPhase& th, bool show_thermo) {
        char p[200];
        string s = "";
        try {
            if (th.name() != "") {
                sprintf(p, " \n  %s:\n", th.name().c_str());
                s += p;
            }
        sprintf(p, " \n       temperature    %12.6g  K\n", th.temperature());
        s += p;
        sprintf(p, "          pressure    %12.6g  Pa\n", th.pressure());
        s += p;
        sprintf(p, "           density    %12.6g  kg/m^3\n", th.density());
        s += p;
        sprintf(p, "  mean mol. weight    %12.6g  amu\n", th.meanMolecularWeight());
        s += p;
#ifdef WITH_PURE_FLUIDS
        if (th.eosType() == cPureFluid) {
            double xx = ((PureFluidPhase*)(&th))->vaporFraction();
            // if (th.temperature() < th.critTemperature()) {
            sprintf(p, "    vapor fraction    %12.6g \n", 
                xx); //th.vaporFraction());
                s += p;
                //}
        }
#endif
        doublereal phi = th.electricPotential();
        if (phi != 0.0) {
            sprintf(p, "         potential    %12.6g  V\n", phi);
            s += p;
        }
        if (show_thermo) {
        sprintf(p, " \n");
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
        try {
            sprintf(p, " heat capacity c_v    %12.6g     %12.4g     J/K\n", 
                th.cv_mass(), th.cv_mole());
            s += p;
        }
        catch(CanteraError) {
            sprintf(p, " heat capacity c_v    <not implemented>       \n");
            s += p;
        }
        }

        int kk = th.nSpecies();
        array_fp x(kk);
        array_fp y(kk);
        array_fp mu(kk);
        th.getMoleFractions(&x[0]);
        th.getMassFractions(&y[0]);
        th.getChemPotentials(&mu[0]);
        doublereal rt = GasConstant * th.temperature(); 
        int k;
        //if (th.nSpecies() > 1) {

            if (show_thermo) {
                sprintf(p, " \n                           X     "
                    "            Y          Chem. Pot. / RT    \n");
                s += p;
                sprintf(p, "                     -------------     "
                    "------------     ------------\n");
                s += p;
                for (k = 0; k < kk; k++) {
                    if (x[k] > SmallNumber) {
                        sprintf(p, "%18s   %12.6g     %12.6g     %12.6g\n", 
                            th.speciesName(k).c_str(), x[k], y[k], mu[k]/rt);
                    }
                    else {
                        sprintf(p, "%18s   %12.6g     %12.6g     \n", 
                            th.speciesName(k).c_str(), x[k], y[k]);
                    }
                    s += p;
                }
            }
            else {
                sprintf(p, " \n                           X"
                    "Y\n");
                s += p;
                sprintf(p, "                     -------------"
                    "     ------------\n");
                s += p;
                for (k = 0; k < kk; k++) {
                    sprintf(p, "%18s   %12.6g     %12.6g\n", 
                        th.speciesName(k).c_str(), x[k], y[k]);
                    s += p;
                }
            }
        }
        //}
        catch (CanteraError) {
            ;
        }
        return s;
    }

    void writephase(const ThermoPhase& th, bool show_thermo) {
        string s = report(th, show_thermo);
        writelog(s+"\n");
    }

    /**
     * Format a composition list for output.
     */           
    string formatCompList(const Phase& mix, int xyc) {

        const doublereal Threshold = 1.e-20;

        char p[200];
        string s = "";
        int kk = mix.nSpecies();
        array_fp zz(kk);
        switch (xyc) {
        case 0: mix.getMoleFractions(&zz[0]); break;
        case 1: mix.getMassFractions(&zz[0]); break;
        case 2: mix.getConcentrations(&zz[0]); break;
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


