/*
 * CLib Demo
 * =========
 *
 * This program illustrates using Cantera's C-library interface to compute
 * thermodynamic, kinetic, and transport properties of a gas mixture. In addition,
 * a simple reactor network simulation is illustrated.
 *
 * .. tags:: C, tutorial, equilibrium, thermodynamics, kinetics, transport,
 *           reactor network
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/clib/ct.h"
#include "cantera/clib/ctreactor.h"

#include <stdio.h>

// sphinx_gallery_start_ignore
// The following header files are not used by this example, but are nevertheless added
// here to ensure C-compatibility of Cantera's clib includes in continuous testing.
#include "cantera/clib/ctfunc.h"
#include "cantera/clib/ctmultiphase.h"
#include "cantera/clib/ctonedim.h"
#include "cantera/clib/ctrpath.h"
#include "cantera/clib/ctsurf.h"
// sphinx_gallery_end_ignore

void exit_with_error()
{
    int buflen = 0;
    char* output_buf = 0;
    buflen = ct_getCanteraError(buflen, output_buf) + 1;
    output_buf = malloc(sizeof(char) * buflen);
    ct_getCanteraError(buflen, output_buf);
    printf("%s", output_buf);
    free(output_buf);
    exit(1);
}

int main(int argc, char** argv)
{
    ct_suppress_deprecation_warnings(0); // throw errors for deprecated functions

    int soln = soln_newSolution("gri30.yaml", "gri30", "default");
    // In principle, one ought to check for errors after every Cantera call. But this
    // is the only one likely to occur in part of this example, due to the input file
    // not being found.
    if (soln < 0) {
        exit_with_error();
    }
    int thermo = soln_thermo(soln);

    thermo_setTemperature(thermo, 500);
    thermo_setPressure(thermo, 5 * 101325);
    thermo_setMoleFractionsByName(thermo, "CH4:1.0, O2:2.0, N2:7.52");
    thermo_equilibrate(thermo, "HP", "auto", 1e-9, 50000, 1000, 0);
    thermo_print(thermo, 1, 0);

    int kin = soln_kinetics(soln);
    size_t nr = kin_nReactions(kin);
    double T = thermo_temperature(thermo);
    thermo_setTemperature(thermo, T - 200);

    char buf [1000];
    double ropf[325];
    printf("\n                   Reaction           Forward ROP\n");
    kin_getFwdRatesOfProgress(kin, 325, ropf);
    for (size_t n = 0; n < nr; n++) {
        kin_getReactionString(kin, n, 1000, buf);
        printf("%35s   %8.6e\n", buf, ropf[n]);
    }

    int tran = soln_transport(soln);
    size_t nsp = thermo_nSpecies(thermo);
    printf("\n  Species    Mix diff coeff\n");
    double dkm[53];
    trans_getMixDiffCoeffs(tran, 53, dkm);
    for (size_t k = 0; k < nsp; k++) {
        thermo_getSpeciesName(thermo, k, 1000, buf);
        printf("%10s   %8.6e\n", buf, dkm[k]);
    }

    thermo_setMoleFractionsByName(thermo, "CH4:1.0, O2:2.0, N2:7.52");
    thermo_setTemperature(thermo, 1050);
    thermo_setPressure(thermo, 5 * 101325);
    thermo_print(thermo, 1, 1e-6);

    printf("\ntime       Temperature\n");
    int reactor = reactor_new("IdealGasReactor", soln, "test");
    int net = reactornet_new();
    ct_suppress_deprecation_warnings(1);
    reactornet_addreactor(net, reactor);
    ct_suppress_deprecation_warnings(0);

    double t = 0.0;
    int ret = 0;
    while (t < 0.1 && ret == 0) {
        double T = reactor_temperature(reactor);
        t = reactornet_time(net);
        printf("%.2e   %.3f\n", t, T);
        ret = reactornet_advance(net, t + 5e-3);
    }

    ct_appdelete();
    return 0;
}
