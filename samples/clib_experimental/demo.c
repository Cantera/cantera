/*
 * Experimental CLib Demo
 * ======================
 *
 * This program illustrates using Cantera's experimental C-library interface to compute
 * thermodynamic, kinetic, and transport properties of a gas mixture. In addition,
 * a simple reactor network simulation is illustrated.
 *
 * .. tags:: C, tutorial, equilibrium, thermodynamics, kinetics, transport,
 *           reactor network
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/clib_experimental/ct3.h"
#include "cantera/clib_experimental/ctkin3.h"
#include "cantera/clib_experimental/ctrxn3.h"
#include "cantera/clib_experimental/ctsol3.h"
#include "cantera/clib_experimental/ctthermo3.h"
#include "cantera/clib_experimental/cttrans3.h"
#include "cantera/clib_experimental/ctreactor3.h"
#include "cantera/clib_experimental/ctreactornet3.h"

#include <stdio.h>

// sphinx_gallery_start_ignore
// The following header files are not used by this example, but are nevertheless added
// here to ensure C-compatibility of Cantera's experimental clib includes in continuous
// testing.
#include "cantera/clib_experimental/ctconnector3.h"
#include "cantera/clib_experimental/ctdomain3.h"
#include "cantera/clib_experimental/ctfunc3.h"
#include "cantera/clib_experimental/ctmix3.h"
#include "cantera/clib_experimental/ctonedim3.h"
#include "cantera/clib_experimental/ctrdiag3.h"

// sphinx_gallery_end_ignore

void exit_with_error()
{
    int buflen = 0;
    char* output_buf = 0;
    buflen = ct3_getCanteraError(buflen, output_buf) + 1;
    output_buf = malloc(sizeof(char) * buflen);
    ct3_getCanteraError(buflen, output_buf);
    printf("%s", output_buf);
    free(output_buf);
    exit(1);
}

int main(int argc, char** argv)
{
    ct3_make_deprecation_warnings_fatal(); // throw errors for deprecated functions

    int sol = sol3_newSolution("gri30.yaml", "gri30", "default");
    // In principle, one ought to check for errors after every Cantera call. But this
    // is the only one likely to occur in part of this example, due to the input file
    // not being found.
    if (sol < 0) {
        exit_with_error();
    }
    int thermo = sol3_thermo(sol);

    thermo3_setTemperature(thermo, 500);
    thermo3_setPressure(thermo, 5 * 101325);
    thermo3_setMoleFractionsByName(thermo, "CH4:1.0, O2:2.0, N2:7.52");
    thermo3_equilibrate(thermo, "HP", "auto", 1e-9, 50000, 1000, 0);
    thermo3_print(thermo, 1, 0);

    int kin = sol3_kinetics(sol);
    size_t nr = kin3_nReactions(kin);
    double T = thermo3_temperature(thermo);
    thermo3_setTemperature(thermo, T - 200);

    char buf[1000];
    double ropf[325];
    printf("\n                   Reaction           Forward ROP\n");
    kin3_getFwdRatesOfProgress(kin, 325, ropf);
    for (size_t n = 0; n < nr; n++) {
        int rxn = kin3_reaction(kin, n);
        rxn3_equation(rxn, 1000, buf);
        printf("%35s   %8.6e\n", buf, ropf[n]);
        rxn3_del(rxn);
    }

    int tran = sol3_transport(sol);
    size_t nsp = thermo3_nSpecies(thermo);
    printf("\n  Species    Mix diff coeff\n");
    double dkm[53];
    trans3_getMixDiffCoeffs(tran, 53, dkm);
    for (size_t k = 0; k < nsp; k++) {
        thermo3_getSpeciesName(thermo, k, 1000, buf);
        printf("%10s   %8.6e\n", buf, dkm[k]);
    }

    thermo3_setMoleFractionsByName(thermo, "CH4:1.0, O2:2.0, N2:7.52");
    thermo3_setTemperature(thermo, 1050);
    thermo3_setPressure(thermo, 5 * 101325);
    thermo3_print(thermo, 1, 1e-6);

    printf("\ntime       Temperature\n");
    int reactor = reactor3_new("IdealGasReactor", sol, "test");
    int reactors[1];
    reactors[0] = reactor;
    int net = reactornet3_new(1, &reactors[0]);

    double t = 0.0;
    int ret = 0;
    while (t < 0.1 && ret == 0) {
        double T = reactor3_temperature(reactor);
        t = reactornet3_time(net);
        printf("%.2e   %.3f\n", t, T);
        ret = reactornet3_advance(net, t + 5e-3);
    }

    ct3_appdelete();
    return 0;
}
