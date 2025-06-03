/*
 * CLib Demo
 * =========
 *
 * This program illustrates using Cantera's generated C-library interface to
 * compute thermodynamic, kinetic, and transport properties of a gas mixture. In
 * addition, a simple reactor network simulation is illustrated.
 *
 * .. tags:: C, tutorial, equilibrium, thermodynamics, kinetics, transport,
 *           reactor network
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera_clib/ct.h"
#include "cantera_clib/ctkin.h"
#include "cantera_clib/ctrxn.h"
#include "cantera_clib/ctsol.h"
#include "cantera_clib/ctthermo.h"
#include "cantera_clib/cttrans.h"
#include "cantera_clib/ctreactor.h"
#include "cantera_clib/ctreactornet.h"

#include <stdio.h>

// sphinx_gallery_start_ignore
// The following header files are not used by this example, but are nevertheless added
// here to ensure C-compatibility of Cantera's generated CLib includes in
// continuous testing.
#include "cantera_clib/ctconnector.h"
#include "cantera_clib/ctdomain.h"
#include "cantera_clib/ctfunc.h"
#include "cantera_clib/ctmix.h"
#include "cantera_clib/ctonedim.h"
#include "cantera_clib/ctrdiag.h"

// sphinx_gallery_end_ignore

void exit_with_error()
{
    int32_t buflen = 0;
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
    ct_make_deprecation_warnings_fatal(); // throw errors for deprecated functions

    int32_t sol = sol_newSolution("gri30.yaml", "gri30", "default");
    // In principle, one ought to check for errors after every Cantera call. But this
    // is the only one likely to occur in part of this example, due to the input file
    // not being found.
    if (sol < 0) {
        exit_with_error();
    }
    int32_t thermo = sol_thermo(sol);

    thermo_setTemperature(thermo, 500);
    thermo_setPressure(thermo, 5 * 101325);
    thermo_setMoleFractionsByName(thermo, "CH4:1.0, O2:2.0, N2:7.52");
    thermo_equilibrate(thermo, "HP", "auto", 1e-9, 50000, 1000, 0);
    thermo_print(thermo, 1, 0);

    int32_t kin = sol_kinetics(sol);
    int32_t nr = kin_nReactions(kin);
    double T = thermo_temperature(thermo);
    thermo_setTemperature(thermo, T - 200);

    char buf[1000];
    double ropf[325];
    printf("\n                   Reaction           Forward ROP\n");
    kin_getFwdRatesOfProgress(kin, 325, ropf);
    for (int32_t n = 0; n < nr; n++) {
        int32_t rxn = kin_reaction(kin, n);
        rxn_equation(rxn, 1000, buf);
        printf("%35s   %8.6e\n", buf, ropf[n]);
        rxn_del(rxn);
    }

    int32_t tran = sol_transport(sol);
    int32_t nsp = thermo_nSpecies(thermo);
    printf("\n  Species    Mix diff coeff\n");
    double dkm[53];
    trans_getMixDiffCoeffs(tran, 53, dkm);
    for (int32_t k = 0; k < nsp; k++) {
        thermo_speciesName(thermo, k, 1000, buf);
        printf("%10s   %8.6e\n", buf, dkm[k]);
    }

    thermo_setMoleFractionsByName(thermo, "CH4:1.0, O2:2.0, N2:7.52");
    thermo_setTemperature(thermo, 1050);
    thermo_setPressure(thermo, 5 * 101325);
    thermo_print(thermo, 1, 1e-6);

    printf("\ntime       Temperature\n");
    int32_t reactor = reactor_new("IdealGasReactor", sol, "test");
    int32_t reactors[1];
    reactors[0] = reactor;
    int32_t net = reactornet_new(1, &reactors[0]);

    double t = 0.0;
    int32_t ret = 0;
    while (t < 0.1 && ret == 0) {
        double T = reactor_temperature(reactor);
        t = reactornet_time(net);
        printf("%.2e   %.3f\n", t, T);
        ret = reactornet_advance(net, t + 5e-3);
    }

    ct_appdelete();
    return 0;
}
