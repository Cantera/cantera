// Include all clib headers to make sure all of them are C-compatible, even if
// we don't actually use all of them in this test.
#include "cantera/clib/ct.h"
#include "cantera/clib/ctfunc.h"
#include "cantera/clib/ctmultiphase.h"
#include "cantera/clib/ctonedim.h"
#include "cantera/clib/ctreactor.h"
#include "cantera/clib/ctrpath.h"
#include "cantera/clib/ctsurf.h"
#include "cantera/clib/ctxml.h"

#include <stdio.h>
#include <assert.h>

int main(int argc, char** argv)
{
    int ret;
    int thermo = thermo_newFromFile("gri30.xml", "gri30_mix");
    assert(thermo > 0);
    size_t nsp = thermo_nSpecies(thermo);
    assert(nsp == 53);

    ret = thermo_setTemperature(thermo, 500);
    assert(ret == 0);
    ret = thermo_setPressure(thermo, 5 * 101325);
    assert(ret == 0);
    ret = thermo_setMoleFractionsByName(thermo, "CH4:1.0, O2:2.0, N2:7.52");
    assert(ret == 0);

    ret = thermo_equilibrate(thermo, "HP", 0, 1e-9, 50000, 1000, 0);
    assert(ret == 0);
    double T = thermo_temperature(thermo);
    assert(T > 2200 && T < 2300);

    ret = thermo_print(thermo, 1, 0);
    assert(ret == 0);

    int kin = kin_newFromFile("gri30.xml", "gri30_mix", thermo, 0, 0, 0, 0);
    assert(kin > 0);

    size_t nr = kin_nReactions(kin);
    assert(nr == 325 );

    ret = thermo_setTemperature(thermo, T - 200);
    assert(ret == 0);

    char buf [1000];
    double ropf[325];
    printf("\n                   Reaction           Forward ROP\n");
    kin_getFwdRatesOfProgress(kin, 325, ropf);
    int n; // declare this here for C89 compatibility
    for (n = 0; n < nr; n++) {
        kin_getReactionString(kin, n, 1000, buf);
        printf("%35s   %8.6e\n", buf, ropf[n]);
    }

    printf("\n  Species    Mix diff coeff\n");
    int tran = trans_newDefault(thermo, 0);
    double dkm[53];
    trans_getMixDiffCoeffs(tran, 53, dkm);
    int k; // declare this here for C89 compatibility
    for (k = 0; k < nsp; k++) {
        thermo_getSpeciesName(thermo, k, 1000, buf);
        printf("%10s   %8.6e\n", buf, dkm[k]);
    }

    ret = thermo_setTemperature(thermo, 1050);
    assert(ret == 0);
    ret = thermo_setPressure(thermo, 5 * 101325);
    assert(ret == 0);
    ret = thermo_setMoleFractionsByName(thermo, "CH4:1.0, O2:2.0, N2:7.52");
    assert(ret == 0);

    printf("\ntime       Temperature\n");
    int reactor = reactor_new2("IdealGasReactor");
    int net = reactornet_new();
    ret = reactor_setThermoMgr(reactor, thermo);
    assert(ret == 0);
    ret = reactor_setKineticsMgr(reactor, kin);
    assert(ret == 0);
    ret = reactornet_addreactor(net, reactor);
    assert(ret == 0);

    double t = 0.0;
    while (t < 0.1 && ret == 0) {
        double T = reactor_temperature(reactor);
        t = reactornet_time(net);
        printf("%.2e   %.3f\n", t, T);
        ret = reactornet_advance(net, t + 5e-3);
        assert(ret == 0);
    }
    ct_appdelete();
    return 0;
}
