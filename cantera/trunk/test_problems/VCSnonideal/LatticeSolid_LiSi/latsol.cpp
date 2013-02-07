#include "cantera/thermo.h"
#include "cantera/thermo/MargulesVPSSTP.h"
#include "cantera/thermo/FixedChemPotSSTP.h"
#include "cantera/thermo/LatticeSolidPhase.h"
#include "cantera/equilibrium.h"
#include "cantera/thermo.h"

#include <stdio.h>

using namespace Cantera;

void testProblem()
{
    double T = 273.15 + 352.0;
    VCSnonideal::vcs_timing_print_lvl = 0;

    // Create the phases
    std::auto_ptr<ThermoPhase> LiSi_solid(newPhase("Li7Si3_ls.xml",
                                          "Li7Si3_and_Interstitials(S)"));
    FixedChemPotSSTP LiFixed("Li", -2.3E7);
    MargulesVPSSTP salt(1);

    // set states
    vector_fp x(salt.nSpecies(), 0);
    x[0] = 0.7;
    x[1] = 1.0 - x[0];

    salt.setState_TPX(T, OneAtm, &x[0]);
    LiSi_solid->setState_TP(T, OneAtm);

    int ee = LiSi_solid->nElements();
    printf("Number of elements = %d\n", ee);

    LiFixed.setState_TP(T, OneAtm);

    double um;
    LiFixed.getChemPotentials(&um);
    printf(" chem pot = %g\n", um);

    double volts = 1.635;     // has some Fe in it // test suite
    double dg_corr =  - volts * Faraday;
    printf("dg_corr = %g\n", dg_corr);

    double um_li_chempot = um + dg_corr;
    printf("um_li_chempot = %g\n", um_li_chempot);
    LiFixed.setChemicalPotential(um_li_chempot);

    Cantera::MultiPhase mmm;

    mmm.addPhase(&salt, 10.);
    mmm.addPhase(LiSi_solid.get(), 1.);
    mmm.addPhase(&LiFixed, 100.);

    int solver = 2;
    int printLvl = 10;
    int estimateEquil = 0;

    vcs_equilibrate(mmm, "TP", estimateEquil, printLvl, solver);
    std::cout << mmm << std::endl;

    Cantera::appdelete();
}

int main(int argc, char** argv)
{
    try {
        testProblem();
        return 0;
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return -1;
    }
}
