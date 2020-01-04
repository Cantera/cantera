#include "cantera/thermo.h"
#include "cantera/thermo/MargulesVPSSTP.h"
#include "cantera/thermo/FixedChemPotSSTP.h"
#include "cantera/thermo/LatticeSolidPhase.h"
#include "cantera/thermo/ConstCpPoly.h"
#include "cantera/equil/vcs_MultiPhaseEquil.h"

#include <stdio.h>

using namespace Cantera;

void testProblem(int printLvl)
{
    double T = 273.15 + 352.0;
    VCS_SOLVE::disableTiming();

    // Create the phases
    std::unique_ptr<ThermoPhase> LiSi_solid(newPhase("Li7Si3_ls.xml",
                                                     "Li7Si3_and_Interstitials(S)"));
    std::unique_ptr<ThermoPhase> Li_liq(newPhase("Li_Liquid.xml", "Li(L)"));
    std::unique_ptr<ThermoPhase> LiFixed(newPhase("Li7Si3_ls.xml", "LiFixed"));
    MargulesVPSSTP salt("LiKCl_liquid.xml", "MoltenSalt_electrolyte");

    // set states
    vector_fp x(salt.nSpecies(), 0);
    x[0] = 0.7;
    x[1] = 1.0 - x[0];

    salt.setState_TPX(T, OneAtm, &x[0]);
    LiSi_solid->setState_TP(T, OneAtm);

    int ee = static_cast<int>(LiSi_solid->nElements());
    printf("Number of elements = %d\n", ee);

    LiFixed->setState_TP(T, OneAtm);

    double um[20];
    LiFixed->getChemPotentials(um);
    printf(" chem pot = %g\n", um[0]);

    double volts = 1.635;     // has some Fe in it // test suite
    double dg_corr =  - volts * Faraday;
    printf("dg_corr = %g\n", dg_corr);

    Li_liq->getChemPotentials(um);
    double um_li_chempot = um[0] + dg_corr;
    printf("um_li_chempot = %g\n", um_li_chempot);
    auto LiFixed_species = LiFixed->species(0);
    auto spthermo = std::dynamic_pointer_cast<ConstCpPoly>(LiFixed_species->thermo);
    spthermo->setParameters(298.15, um_li_chempot, 0, 0);
    LiFixed->modifySpecies(0, LiFixed_species);

    MultiPhase mmm;

    mmm.addPhase(&salt, 10.);
    mmm.addPhase(LiSi_solid.get(), 1.);
    mmm.addPhase(LiFixed.get(), 100.);

    int estimateEquil = 0;

    mmm.equilibrate("TP", "vcs", 1e-9, 50000, 100, estimateEquil, printLvl);
    std::cout << mmm << std::endl;

    appdelete();
}

int main(int argc, char** argv)
{
    int printLvl = 3;
    if (argc > 1) {
        printLvl = atoi(argv[1]);
    }
    try {
        testProblem(printLvl);
        return 0;
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return -1;
    }
}
