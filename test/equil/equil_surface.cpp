// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/equil/MultiPhase.h"
#include "cantera/equil/MultiPhaseEquil.h"
#include "cantera/equil/vcs_MultiPhaseEquil.h"

using namespace Cantera;

TEST(SurfaceEquilibrium, directNonTPDriversReject)
{
    auto surface = newThermo("diamond.yaml", "diamond_100");
    MultiPhase mix;
    mix.addPhase(surface, 1e-5);
    mix.init();
    vector<double> initial(mix.nSpecies());
    mix.getMoles(initial);
    // Construction must succeed; only the requested non-TP modes should fail.
    vcs_MultiPhaseEquil solver(&mix, 0);
    EXPECT_THROW(solver.equilibrate_HP(mix.enthalpy(), HP, 300, 2000), CanteraError);
    EXPECT_THROW(solver.equilibrate_HP(mix.IntEnergy(), UP, 300, 2000), CanteraError);
    EXPECT_THROW(solver.equilibrate_SP(mix.entropy(), 300, 2000), CanteraError);
    for (int mode : {TV, HV, UV, SV}) {
        EXPECT_THROW(solver.equilibrate_TV(mode, 1200), CanteraError);
    }
    vector<double> final(mix.nSpecies());
    mix.getMoles(final);
    EXPECT_EQ(initial, final);
    EXPECT_DOUBLE_EQ(mix.temperature(), 1200);
}

TEST(SurfaceEquilibrium, directGibbsRejectsSurface)
{
    auto surface = newThermo("diamond.yaml", "diamond_100");
    MultiPhase mix;
    mix.addPhase(surface, 1e-5);
    mix.init();
    EXPECT_THROW(MultiPhaseEquil solver(&mix), CanteraError);
}

TEST(SurfaceEquilibrium, surfaceAddsNoVolume)
{
    auto gas = newThermo("diamond.yaml", "gas");
    auto surface = newThermo("diamond.yaml", "diamond_100");
    gas->setState_TP(1200, OneAtm);
    surface->setState_TP(1200, OneAtm);
    MultiPhase mix;
    mix.addPhase(surface, 1e-5);
    mix.addPhase(gas, 0.4);
    mix.init();
    EXPECT_NEAR(mix.volume(), 0.4 * gas->molarVolume(), 1e-12);
}
