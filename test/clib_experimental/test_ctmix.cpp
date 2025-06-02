#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <fstream>

#include "cantera/core.h"
#include "cantera_clib/ct.h"
#include "cantera_clib/ctsol.h"
#include "cantera_clib/ctmix.h"

using namespace Cantera;


TEST(ctmix, new)
{
    int sol0 = sol_newSolution("h2o2.yaml", "", "none");
    int thermo = sol_thermo(sol0);
    int mix = mix_new();
    ASSERT_GE(mix, 0);
    int ret = mix_addPhase(mix, thermo, 1.);
    ASSERT_EQ(ret, 0);

    ASSERT_EQ(mix_nPhases(mix), 1);
    ASSERT_EQ(mix_nElements(mix), 4);
}
