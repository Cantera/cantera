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
    int32_t sol0 = sol_newSolution("h2o2.yaml", "", "none");
    int32_t thermo = sol_thermo(sol0);
    int32_t mix = mix_new();
    ASSERT_GE(mix, 0);
    int32_t ret = mix_addPhase(mix, thermo, 1.);
    ASSERT_EQ(ret, 0);

    ASSERT_EQ(mix_nPhases(mix), 1);
    ASSERT_EQ(mix_nElements(mix), 4);
}
