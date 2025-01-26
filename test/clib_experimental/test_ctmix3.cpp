#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <fstream>

#include "cantera/core.h"
#include "cantera/clib/clib_defs.h"
#include "cantera/clib_experimental/ct3.h"
#include "cantera/clib_experimental/ctsol3.h"
#include "cantera/clib_experimental/ctmix3.h"

using namespace Cantera;


TEST(ctmix3, new)
{
    int sol0 = sol3_newSolution("h2o2.yaml", "", "none");
    int thermo = sol3_thermo(sol0);
    int mix = mix3_new();
    ASSERT_GE(mix, 0);
    int ret = mix3_addPhase(mix, thermo, 1.);
    ASSERT_EQ(ret, 0);

    ASSERT_EQ(mix3_nPhases(mix), 1);
    ASSERT_EQ(mix3_nElements(mix), 4);
}
