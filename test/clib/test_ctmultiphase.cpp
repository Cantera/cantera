#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <fstream>

#include "cantera/core.h"
#include "cantera/clib/ct.h"
#include "cantera/clib/ctmultiphase.h"

using namespace Cantera;

TEST(ctmix, new)
{
    int thermo = thermo_newFromFile("gri30.yaml", "gri30");
    int mix = mix_new();
    ASSERT_GE(mix, 0);
    int ret = mix_addPhase(mix, thermo, 1.);
    ASSERT_EQ(ret, 0);
}
