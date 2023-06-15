// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "gtest/gtest.h"
#include "cantera/base/global.h"

using namespace Cantera;

TEST(FatalError, stacktrace) {
    EXPECT_DEATH(std::abort(), "Stack trace");
}
