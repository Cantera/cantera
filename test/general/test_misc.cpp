// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "cantera/base/global.h"
#include "cantera/base/Solution.h"

using namespace Cantera;
using ::testing::HasSubstr;

TEST(FatalError, stacktrace) {
    EXPECT_DEATH(std::abort(), "Stack trace");
}

TEST(CanteraError, stacktrace) {
    bool raised = false;
    try {
        newSolution("xyz567.yaml");
    } catch (CanteraError& err) {
        raised = true;
        // Check for information about an intermediate function call
        EXPECT_THAT(err.what(), testing::HasSubstr("AnyMap::fromYamlFile"));
    }
    EXPECT_TRUE(raised);
}
