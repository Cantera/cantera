// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "cantera/base/global.h"
#include "cantera/base/Solution.h"

#include <cstdlib>

using namespace Cantera;
using ::testing::HasSubstr;

namespace {

void setFatalDeprecationEnv(bool enable)
{
    // The value is ignored by Cantera; only the presence of the variable matters.
#ifdef _WIN32
    _putenv_s("CANTERA_FATAL_DEPRECATION_WARNINGS", enable ? "1" : "");
#else
    if (enable) {
        setenv("CANTERA_FATAL_DEPRECATION_WARNINGS", "1", 1);
    } else {
        unsetenv("CANTERA_FATAL_DEPRECATION_WARNINGS");
    }
#endif
}

} // namespace

// Verify that defining the CANTERA_FATAL_DEPRECATION_WARNINGS environment variable
// makes deprecation warnings fatal. The variable is read once when the Application
// singleton is constructed (see Application::Application), so appdelete() is used to
// force it to be reconstructed and the variable re-read. The configuration set in
// main() (see string_processing.cpp) is restored at the end for subsequent tests.
TEST(Deprecation, fatal_from_environment) {
    setFatalDeprecationEnv(true);
    appdelete();
    EXPECT_THROW(warn_deprecated("Deprecation::fatal_from_environment", "testing"),
                 CanteraError);

    // Without the variable defined, the same warning is non-fatal. Suppress it to
    // avoid polluting the test output.
    setFatalDeprecationEnv(false);
    appdelete();
    suppress_deprecation_warnings();
    EXPECT_NO_THROW(warn_deprecated("Deprecation::fatal_from_environment", "testing"));

    // Restore the configuration established in main() for any subsequent tests.
    appdelete();
    make_deprecation_warnings_fatal();
    addDataDirectory("test/data");
    addDataDirectory("data");
}

TEST(FatalError, stacktrace) {
    EXPECT_DEATH(std::abort(), "Stack trace");
}

TEST(CanteraError, stacktrace) {
    // MinGW requires an external dependency (libbacktrace) to give meaningful
    // stacktraces.
    #ifdef __MINGW32__
        GTEST_SKIP();
    #endif

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
