// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_GTEST_UTILS_H
#define CT_GTEST_UTILS_H

// Wrap a test name in SLOW_TEST to disable it if the "skip_slow_tests=y" option
// is passed to SCons. For example:
//     TEST_F(GriMatrix, SLOW_TEST(VcsNonideal_CH4_O2))

#ifdef CT_SKIP_SLOW
#define SLOW_TEST(x) DISABLED_ ## x
#else
#define SLOW_TEST(x) x
#endif

#endif
