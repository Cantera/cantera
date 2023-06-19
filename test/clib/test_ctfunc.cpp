#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <fstream>

#include "cantera/core.h"
#include "cantera/clib/ctfunc.h"
#include "cantera/numerics/Func1.h"

using namespace Cantera;

TEST(ctfunc, sin)
{
    double omega = 2.;
    int fcn = func_new(SinFuncType, 0, 1, &omega);
    ASSERT_GE(fcn, 0);
    ASSERT_EQ(func_value(fcn, 0.), 0.);
    ASSERT_EQ(func_value(fcn, 0.5), sin(omega * 0.5));
}

TEST(ctfunc, cos)
{
    double omega = 2.;
    int fcn = func_new(CosFuncType, 0, 1, &omega);
    ASSERT_GE(fcn, 0);
    ASSERT_EQ(func_value(fcn, 0.), 1.);
    ASSERT_EQ(func_value(fcn, 0.5), cos(omega * 0.5));
}

TEST(ctfunc, sum)
{
    double omega = 2.;
    int fcn0 = func_new(SinFuncType, 0, 1, &omega);
    int fcn1 = func_new(CosFuncType, 0, 1, &omega);
    int fcn = func_new(SumFuncType, fcn0, fcn1, NULL);
    ASSERT_GE(fcn, 0);
    ASSERT_EQ(func_value(fcn, 0.), 1.);
    ASSERT_EQ(func_value(fcn, 0.5), sin(omega * 0.5) + cos(omega * 0.5));
}

TEST(ctfunc, diff)
{
    double omega = 2.;
    int fcn0 = func_new(SinFuncType, 0, 1, &omega);
    int fcn1 = func_new(CosFuncType, 0, 1, &omega);
    int fcn = func_new(DiffFuncType, fcn0, fcn1, NULL);
    ASSERT_GE(fcn, 0);
    ASSERT_EQ(func_value(fcn, 0.), -1.);
    ASSERT_EQ(func_value(fcn, 0.5), sin(omega * 0.5) - cos(omega * 0.5));
}

TEST(ctfunc, prod)
{
    double omega = 2.;
    int fcn0 = func_new(SinFuncType, 0, 1, &omega);
    int fcn1 = func_new(CosFuncType, 0, 1, &omega);
    int fcn = func_new(ProdFuncType, fcn0, fcn1, NULL);
    ASSERT_GE(fcn, 0);
    ASSERT_EQ(func_value(fcn, 0.), 0);
    ASSERT_EQ(func_value(fcn, 0.5), sin(omega * 0.5) * cos(omega * 0.5));
}

TEST(ctfunc, ratio)
{
    double omega = 2.;
    int fcn0 = func_new(SinFuncType, 0, 1, &omega);
    int fcn1 = func_new(CosFuncType, 0, 1, &omega);
    int fcn = func_new(RatioFuncType, fcn0, fcn1, NULL);
    ASSERT_GE(fcn, 0);
    ASSERT_EQ(func_value(fcn, 0.), 0.);
    ASSERT_EQ(func_value(fcn, 0.5), sin(omega * 0.5) / cos(omega * 0.5));
}

TEST(ctfunc, composite)
{
    double omega = 2.;
    int fcn0 = func_new(SinFuncType, 0, 1, &omega);
    int fcn1 = func_new(CosFuncType, 0, 1, &omega);
    int fcn = func_new(CompositeFuncType, fcn0, fcn1, NULL);
    ASSERT_GE(fcn, 0);
    ASSERT_EQ(func_value(fcn, 0.), sin(omega));
    ASSERT_EQ(func_value(fcn, 0.5), sin(omega * cos(omega * 0.5)));
}

TEST(ctfunc, times_constant)
{
    double omega = 2.;
    int fcn0 = func_new(SinFuncType, 0, 1, &omega);
    double A = 1.234;
    int fcn = func_new(TimesConstantFuncType, fcn0, 1, &A);
    ASSERT_GE(fcn, 0);
    ASSERT_EQ(func_value(fcn, 0.), 0.);
    ASSERT_EQ(func_value(fcn, 0.5), sin(omega * 0.5) * A);
}

TEST(ctfunc, plus_constant)
{
    double omega = 2.;
    int fcn0 = func_new(SinFuncType, 0, 1, &omega);
    double A = 1.234;
    int fcn = func_new(PlusConstantFuncType, fcn0, 1, &A);
    ASSERT_GE(fcn, 0);
    ASSERT_EQ(func_value(fcn, 0.), A);
    ASSERT_EQ(func_value(fcn, 0.5), sin(omega * 0.5) + A);
}
