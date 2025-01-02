#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <fstream>

#include "cantera/core.h"
#include "cantera/clib_experimental/ctfunc3.h"
#include "cantera/numerics/Func1Factory.h"

using namespace Cantera;


TEST(ctfunc3, invalid)
{
    // exceptions return -1
    ASSERT_EQ(func13_newBasic("spam", 0.), -2);
    ASSERT_EQ(func13_newAdvanced("eggs", 0, NULL), -2);
}

TEST(ctfunc3, sin)
{
    double omega = 2.1;
    int fcn = func13_newBasic("sin", omega);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func13_eval(fcn, 0.5), sin(omega * 0.5));

    int dfcn = func13_newDerivative(fcn);
    EXPECT_DOUBLE_EQ(func13_eval(dfcn, 0.5), omega * cos(omega * 0.5));

    int buflen = func13_write(fcn, "x", 0, 0);
    vector<char> buf(buflen);
    func13_write(fcn, "x", buflen, buf.data());
    string rep(buf.data());
    ASSERT_EQ(rep, "\\sin(2.1x)");
}

TEST(ctfunc3, cos)
{
    double omega = 2.;
    int fcn = func13_newBasic("cos", omega);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func13_eval(fcn, 0.5), cos(omega * 0.5));

    int dfcn = func13_newDerivative(fcn);
    EXPECT_DOUBLE_EQ(func13_eval(dfcn, 0.5), -omega * sin(omega * 0.5));
}

TEST(ctfunc3, exp)
{
    double omega = 2.;
    int fcn = func13_newBasic("exp", omega);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func13_eval(fcn, 0.5), exp(omega * 0.5));

    int dfcn = func13_newDerivative(fcn);
    EXPECT_DOUBLE_EQ(func13_eval(dfcn, 0.5), omega * exp(omega * 0.5));
}

TEST(ctfunc3, log)
{
    double omega = 2.;
    int fcn = func13_newBasic("log", omega);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func13_eval(fcn, 1. / omega), 0.);

    int dfcn = func13_newDerivative(fcn);
    EXPECT_DOUBLE_EQ(func13_eval(dfcn, .5), omega / .5);
}

TEST(ctfunc3, pow)
{
    double exp = .5;
    int fcn = func13_newBasic("pow", exp);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func13_eval(fcn, 0.5), pow(0.5, exp));

    int dfcn = func13_newDerivative(fcn);
    EXPECT_DOUBLE_EQ(func13_eval(dfcn, 0.5), exp * pow(0.5, exp - 1));
}

TEST(ctfunc3, constant)
{
    double a = .5;
    int fcn = func13_newBasic("constant", a);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func13_eval(fcn, 0.5), a);

    int dfcn = func13_newDerivative(fcn);
    EXPECT_DOUBLE_EQ(func13_eval(dfcn, .5), 0.);
}

TEST(ctfunc3, tabulated_linear)
{
    vector<double> params = {0., 1., 2., 1., 0., 1.};

    int fcn = func13_newAdvanced("tabulated-linear", params.size(), params.data());
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func13_eval(fcn, 0.5), 0.5);

    // exceptions return -1
    EXPECT_EQ(func13_newAdvanced("tabulated-linear", 5, params.data()), -2);
}

TEST(ctfunc3, tabulated_previous)
{
    vector<double> params = {0., 1., 2., 1., 0., 1.};

    int fcn = func13_newAdvanced("tabulated-previous", params.size(), params.data());
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func13_eval(fcn, 0.5), 1.);
}

TEST(ctfunc3, poly)
{
    double a0 = .5;
    double a1 = .25;
    double a2 = .125;
    vector<double> params = {a2, a1, a0};
    int fcn = func13_newAdvanced("polynomial3", params.size(), params.data());
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func13_eval(fcn, .5), (a2 * .5 + a1) * .5 + a0);

    params = {1, 0, -2.2, 3.1};
    fcn = func13_newAdvanced("polynomial3", params.size(), params.data());
    int buflen = func13_write(fcn, "x", 0, 0);
    vector<char> buf(buflen);
    func13_write(fcn, "x", buflen, buf.data());
    string rep(buf.data());
    ASSERT_EQ(rep, "x^3 - 2.2x + 3.1");
}

TEST(ctfunc3, Fourier)
{
    double a0 = .5;
    double a1 = .25;
    double b1 = .125;
    double omega = 2.;
    vector<double> params = {a0, a1, omega, b1};
    int fcn = func13_newAdvanced("Fourier", params.size(), params.data());
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(
        func13_eval(fcn, .5), .5 * a0 + a1 * cos(omega * .5) + b1 * sin(omega * .5));
}

TEST(ctfunc3, Gaussian)
{
    double A = .5;
    double t0 = .6;
    double fwhm = .25;
    vector<double> params = {A, t0, fwhm};
    int fcn = func13_newAdvanced("Gaussian", params.size(), params.data());
    ASSERT_GE(fcn, 0);
    double tau = fwhm / (2. * sqrt(log(2.)));
    double x = - t0 / tau;
    EXPECT_DOUBLE_EQ(func13_eval(fcn, 0.), A * exp(-x * x));

    // exceptions return -1
    EXPECT_EQ(func13_newAdvanced("Gaussian", 2, params.data()), -2);
}

TEST(ctfunc3, Arrhenius)
{
    double A = 38.7;
    double b = 2.7;
    double E = 2.619184e+07 / GasConstant;
    vector<double> params = {A, b, E};
    int fcn = func13_newAdvanced("Arrhenius", params.size(), params.data());
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func13_eval(fcn, 1000.), A * pow(1000., b) * exp(-E/1000.));
}

TEST(ctmath, invalid)
{
    // exceptions return -1
    int fcn0 = func13_newBasic("sin", 1.);
    int fcn1 = func13_newBasic("cos", 1.);
    ASSERT_EQ(func13_newCompound("foo", fcn0, fcn1), -2);
    ASSERT_EQ(func13_newModified("bar", fcn0, 0.), -2);
}

TEST(ctmath, sum)
{
    double omega = 2.;
    int fcn0 = func13_newBasic("sin", omega);
    int fcn1 = func13_newBasic("cos", omega);
    int fcn = func13_newCompound("sum", fcn0, fcn1);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func13_eval(fcn, 0.5), sin(omega * 0.5) + cos(omega * 0.5));
    int fcn2 = func13_newSum(fcn0, fcn1);
    EXPECT_DOUBLE_EQ(func13_eval(fcn, 0.5), func13_eval(fcn2, 0.5));
}

TEST(ctmath, diff)
{
    double omega = 2.;
    int fcn0 = func13_newBasic("sin", omega);
    int fcn1 = func13_newBasic("cos", omega);
    int fcn = func13_newCompound("diff", fcn0, fcn1);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func13_eval(fcn, 0.5), sin(omega * 0.5) - cos(omega * 0.5));
    int fcn2 = func13_newDiff(fcn0, fcn1);
    EXPECT_DOUBLE_EQ(func13_eval(fcn, 0.5), func13_eval(fcn2, 0.5));
}

TEST(ctmath, prod)
{
    double omega = 2.;
    int fcn0 = func13_newBasic("sin", omega);
    int fcn1 = func13_newBasic("cos", omega);
    int fcn = func13_newCompound("product", fcn0, fcn1);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func13_eval(fcn, 0.5), sin(omega * 0.5) * cos(omega * 0.5));
    int fcn2 = func13_newProd(fcn0, fcn1);
    EXPECT_DOUBLE_EQ(func13_eval(fcn, 0.5), func13_eval(fcn2, 0.5));
}

TEST(ctmath, ratio)
{
    double omega = 2.;
    int fcn0 = func13_newBasic("sin", omega);
    int fcn1 = func13_newBasic("cos", omega);
    int fcn = func13_newCompound("ratio", fcn0, fcn1);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func13_eval(fcn, 0.5), sin(omega * 0.5) / cos(omega * 0.5));
    int fcn2 = func13_newRatio(fcn0, fcn1);
    EXPECT_DOUBLE_EQ(func13_eval(fcn, 0.5), func13_eval(fcn2, 0.5));
}

TEST(ctmath, composite)
{
    double omega = 2.;
    int fcn0 = func13_newBasic("sin", omega);
    int fcn1 = func13_newBasic("cos", omega);
    int fcn = func13_newCompound("composite", fcn0, fcn1);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func13_eval(fcn, 0.5), sin(omega * cos(omega * 0.5)));
}

TEST(ctmath, times_constant)
{
    double omega = 2.;
    int fcn0 = func13_newBasic("sin", omega);
    double A = 1.234;
    int fcn = func13_newModified("times-constant", fcn0, A);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func13_eval(fcn, 0.5), sin(omega * 0.5) * A);
}

TEST(ctmath, plus_constant)
{
    double omega = 2.;
    int fcn0 = func13_newBasic("sin", omega);
    double A = 1.234;
    int fcn = func13_newModified("plus-constant", fcn0, A);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func13_eval(fcn, 0.5), sin(omega * 0.5) + A);
}

TEST(ctmath, periodic)
{
    double omega = 2.;
    int fcn0 = func13_newBasic("sin", omega);
    double A = 1.234;
    int fcn = func13_newModified("periodic", fcn0, A);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func13_eval(fcn, 0.5), func13_eval(fcn, 0.5 + A));
}

TEST(ctmath, generic)
{
    // Composite function reported in issue #752

    vector<double> params = {0., 0., 1., 1.};
    int fs = func13_newAdvanced("Fourier", params.size(), params.data()); // sin(x)
    auto func13_s = newFunc1("sin", 1.);
    EXPECT_DOUBLE_EQ(func13_eval(fs, 0.5), func13_s->eval(0.5));

    int fs2 = func13_newCompound("product", fs, fs);  // (sin(x)^2)
    auto func13_s2 = newFunc1("product", func13_s, func13_s);
    EXPECT_DOUBLE_EQ(func13_eval(fs2, 0.5), func13_s2->eval(0.5));

    double one = 1.;
    int fs1 = func13_newAdvanced("polynomial3", 1, &one); // 1.
    auto func13_s1 = newFunc1("constant", 1.);
    EXPECT_DOUBLE_EQ(func13_eval(fs1, 0.5), func13_s1->eval(0.5));
    EXPECT_DOUBLE_EQ(func13_eval(fs1, 0.5), 1.);

    int fs2_1 = func13_newCompound("diff", fs1, fs2);  // 1-(sin(x))^2
    auto func13_s2_1 = newFunc1("diff", func13_s1, func13_s2);
    EXPECT_DOUBLE_EQ(func13_eval(fs2_1, 0.5), func13_s2_1->eval(0.5));

    vector<double> p_arr = {1., .5, 0.};
    int fs3 = func13_newAdvanced("Arrhenius", 3, p_arr.data());  // sqrt function
    auto func13_s3 = newFunc1("Arrhenius", p_arr);
    EXPECT_DOUBLE_EQ(func13_eval(fs3, 0.5), func13_s3->eval(0.5));

    // overall composite function
    int fs4 = func13_newCompound("composite", fs3, fs2_1);  // sqrt(1-(sin(x))^2)
    auto func13_s4 = newFunc1("composite", func13_s3, func13_s2_1);
    EXPECT_DOUBLE_EQ(func13_eval(fs4, 0.5), func13_s4->eval(0.5));

    // an easier equivalent expression (using trigonometry)
    auto func13_s5 = newFunc1("cos", 1.);  // missing the absolute value
    EXPECT_DOUBLE_EQ(func13_eval(fs4, 0.5), func13_s5->eval(0.5));
    EXPECT_DOUBLE_EQ(func13_eval(fs4, 3.5), -func13_s5->eval(3.5));
}
