#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <fstream>

#include "cantera/core.h"
#include "cantera_clib/ctfunc.h"
#include "cantera/numerics/Func1Factory.h"

using namespace Cantera;


TEST(ctfunc, invalid)
{
    // exceptions return -1
    ASSERT_EQ(func1_newBasic("spam", 0.), -2);
    ASSERT_EQ(func1_newAdvanced("eggs", 0, NULL), -2);
}

TEST(ctfunc, sin)
{
    double omega = 2.1;
    int32_t fcn = func1_newBasic("sin", omega);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func1_eval(fcn, 0.5), sin(omega * 0.5));

    int32_t dfcn = func1_derivative(fcn);
    EXPECT_DOUBLE_EQ(func1_eval(dfcn, 0.5), omega * cos(omega * 0.5));

    int32_t buflen = func1_write(fcn, "x", 0, 0);
    vector<char> buf(buflen);
    func1_write(fcn, "x", buflen, buf.data());
    string rep(buf.data());
    ASSERT_EQ(rep, "\\sin(2.1x)");
}

TEST(ctfunc, cos)
{
    double omega = 2.;
    int32_t fcn = func1_newBasic("cos", omega);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func1_eval(fcn, 0.5), cos(omega * 0.5));

    int32_t dfcn = func1_derivative(fcn);
    EXPECT_DOUBLE_EQ(func1_eval(dfcn, 0.5), -omega * sin(omega * 0.5));
}

TEST(ctfunc, exp)
{
    double omega = 2.;
    int32_t fcn = func1_newBasic("exp", omega);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func1_eval(fcn, 0.5), exp(omega * 0.5));

    int32_t dfcn = func1_derivative(fcn);
    EXPECT_DOUBLE_EQ(func1_eval(dfcn, 0.5), omega * exp(omega * 0.5));
}

TEST(ctfunc, log)
{
    double omega = 2.;
    int32_t fcn = func1_newBasic("log", omega);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func1_eval(fcn, 1. / omega), 0.);

    int32_t dfcn = func1_derivative(fcn);
    EXPECT_DOUBLE_EQ(func1_eval(dfcn, .5), omega / .5);
}

TEST(ctfunc, pow)
{
    double exp = .5;
    int32_t fcn = func1_newBasic("pow", exp);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func1_eval(fcn, 0.5), pow(0.5, exp));

    int32_t dfcn = func1_derivative(fcn);
    EXPECT_DOUBLE_EQ(func1_eval(dfcn, 0.5), exp * pow(0.5, exp - 1));
}

TEST(ctfunc, constant)
{
    double a = .5;
    int32_t fcn = func1_newBasic("constant", a);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func1_eval(fcn, 0.5), a);

    int32_t dfcn = func1_derivative(fcn);
    EXPECT_DOUBLE_EQ(func1_eval(dfcn, .5), 0.);
}

TEST(ctfunc, tabulated_linear)
{
    vector<double> params = {0., 1., 2., 1., 0., 1.};

    int32_t fcn = func1_newAdvanced("tabulated-linear", params.size(), params.data());
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func1_eval(fcn, 0.5), 0.5);

    // exceptions return -1
    EXPECT_EQ(func1_newAdvanced("tabulated-linear", 5, params.data()), -2);
}

TEST(ctfunc, tabulated_previous)
{
    vector<double> params = {0., 1., 2., 1., 0., 1.};

    int32_t fcn = func1_newAdvanced("tabulated-previous", params.size(), params.data());
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func1_eval(fcn, 0.5), 1.);
}

TEST(ctfunc, poly)
{
    double a0 = .5;
    double a1 = .25;
    double a2 = .125;
    vector<double> params = {a2, a1, a0};
    int32_t fcn = func1_newAdvanced("polynomial3", params.size(), params.data());
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func1_eval(fcn, .5), (a2 * .5 + a1) * .5 + a0);

    params = {1, 0, -2.2, 3.1};
    fcn = func1_newAdvanced("polynomial3", params.size(), params.data());
    int32_t buflen = func1_write(fcn, "x", 0, 0);
    vector<char> buf(buflen);
    func1_write(fcn, "x", buflen, buf.data());
    string rep(buf.data());
    ASSERT_EQ(rep, "x^3 - 2.2x + 3.1");
}

TEST(ctfunc, Fourier)
{
    double a0 = .5;
    double a1 = .25;
    double b1 = .125;
    double omega = 2.;
    vector<double> params = {a0, a1, omega, b1};
    int32_t fcn = func1_newAdvanced("Fourier", params.size(), params.data());
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(
        func1_eval(fcn, .5), .5 * a0 + a1 * cos(omega * .5) + b1 * sin(omega * .5));
}

TEST(ctfunc, Gaussian)
{
    double A = .5;
    double t0 = .6;
    double fwhm = .25;
    vector<double> params = {A, t0, fwhm};
    int32_t fcn = func1_newAdvanced("Gaussian", params.size(), params.data());
    ASSERT_GE(fcn, 0);
    double tau = fwhm / (2. * sqrt(log(2.)));
    double x = - t0 / tau;
    EXPECT_DOUBLE_EQ(func1_eval(fcn, 0.), A * exp(-x * x));

    // exceptions return -1
    EXPECT_EQ(func1_newAdvanced("Gaussian", 2, params.data()), -2);
}

TEST(ctfunc, Arrhenius)
{
    double A = 38.7;
    double b = 2.7;
    double E = 2.619184e+07 / GasConstant;
    vector<double> params = {A, b, E};
    int32_t fcn = func1_newAdvanced("Arrhenius", params.size(), params.data());
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func1_eval(fcn, 1000.), A * pow(1000., b) * exp(-E/1000.));
}

TEST(ctmath, invalid)
{
    // exceptions return -1
    int32_t fcn0 = func1_newBasic("sin", 1.);
    int32_t fcn1 = func1_newBasic("cos", 1.);
    ASSERT_EQ(func1_newCompound("foo", fcn0, fcn1), -2);
    ASSERT_EQ(func1_newModified("bar", fcn0, 0.), -2);
}

TEST(ctmath, sum)
{
    double omega = 2.;
    int32_t fcn0 = func1_newBasic("sin", omega);
    int32_t fcn1 = func1_newBasic("cos", omega);
    int32_t fcn = func1_newCompound("sum", fcn0, fcn1);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func1_eval(fcn, 0.5), sin(omega * 0.5) + cos(omega * 0.5));
    int32_t fcn2 = func1_newSumFunction(fcn0, fcn1);
    EXPECT_DOUBLE_EQ(func1_eval(fcn, 0.5), func1_eval(fcn2, 0.5));
}

TEST(ctmath, diff)
{
    double omega = 2.;
    int32_t fcn0 = func1_newBasic("sin", omega);
    int32_t fcn1 = func1_newBasic("cos", omega);
    int32_t fcn = func1_newCompound("diff", fcn0, fcn1);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func1_eval(fcn, 0.5), sin(omega * 0.5) - cos(omega * 0.5));
    int32_t fcn2 = func1_newDiffFunction(fcn0, fcn1);
    EXPECT_DOUBLE_EQ(func1_eval(fcn, 0.5), func1_eval(fcn2, 0.5));
}

TEST(ctmath, prod)
{
    double omega = 2.;
    int32_t fcn0 = func1_newBasic("sin", omega);
    int32_t fcn1 = func1_newBasic("cos", omega);
    int32_t fcn = func1_newCompound("product", fcn0, fcn1);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func1_eval(fcn, 0.5), sin(omega * 0.5) * cos(omega * 0.5));
    int32_t fcn2 = func1_newProdFunction(fcn0, fcn1);
    EXPECT_DOUBLE_EQ(func1_eval(fcn, 0.5), func1_eval(fcn2, 0.5));
}

TEST(ctmath, ratio)
{
    double omega = 2.;
    int32_t fcn0 = func1_newBasic("sin", omega);
    int32_t fcn1 = func1_newBasic("cos", omega);
    int32_t fcn = func1_newCompound("ratio", fcn0, fcn1);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func1_eval(fcn, 0.5), sin(omega * 0.5) / cos(omega * 0.5));
    int32_t fcn2 = func1_newRatioFunction(fcn0, fcn1);
    EXPECT_DOUBLE_EQ(func1_eval(fcn, 0.5), func1_eval(fcn2, 0.5));
}

TEST(ctmath, composite)
{
    double omega = 2.;
    int32_t fcn0 = func1_newBasic("sin", omega);
    int32_t fcn1 = func1_newBasic("cos", omega);
    int32_t fcn = func1_newCompound("composite", fcn0, fcn1);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func1_eval(fcn, 0.5), sin(omega * cos(omega * 0.5)));
}

TEST(ctmath, times_constant)
{
    double omega = 2.;
    int32_t fcn0 = func1_newBasic("sin", omega);
    double A = 1.234;
    int32_t fcn = func1_newModified("times-constant", fcn0, A);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func1_eval(fcn, 0.5), sin(omega * 0.5) * A);
}

TEST(ctmath, plus_constant)
{
    double omega = 2.;
    int32_t fcn0 = func1_newBasic("sin", omega);
    double A = 1.234;
    int32_t fcn = func1_newModified("plus-constant", fcn0, A);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func1_eval(fcn, 0.5), sin(omega * 0.5) + A);
}

TEST(ctmath, periodic)
{
    double omega = 2.;
    int32_t fcn0 = func1_newBasic("sin", omega);
    double A = 1.234;
    int32_t fcn = func1_newModified("periodic", fcn0, A);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func1_eval(fcn, 0.5), func1_eval(fcn, 0.5 + A));
}

TEST(ctmath, generic)
{
    // Composite function reported in issue #752

    vector<double> params = {0., 0., 1., 1.};
    int32_t fs = func1_newAdvanced("Fourier", params.size(), params.data()); // sin(x)
    auto func1_s = newFunc1("sin", 1.);
    EXPECT_DOUBLE_EQ(func1_eval(fs, 0.5), func1_s->eval(0.5));

    int32_t fs2 = func1_newCompound("product", fs, fs);  // (sin(x)^2)
    auto func1_s2 = newFunc1("product", func1_s, func1_s);
    EXPECT_DOUBLE_EQ(func1_eval(fs2, 0.5), func1_s2->eval(0.5));

    double one = 1.;
    int32_t fs1 = func1_newAdvanced("polynomial3", 1, &one); // 1.
    auto func1_s1 = newFunc1("constant", 1.);
    EXPECT_DOUBLE_EQ(func1_eval(fs1, 0.5), func1_s1->eval(0.5));
    EXPECT_DOUBLE_EQ(func1_eval(fs1, 0.5), 1.);

    int32_t fs2_1 = func1_newCompound("diff", fs1, fs2);  // 1-(sin(x))^2
    auto func1_s2_1 = newFunc1("diff", func1_s1, func1_s2);
    EXPECT_DOUBLE_EQ(func1_eval(fs2_1, 0.5), func1_s2_1->eval(0.5));

    vector<double> p_arr = {1., .5, 0.};
    int32_t fs3 = func1_newAdvanced("Arrhenius", 3, p_arr.data());  // sqrt function
    auto func1_s3 = newFunc1("Arrhenius", p_arr);
    EXPECT_DOUBLE_EQ(func1_eval(fs3, 0.5), func1_s3->eval(0.5));

    // overall composite function
    int32_t fs4 = func1_newCompound("composite", fs3, fs2_1);  // sqrt(1-(sin(x))^2)
    auto func1_s4 = newFunc1("composite", func1_s3, func1_s2_1);
    EXPECT_DOUBLE_EQ(func1_eval(fs4, 0.5), func1_s4->eval(0.5));

    // an easier equivalent expression (using trigonometry)
    auto func1_s5 = newFunc1("cos", 1.);  // missing the absolute value
    EXPECT_DOUBLE_EQ(func1_eval(fs4, 0.5), func1_s5->eval(0.5));
    EXPECT_DOUBLE_EQ(func1_eval(fs4, 3.5), -func1_s5->eval(3.5));
}
