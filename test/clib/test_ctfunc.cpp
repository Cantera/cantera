#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <fstream>

#include "cantera/core.h"
#include "cantera/clib/ctfunc.h"
#include "cantera/numerics/Func1Factory.h"

using namespace Cantera;


TEST(ctfunc, invalid)
{
    // exceptions return -1
    ASSERT_EQ(func_new_basic("spam", 0.), -1);
    ASSERT_EQ(func_new_advanced("eggs", 0, NULL), -1);
}

TEST(ctfunc, sin)
{
    double omega = 2.1;
    int fcn = func_new_basic("sin", omega);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func_value(fcn, 0.5), sin(omega * 0.5));

    int dfcn = func_derivative(fcn);
    EXPECT_DOUBLE_EQ(func_value(dfcn, 0.5), omega * cos(omega * 0.5));

    int buflen = func_write(fcn, "x", 0, 0);
    vector<char> buf(buflen);
    func_write(fcn, "x", buflen, buf.data());
    string rep(buf.data());
    ASSERT_EQ(rep, "\\sin(2.1x)");
}

TEST(ctfunc, cos)
{
    double omega = 2.;
    int fcn = func_new_basic("cos", omega);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func_value(fcn, 0.5), cos(omega * 0.5));

    int dfcn = func_derivative(fcn);
    EXPECT_DOUBLE_EQ(func_value(dfcn, 0.5), -omega * sin(omega * 0.5));
}

TEST(ctfunc, exp)
{
    double omega = 2.;
    int fcn = func_new_basic("exp", omega);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func_value(fcn, 0.5), exp(omega * 0.5));

    int dfcn = func_derivative(fcn);
    EXPECT_DOUBLE_EQ(func_value(dfcn, 0.5), omega * exp(omega * 0.5));
}

TEST(ctfunc, log)
{
    double omega = 2.;
    int fcn = func_new_basic("log", omega);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func_value(fcn, 1. / omega), 0.);

    int dfcn = func_derivative(fcn);
    EXPECT_DOUBLE_EQ(func_value(dfcn, .5), omega / .5);
}

TEST(ctfunc, pow)
{
    double exp = .5;
    int fcn = func_new_basic("pow", exp);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func_value(fcn, 0.5), pow(0.5, exp));

    int dfcn = func_derivative(fcn);
    EXPECT_DOUBLE_EQ(func_value(dfcn, 0.5), exp * pow(0.5, exp - 1));
}

TEST(ctfunc, constant)
{
    double a = .5;
    int fcn = func_new_basic("constant", a);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func_value(fcn, 0.5), a);

    int dfcn = func_derivative(fcn);
    EXPECT_DOUBLE_EQ(func_value(dfcn, .5), 0.);
}

TEST(ctfunc, tabulated_linear)
{
    vector<double> params = {0., 1., 2., 1., 0., 1.};

    int fcn = func_new_advanced("tabulated-linear", params.size(), params.data());
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func_value(fcn, 0.5), 0.5);

    // exceptions return -1
    EXPECT_EQ(func_new_advanced("tabulated-linear", 5, params.data()), -1);
}

TEST(ctfunc, tabulated_previous)
{
    vector<double> params = {0., 1., 2., 1., 0., 1.};

    int fcn = func_new_advanced("tabulated-previous", params.size(), params.data());
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func_value(fcn, 0.5), 1.);
}

TEST(ctfunc, poly)
{
    double a0 = .5;
    double a1 = .25;
    double a2 = .125;
    vector<double> params = {a2, a1, a0};
    int fcn = func_new_advanced("polynomial3", params.size(), params.data());
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func_value(fcn, .5), (a2 * .5 + a1) * .5 + a0);

    params = {1, 0, -2.2, 3.1};
    fcn = func_new_advanced("polynomial3", params.size(), params.data());
    int buflen = func_write(fcn, "x", 0, 0);
    vector<char> buf(buflen);
    func_write(fcn, "x", buflen, buf.data());
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
    int fcn = func_new_advanced("Fourier", params.size(), params.data());
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(
        func_value(fcn, .5), .5 * a0 + a1 * cos(omega * .5) + b1 * sin(omega * .5));
}

TEST(ctfunc, Gaussian)
{
    double A = .5;
    double t0 = .6;
    double fwhm = .25;
    vector<double> params = {A, t0, fwhm};
    int fcn = func_new_advanced("Gaussian", params.size(), params.data());
    ASSERT_GE(fcn, 0);
    double tau = fwhm / (2. * sqrt(log(2.)));
    double x = - t0 / tau;
    EXPECT_DOUBLE_EQ(func_value(fcn, 0.), A * exp(-x * x));

    // exceptions return -1
    EXPECT_EQ(func_new_advanced("Gaussian", 2, params.data()), -1);
}

TEST(ctfunc, Arrhenius)
{
    double A = 38.7;
    double b = 2.7;
    double E = 2.619184e+07 / GasConstant;
    vector<double> params = {A, b, E};
    int fcn = func_new_advanced("Arrhenius", params.size(), params.data());
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func_value(fcn, 1000.), A * pow(1000., b) * exp(-E/1000.));
}

TEST(ctmath, invalid)
{
    // exceptions return -1
    int fcn0 = func_new_basic("sin", 1.);
    int fcn1 = func_new_basic("cos", 1.);
    ASSERT_EQ(func_new_compound("foo", fcn0, fcn1), -1);
    ASSERT_EQ(func_new_modified("bar", fcn0, 0.), -1);
}

TEST(ctmath, sum)
{
    double omega = 2.;
    int fcn0 = func_new_basic("sin", omega);
    int fcn1 = func_new_basic("cos", omega);
    int fcn = func_new_compound("sum", fcn0, fcn1);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func_value(fcn, 0.5), sin(omega * 0.5) + cos(omega * 0.5));
    int fcn2 = func_new_sum(fcn0, fcn1);
    EXPECT_DOUBLE_EQ(func_value(fcn, 0.5), func_value(fcn2, 0.5));
}

TEST(ctmath, diff)
{
    double omega = 2.;
    int fcn0 = func_new_basic("sin", omega);
    int fcn1 = func_new_basic("cos", omega);
    int fcn = func_new_compound("diff", fcn0, fcn1);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func_value(fcn, 0.5), sin(omega * 0.5) - cos(omega * 0.5));
    int fcn2 = func_new_diff(fcn0, fcn1);
    EXPECT_DOUBLE_EQ(func_value(fcn, 0.5), func_value(fcn2, 0.5));
}

TEST(ctmath, prod)
{
    double omega = 2.;
    int fcn0 = func_new_basic("sin", omega);
    int fcn1 = func_new_basic("cos", omega);
    int fcn = func_new_compound("product", fcn0, fcn1);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func_value(fcn, 0.5), sin(omega * 0.5) * cos(omega * 0.5));
    int fcn2 = func_new_prod(fcn0, fcn1);
    EXPECT_DOUBLE_EQ(func_value(fcn, 0.5), func_value(fcn2, 0.5));
}

TEST(ctmath, ratio)
{
    double omega = 2.;
    int fcn0 = func_new_basic("sin", omega);
    int fcn1 = func_new_basic("cos", omega);
    int fcn = func_new_compound("ratio", fcn0, fcn1);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func_value(fcn, 0.5), sin(omega * 0.5) / cos(omega * 0.5));
    int fcn2 = func_new_ratio(fcn0, fcn1);
    EXPECT_DOUBLE_EQ(func_value(fcn, 0.5), func_value(fcn2, 0.5));
}

TEST(ctmath, composite)
{
    double omega = 2.;
    int fcn0 = func_new_basic("sin", omega);
    int fcn1 = func_new_basic("cos", omega);
    int fcn = func_new_compound("composite", fcn0, fcn1);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func_value(fcn, 0.5), sin(omega * cos(omega * 0.5)));
}

TEST(ctmath, times_constant)
{
    double omega = 2.;
    int fcn0 = func_new_basic("sin", omega);
    double A = 1.234;
    int fcn = func_new_modified("times-constant", fcn0, A);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func_value(fcn, 0.5), sin(omega * 0.5) * A);
}

TEST(ctmath, plus_constant)
{
    double omega = 2.;
    int fcn0 = func_new_basic("sin", omega);
    double A = 1.234;
    int fcn = func_new_modified("plus-constant", fcn0, A);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func_value(fcn, 0.5), sin(omega * 0.5) + A);
}

TEST(ctmath, periodic)
{
    double omega = 2.;
    int fcn0 = func_new_basic("sin", omega);
    double A = 1.234;
    int fcn = func_new_modified("periodic", fcn0, A);
    ASSERT_GE(fcn, 0);
    EXPECT_DOUBLE_EQ(func_value(fcn, 0.5), func_value(fcn, 0.5 + A));
}

TEST(ctmath, generic)
{
    // Composite function reported in issue #752

    vector<double> params = {0., 0., 1., 1.};
    int fs = func_new_advanced("Fourier", params.size(), params.data()); // sin(x)
    auto func_s = newFunc1("sin", 1.);
    EXPECT_DOUBLE_EQ(func_value(fs, 0.5), func_s->eval(0.5));

    int fs2 = func_new_compound("product", fs, fs);  // (sin(x)^2)
    auto func_s2 = newFunc1("product", func_s, func_s);
    EXPECT_DOUBLE_EQ(func_value(fs2, 0.5), func_s2->eval(0.5));

    double one = 1.;
    int fs1 = func_new_advanced("polynomial3", 1, &one); // 1.
    auto func_s1 = newFunc1("constant", 1.);
    EXPECT_DOUBLE_EQ(func_value(fs1, 0.5), func_s1->eval(0.5));
    EXPECT_DOUBLE_EQ(func_value(fs1, 0.5), 1.);

    int fs2_1 = func_new_compound("diff", fs1, fs2);  // 1-(sin(x))^2
    auto func_s2_1 = newFunc1("diff", func_s1, func_s2);
    EXPECT_DOUBLE_EQ(func_value(fs2_1, 0.5), func_s2_1->eval(0.5));

    vector<double> p_arr = {1., .5, 0.};
    int fs3 = func_new_advanced("Arrhenius", 3, p_arr.data());  // sqrt function
    auto func_s3 = newFunc1("Arrhenius", p_arr);
    EXPECT_DOUBLE_EQ(func_value(fs3, 0.5), func_s3->eval(0.5));

    // overall composite function
    int fs4 = func_new_compound("composite", fs3, fs2_1);  // sqrt(1-(sin(x))^2)
    auto func_s4 = newFunc1("composite", func_s3, func_s2_1);
    EXPECT_DOUBLE_EQ(func_value(fs4, 0.5), func_s4->eval(0.5));

    // an easier equivalent expression (using trigonometry)
    auto func_s5 = newFunc1("cos", 1.);  // missing the absolute value
    EXPECT_DOUBLE_EQ(func_value(fs4, 0.5), func_s5->eval(0.5));
    EXPECT_DOUBLE_EQ(func_value(fs4, 3.5), -func_s5->eval(3.5));
}
