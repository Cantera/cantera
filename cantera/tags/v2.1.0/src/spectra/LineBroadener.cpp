#include "cantera/base/ct_defs.h"
#include <math.h>
#include <iostream>

#ifdef USE_BOOST_MATH
#include <boost/math/special_functions/erf.hpp>
using boost::math::erf;
#endif

#include "cantera/spectra/LineBroadener.h"

namespace Cantera
{

LorentzianProfile::LorentzianProfile(doublereal gamma)
{
    m_hwhm = gamma;
    m_hwhm2 = m_hwhm*m_hwhm;
}

/**
 * The Lorentzian profile for collision-broadened lines.
 *
 *\f[
 * \frac{1}{\pi} \frac{\gamma}{ (\Delta\nu)^2 + \gamma^2}
 *\f]
 * Units: 1/wavenumber (or cm).
 */
doublereal LorentzianProfile::profile(doublereal deltaFreq)
{
    return (1.0/Cantera::Pi) *m_hwhm/(deltaFreq*deltaFreq + m_hwhm2);
}

/**
 *
 * The cumulative profile, given by
 * \f[
 * \frac{1}{\pi} \tan^{-1}\left(\frac{\Delta\nu}{gamma}\right) + 0.5
 * \f]
 */
doublereal LorentzianProfile::cumulative(doublereal deltaFreq)
{
    return (1.0/Pi) * atan(deltaFreq/m_hwhm) + 0.5;
}

doublereal LorentzianProfile::width()
{
    return 2.0*m_hwhm;
}

GaussianProfile::GaussianProfile(doublereal sigma)
{
    m_sigma = sigma;
    m_sigma2 = m_sigma*m_sigma;
}

doublereal GaussianProfile::profile(doublereal deltaFreq)
{
    //cout << "entered Gaussian::profile" << endl;
    //cout << "deltaFreq = " << deltaFreq << endl;
    //cout << "m_sigma = " << m_sigma << endl;
    return 1.0/(m_sigma * Cantera::SqrtTwo *Cantera::SqrtPi) *
           exp(-deltaFreq*deltaFreq/(2.0*m_sigma2));
}

doublereal GaussianProfile::cumulative(doublereal deltaFreq)
{
    return 0.5*(1.0 + erf(deltaFreq/(m_sigma*SqrtTwo)));
}

doublereal GaussianProfile::width()
{
    return 2.0*m_sigma*sqrt(log(4.0));
}



/**
 * @param sigma The standard deviation of the Gaussian
 * @param gamma The half-width of the Lorentzian.
 */
Voigt::Voigt(doublereal sigma, doublereal gamma)
{
    m_sigma = sigma;
    m_sigma2 = m_sigma*m_sigma;
    m_gamma_lor = gamma;
    m_sigsqrt2 = SqrtTwo*m_sigma;
    m_gamma = gamma/m_sigsqrt2;
    m_eps = 1.0e-20;
}

void Voigt::testv()
{
    m_gamma = 1.0e1;
    std::cout << F(1.0) << std::endl;
    m_gamma = 0.5;
    std::cout << F(1.0) << std::endl;
    m_gamma = 0.0001;
    std::cout << F(10.0) << std::endl;
}
/**
 * This method evaluates the function
 * \f[
 * F(x, y) = \frac{y}{\pi}\int_{-\infty}^{+\infty} \frac{e^{-z^2}}
 * {(x - z)^2 + y^2} dz
 * \f]
 * The algorithm used to cmpute this function is described in the
 * reference below.  @see F. G. Lether and P. R. Wenston, "The
 * numerical computation of the %Voigt function by a corrected
 * midpoint quadrature rule for \f$ (-\infty, \infty) \f$. Journal
 * of Computational and Applied Mathematics}, 34 (1):75--92, 1991.
 */
doublereal Voigt::F(doublereal x)
{

    if (x < 0.0) {
        x = -x;
    }
    double y = m_gamma;

    double c3 = log(Pi*m_eps/2.0);
    double tau = sqrt(-log(y) - c3);
    double b = (tau + x)/y;
    double t = b*y;
    double f1, f2, f3;
    const double c0 = 2.0/(Pi*exp(0.0));
    const double c1 = 1.0/SqrtTwo;
    const double c2 = 2.0/SqrtPi;

    if (y > c0/m_eps) {
        return 0.0;
    }
    double f0, ef0;
    while (1 > 0) {
        f0 = Pi*Pi/(t*t);
        ef0 = exp(-f0);
        f1 = c2*y*ef0;
        f2 = fabs(y*y - Pi*Pi/(t*t));
        f3 = 1.0 - ef0*ef0;
        t *= c1;
        if (f1/(f2*f3) < 0.5*m_eps) {
            break;
        }
    }
    double h = t/y;
    int N = int(0.5 + b/h);
    double S = 0.0;
    double u = h/2;
    for (int i = 0; i < N; i++) {
        S += (1.0 + exp(-4.0*x*y*u))*exp(-pow(y*u-x,2))/(u*u+1.0);
        u += h;
    }
    double Q = h*S/Pi;
    double C = 0.0;
    if (y*y < Pi/h) {
        C = 2.0*exp(y*y - x*x)*cos(2*x*y)/(1.0 + exp(2*Pi/h));
    } else {
        return 0.0;
    }
    return Q + C;
}

/**
 * Voigt profile.
 *
 * Not sure that constant is right.
 */
doublereal Voigt::profile(doublereal deltaFreq)
{
    const double ff = 1.0/(m_sigsqrt2*SqrtPi);
    return ff*F(deltaFreq/m_sigsqrt2);
}

}
