/**
 * @file LineBroadener.h
 * Header file for class LineBroadener
 * @ingroup spectroscopy
 */

#include "cantera/base/ct_defs.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/global.h"

namespace Cantera
{

/**
 * Base class for classes implementing line shapes of
 * various types.
 * @ingroup spectroscopy
 * @deprecated incomplete / abandoned
 */
class LineBroadener
{

public:

    /// Default constructor
    LineBroadener() {
        warn_deprecated("class LineBroadener");
    }

    /// Destructor
    virtual ~LineBroadener() {}

    /**
     * The line shape profile
     * \f[
     * P(\Delta\nu)
     *\f]
     * as a function of distance from line
     * center \f$ \Delta\nu \f$.
     * This function must have total area = 1.0.
     * Note that this method must be overloaded in each
     * derived class. If the base class method is called,
     * an exception will be thrown.
     */
    virtual doublereal profile(doublereal deltaFreq) {
        throw CanteraError("LineBroadener::profile",
                           "base class method called!");
    }

    doublereal operator()(doublereal deltaFreq) {
        return profile(deltaFreq);
    }

    /**
     * The cumulative profile, defined as
     * \f[
     * C(\Delta \nu) = \int_{-\infty}^{\Delta \nu} P(x) dx
     * \f]
     */
    virtual doublereal cumulative(doublereal deltaFreq) {
        throw CanteraError("LineBroadener::cumulative",
                           "base class method called!");
    }

    virtual doublereal width() {
        return 0.0;
    }
};

/**
 * The line shape for pure collisional broadening. The Lorentzian line
 * shape is
 * \f[
 * L(\Delta\nu) = \frac{1}{\pi}\frac{\gamma}{\Delta\nu^2 + \gamma^2}
 * \f]
 * where \f$ \gamma = {\mbox{FWHM}/2} \f$.
 */
class LorentzianProfile : public LineBroadener
{
public:
    LorentzianProfile(doublereal FWHM);
    virtual doublereal profile(doublereal deltaFreq);
    virtual doublereal cumulative(doublereal deltaFreq);
    virtual double width();

protected:
    doublereal m_hwhm;
    doublereal m_hwhm2;
};

/**
 * A Gaussian line profile. This profile results when Doppler
 * broadening is dominant.
 */
class GaussianProfile : public LineBroadener
{
public:

    /**
     * Constructor.
     */
    GaussianProfile(doublereal sigma);
    virtual doublereal profile(doublereal deltaFreq);
    virtual doublereal cumulative(doublereal deltaFreq);
    virtual doublereal width();

    doublereal standardDev() {
        return m_sigma;
    }

protected:
    doublereal m_sigma;
    doublereal m_sigma2;
    doublereal m_width;
};


/**
 * A Voigt profile is the convolution of a Lorentzian and a
 * Gaussian profile. This profile results when Doppler
 * broadening and collisional broadening both are important.
 */
class Voigt : public LineBroadener
{
public:

    /**
     * Constructor.
     */
    Voigt(doublereal sigma, doublereal gamma);
    virtual doublereal profile(doublereal deltaFreq);
    //virtual doublereal cumulative(doublereal deltaFreq)
    //virtual doublereal width()

    void testv();

protected:

    doublereal F(doublereal x);

    doublereal m_sigma;
    doublereal m_gamma_lor;
    doublereal m_sigma2;
    doublereal m_width;
    doublereal m_gamma;
    doublereal m_sigsqrt2;
    doublereal m_a;
    doublereal m_eps;
};

}
