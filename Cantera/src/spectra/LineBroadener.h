/**
 * @file LinerBoadener.h
 * Header file for class LineBroadener
 */

#include "ct_defs.h"
#include "ctexceptions.h"

namespace Cantera {

    /**
     * Base class for classes implementing line shapes of
     * various types.
     */
    class LineBroadener {

    public:        

        LineBroadener() {}

        virtual ~LineBroadener() {}

        /**
         * The line shape profile
         * \f[
         * P(\Delta\nu)
         *\f]
         * as a function of distance from line
         * center. This function must have total area = 1.0.
         * Note that this method must be overloaded in each
         * derived class. If the base class method is called,
         * an exception will be thrown.
         */
        virtual doublereal profile(doublereal deltaFreq) {
            throw CanteraError("LineBroadener::profile",
                "base class method called!");
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

    };

    /**
     * The line shape for pure collisional broadening. The Lorentzian line
     * shape is
     * \f[
     * L(\Delta\nu) = \frac{1}{\pi}\frac{\gamma}{\Delta\nu^2 + \gamma^2}
     * \f]
     * where \f$ \gamma = {\mbox{FWHM}/2 \f$.  
     */
    class Lorentzian : public LineBroadener {
    public:
        Lorentzian(doublereal FWHM);
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
    class Gaussian : public LineBroadener {
    public:

        /**
         * Constructor.
         * @param FWHM Full width at half-maximum.
         */
        Gaussian(doublereal FWHM);
        virtual doublereal profile(doublereal deltaFreq);
        virtual doublereal cumulative(doublereal deltaFreq);
        virtual doublereal width();

    protected:
        doublereal m_sigma;
        doublereal m_sigma2;
        doublereal m_width;
    };

}
