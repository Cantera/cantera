#include "ct_defs.h"
#include "ctexceptions.h"

namespace Cantera {

    class LineBroadener {

    public:        

        LineBroadener() {}

        virtual ~LineBroadener() {}

        /**
         * the line shape profile, as a function of distance from line
         * center. This function must have total area = 1.0.
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

    class Gaussian : public LineBroadener {
    public:
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
