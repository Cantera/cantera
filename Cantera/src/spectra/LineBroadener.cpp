#include "ct_defs.h"
#include "LineBroadener.h"

using namespace std;

namespace Cantera {

    doublereal to_wavenumbers(doublereal freq) {
        return freq/(100.0*lightSpeed);
    }

    doublereal from_wavenumbers(doublereal omega) {
        return omega*(100.0*lightSpeed);
    }

    Lorentzian::Lorentzian(doublereal FWHM) {
        m_hwhm = 0.5*FWHM;
        m_hwhm2 = m_hwhm*m_hwhm;
    }

    doublereal Lorentzian::profile(doublereal deltaFreq) {
        return (1.0/Pi) *m_hwhm/(deltaFreq*deltaFreq + m_hwhm2);
    }
        
    doublereal Lorentzian::cumulative(doublereal deltaFreq) {
        return (1.0/Pi) * atan(deltaFreq/m_hwhm) + 0.5;
    }

    doublereal Lorentzian::width() {
        return 2.0*m_hwhm;
    }

    Gaussian::Gaussian(doublereal FWHM) {
        m_width = FWHM;
        m_sigma = 0.5*FWHM / sqrt(2.0 * log(2.0));
        m_sigma2 = m_sigma*m_sigma;
    }

    doublereal Gaussian::profile(doublereal deltaFreq) {
        return 1.0/(m_sigma*SqrtTwo*SqrtPi) * 
            exp(-deltaFreq*deltaFreq/(2.0*m_sigma2));
    }
        
    doublereal Gaussian::cumulative(doublereal deltaFreq) {
        return 0.5*(1.0 + erf(deltaFreq/(m_sigma*SqrtTwo)));
    }

    doublereal Gaussian::width() {
        return m_width;
    }

}
