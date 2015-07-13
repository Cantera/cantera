#include "cantera/thermo/NasaPoly2.h"
#include "cantera/base/global.h"
#include "cantera/base/stringUtils.h"

namespace Cantera {

void NasaPoly2::validate(const std::string& name)
{
    double cp_low, h_low, s_low;
    double cp_high, h_high, s_high;
    mnp_low.updatePropertiesTemp(m_midT, &cp_low, &h_low, &s_low);
    mnp_high.updatePropertiesTemp(m_midT, &cp_high, &h_high, &s_high);

    double delta = cp_low - cp_high;
    if (fabs(delta/(fabs(cp_low)+1.0E-4)) > 0.001) {
        writelog("\n\n**** WARNING ****\nFor species "+name+
                 ", discontinuity in cp/R detected at Tmid = "
                 +fp2str(m_midT)+"\n");
        writelog("\tValue computed using low-temperature polynomial:  "
                 +fp2str(cp_low)+".\n");
        writelog("\tValue computed using high-temperature polynomial: "
                 +fp2str(cp_high)+".\n");
    }

    // enthalpy
    delta = h_low - h_high;
    if (fabs(delta/(fabs(h_low)+cp_low*m_midT)) > 0.001) {
        writelog("\n\n**** WARNING ****\nFor species "+name+
                 ", discontinuity in h/RT detected at Tmid = "
                 +fp2str(m_midT)+"\n");
        writelog("\tValue computed using low-temperature polynomial:  "
                 +fp2str(h_low)+".\n");
        writelog("\tValue computed using high-temperature polynomial: "
                 +fp2str(h_high)+".\n");
    }

    // entropy
    delta = s_low - s_high;
    if (fabs(delta/(fabs(s_low)+cp_low)) > 0.001) {
        writelog("\n\n**** WARNING ****\nFor species "+name+
                 ", discontinuity in s/R detected at Tmid = "
                 +fp2str(m_midT)+"\n");
        writelog("\tValue computed using low-temperature polynomial:  "
                 +fp2str(s_low)+".\n");
        writelog("\tValue computed using high-temperature polynomial: "
                 +fp2str(s_high)+".\n");
    }
}

}
