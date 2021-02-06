// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/NasaPoly2.h"
#include "cantera/base/global.h"
#include "cantera/base/stringUtils.h"

namespace Cantera {

NasaPoly2::NasaPoly2()
    : m_midT(0)
{
}

void NasaPoly2::setParameters(double Tmid, const vector_fp& low,
                              const vector_fp& high) {
    m_midT = Tmid;
    mnp_low.setMaxTemp(Tmid);
    mnp_high.setMinTemp(Tmid);
    mnp_low.setParameters(low);
    mnp_high.setParameters(high);
}

void NasaPoly2::validate(const std::string& name)
{
    if (thermo_warnings_suppressed()) {
        return;
    }

    double cp_low, h_low, s_low;
    double cp_high, h_high, s_high;
    mnp_low.updatePropertiesTemp(m_midT, &cp_low, &h_low, &s_low);
    mnp_high.updatePropertiesTemp(m_midT, &cp_high, &h_high, &s_high);

    double delta = cp_low - cp_high;
    if (fabs(delta/(fabs(cp_low)+1.0E-4)) > 0.01) {
        warn_user("NasaPoly2::validate",
            "\nFor species {}, discontinuity in cp/R detected at Tmid = {}\n"
            "\tValue computed using low-temperature polynomial:  {}\n"
            "\tValue computed using high-temperature polynomial: {}\n",
            name, m_midT, cp_low, cp_high);
    }

    // enthalpy
    delta = h_low - h_high;
    if (fabs(delta/cp_low) > 0.001) {
        warn_user("NasaPoly2::validate",
            "\nFor species {}, discontinuity in h/RT detected at Tmid = {}\n"
            "\tValue computed using low-temperature polynomial:  {}\n"
            "\tValue computed using high-temperature polynomial: {}\n",
            name, m_midT, h_low, h_high);
    }

    // entropy
    delta = s_low - s_high;
    if (fabs(delta/(fabs(s_low)+cp_low)) > 0.001) {
        warn_user("NasaPoly2::validate",
            "\nFor species {}, discontinuity in s/R detected at Tmid = {}\n"
            "\tValue computed using low-temperature polynomial:  {}\n"
            "\tValue computed using high-temperature polynomial: {}\n",
            name, m_midT, s_low, s_high);
    }
}

}
