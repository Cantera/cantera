/*!
 * @file NasaThermo.cpp Implementation of class Cantera::NasaThermo
 */
#include "NasaThermo.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

NasaThermo::NasaThermo() :
    ID(NASA),
    m_tlow_max(0.0),
    m_thigh_min(1.e30),
    m_p0(-1.0),
    m_ngroups(0) {
    m_t.resize(6);
    }

NasaThermo::NasaThermo(const NasaThermo& right) :
    ID(NASA),
    m_tlow_max(0.0),
    m_thigh_min(1.e30),
    m_p0(-1.0),
    m_ngroups(0) {
    *this = operator=(right);
}

NasaThermo& NasaThermo::operator=(const NasaThermo& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    m_high           = right.m_high;
    m_low            = right.m_low;
    m_index          = right.m_index;
    m_tmid           = right.m_tmid;
    m_tlow_max       = right.m_tlow_max;
    m_thigh_min      = right.m_thigh_min;
    m_tlow           = right.m_tlow;
    m_thigh          = right.m_thigh;
    m_p0             = right.m_p0;
    m_ngroups        = right.m_ngroups;
    m_t              = right.m_t;
    m_group_map      = right.m_group_map;
    m_posInGroup_map = right.m_posInGroup_map;
    m_name           = right.m_name;

    return *this;
}

void NasaThermo::install(const std::string& name, size_t index, int type,
                         const doublereal* c,
                         doublereal min_temp, doublereal max_temp,
                         doublereal ref_pressure)
{
    m_name[index] = name;
    int imid = int(c[0]);       // midpoint temp converted to integer
    int igrp = m_index[imid];   // has this value been seen before?
    if (igrp == 0) {            // if not, prepare new group
        std::vector<NasaPoly1> v;
        m_high.push_back(v);
        m_low.push_back(v);
        m_tmid.push_back(c[0]);
        m_index[imid] = igrp = static_cast<int>(m_high.size());
        m_ngroups++;
    }

    m_group_map[index] = igrp;
    m_posInGroup_map[index] = (int) m_low[igrp-1].size();

    doublereal tlow  = min_temp;
    doublereal tmid  = c[0];
    doublereal thigh = max_temp;

    vector_fp chigh(c+8, c+15);
    vector_fp clow(c+1, c+8);

    ensureContinuity(name, tmid, &clow[0], &chigh[0]);

    m_high[igrp-1].push_back(NasaPoly1(index, tmid, thigh,
                                       ref_pressure, &chigh[0]));
    m_low[igrp-1].push_back(NasaPoly1(index, tlow, tmid,
                                      ref_pressure, &clow[0]));

    if (tlow > m_tlow_max) {
        m_tlow_max = tlow;
    }
    if (thigh < m_thigh_min) {
        m_thigh_min = thigh;
    }
    if (m_tlow.size() < index + 1) {
        m_tlow.resize(index + 1,  tlow);
        m_thigh.resize(index + 1, thigh);
    }
    m_tlow[index] = tlow;
    m_thigh[index] = thigh;
    if (m_p0 < 0.0) {
        m_p0 = ref_pressure;
    } else if (fabs(m_p0 - ref_pressure) > 0.1) {
        std::string logmsg =  " ERROR NasaThermo: New Species, " + name +  ", has a different reference pressure, "
                              + fp2str(ref_pressure) + ", than existing reference pressure, " + fp2str(m_p0) + "\n";
        writelog(logmsg);
        logmsg = "                  This is now a fatal error\n";
        writelog(logmsg);
        throw CanteraError("install()", "species have different reference pressures");
    }
    m_p0 = ref_pressure;
}

void NasaThermo::update_one(size_t k, doublereal t, doublereal* cp_R,
                            doublereal* h_RT, doublereal* s_R) const
{
    m_t[0] = t;
    m_t[1] = t*t;
    m_t[2] = m_t[1]*t;
    m_t[3] = m_t[2]*t;
    m_t[4] = 1.0/t;
    m_t[5] = log(t);

    size_t grp = m_group_map[k];
    size_t pos = m_posInGroup_map[k];
    const std::vector<NasaPoly1> &mlg = m_low[grp-1];
    const NasaPoly1* nlow = &(mlg[pos]);

    doublereal tmid = nlow->maxTemp();
    if (t < tmid) {
        nlow->updateProperties(&m_t[0], cp_R, h_RT, s_R);
    } else {
        const std::vector<NasaPoly1> &mhg = m_high[grp-1];
        const NasaPoly1* nhigh = &(mhg[pos]);
        nhigh->updateProperties(&m_t[0], cp_R, h_RT, s_R);
    }
}

void NasaThermo::update(doublereal t, doublereal* cp_R,
                        doublereal* h_RT, doublereal* s_R) const
{
    int i;

    // load functions of temperature into m_t vector
    m_t[0] = t;
    m_t[1] = t*t;
    m_t[2] = m_t[1]*t;
    m_t[3] = m_t[2]*t;
    m_t[4] = 1.0/t;
    m_t[5] = log(t);

    // iterate over the groups
    std::vector<NasaPoly1>::const_iterator _begin, _end;
    for (i = 0; i != m_ngroups; i++) {
        if (t > m_tmid[i]) {
            _begin  = m_high[i].begin();
            _end    = m_high[i].end();
        } else {
            _begin  = m_low[i].begin();
            _end    = m_low[i].end();
        }
        for (; _begin != _end; ++_begin) {
            _begin->updateProperties(&m_t[0], cp_R, h_RT, s_R);
        }
    }
}

void NasaThermo::reportParams(size_t index, int& type,
                              doublereal* const c,
                              doublereal& minTemp,
                              doublereal& maxTemp,
                              doublereal& refPressure) const
{
    type = reportType(index);
    if (type == NASA) {
        size_t grp = m_group_map[index];
        size_t pos = m_posInGroup_map[index];
        const std::vector<NasaPoly1> &mlg = m_low[grp-1];
        const std::vector<NasaPoly1> &mhg = m_high[grp-1];
        const NasaPoly1* lowPoly  = &(mlg[pos]);
        const NasaPoly1* highPoly = &(mhg[pos]);
        int itype = NASA;
        doublereal tmid = lowPoly->maxTemp();
        c[0] = tmid;
        size_t n;
        double ttemp;
        lowPoly->reportParameters(n, itype, minTemp, ttemp, refPressure,
                                  c + 1);
        if (n != index) {
            throw CanteraError("  ", "confused");
        }
        if (itype != NASA1) {
            throw CanteraError("  ", "confused");
        }
        highPoly->reportParameters(n, itype, ttemp, maxTemp, refPressure,
                                   c + 8);
        if (n != index) {
            throw CanteraError("  ", "confused");
        }
        if (itype != NASA1) {
            throw CanteraError("  ", "confused");
        }
    } else {
        throw CanteraError(" ", "confused");
    }
}

#ifdef H298MODIFY_CAPABILITY
doublereal NasaThermo::reportOneHf298(const int k) const
{
    int grp = m_group_map[k];
    int pos = m_posInGroup_map[k];
    const std::vector<NasaPoly1> &mlg = m_low[grp-1];
    const NasaPoly1* nlow = &(mlg[pos]);
    doublereal tmid = nlow->maxTemp();
    double h;
    if (298.15 <= tmid) {
        h = nlow->reportHf298(0);
    } else {
        const std::vector<NasaPoly1> &mhg = m_high[grp-1];
        const NasaPoly1* nhigh = &(mhg[pos]);
        h = nhigh->reportHf298(0);
    }
    return h;
}

void NasaThermo::modifyOneHf298(const int k, const doublereal Hf298New)
{
    int grp = m_group_map[k];
    int pos = m_posInGroup_map[k];
    std::vector<NasaPoly1> &mlg = m_low[grp-1];
    NasaPoly1* nlow = &(mlg[pos]);
    std::vector<NasaPoly1> &mhg = m_high[grp-1];
    NasaPoly1* nhigh = &(mhg[pos]);
    doublereal tmid = nlow->maxTemp();

    double hnow = reportOneHf298(k);
    double delH =  Hf298New - hnow;
    if (298.15 <= tmid) {
        nlow->modifyOneHf298(k, Hf298New);
        double h = nhigh->reportHf298(0);
        double hnew = h + delH;
        nhigh->modifyOneHf298(k, hnew);
    } else {
        nhigh->modifyOneHf298(k, Hf298New);
        double h = nlow->reportHf298(0);
        double hnew = h + delH;
        nlow->modifyOneHf298(k, hnew);
    }
}
#endif

doublereal NasaThermo::enthalpy_RT(double t, const doublereal* c) {
    return c[2] + 0.5*c[3]*t + OneThird*c[4]*t*t
           + 0.25*c[5]*t*t*t + 0.2*c[6]*t*t*t*t
           + c[0]/t;
}

doublereal NasaThermo::entropy_R(double t, const doublereal* c) {
    return c[2]*log(t) + c[3]*t + 0.5*c[4]*t*t
           + OneThird*c[5]*t*t*t + 0.25*c[6]*t*t*t*t
           + c[1];
}

void NasaThermo::ensureContinuity(const std::string& name, double tmid,
                                  doublereal* clow, doublereal* chigh)
{
    // heat capacity
    doublereal cplow = poly4(tmid, clow + 2);
    doublereal cphigh = poly4(tmid, chigh + 2);
    doublereal delta = cplow - cphigh;
    if (fabs(delta/(fabs(cplow)+1.0E-4)) > 0.001) {
        writelog("\n\n**** WARNING ****\nFor species "+name+
                 ", discontinuity in cp/R detected at Tmid = "
                 +fp2str(tmid)+"\n");
        writelog("\tValue computed using low-temperature polynomial:  "
                 +fp2str(cplow)+".\n");
        writelog("\tValue computed using high-temperature polynomial: "
                 +fp2str(cphigh)+".\n");
    }

    // Adjust coefficients to eliminate any discontinuity
    chigh[2] += 0.5 * delta;
    clow[2] -= 0.5 * delta;

    AssertThrowMsg(std::abs(poly4(tmid, clow+2) - poly4(tmid, chigh+2)) < 1e-12,
                   "NasaThermo::ensureContinuity", "Cp/R does not match");

    // enthalpy
    doublereal hrtlow = enthalpy_RT(tmid, clow);
    doublereal hrthigh = enthalpy_RT(tmid, chigh);
    delta = hrtlow - hrthigh;
    if (fabs(delta/(fabs(hrtlow)+cplow*tmid)) > 0.001) {
        writelog("\n\n**** WARNING ****\nFor species "+name+
                 ", discontinuity in h/RT detected at Tmid = "
                 +fp2str(tmid)+"\n");
        writelog("\tValue computed using low-temperature polynomial:  "
                 +fp2str(hrtlow)+".\n");
        writelog("\tValue computed using high-temperature polynomial: "
                 +fp2str(hrthigh)+".\n");
    }

    // Adjust coefficients to eliminate any discontinuity
    chigh[0] += 0.5 * delta * tmid;
    clow[0] -= 0.5 * delta * tmid;

    AssertThrowMsg(std::abs(enthalpy_RT(tmid, clow) -
                            enthalpy_RT(tmid, chigh)) < 1e-12,
                   "NasaThermo::ensureContinuity", "H/RT does not match");

    // entropy
    doublereal srlow = entropy_R(tmid, clow);
    doublereal srhigh = entropy_R(tmid, chigh);
    delta = srlow - srhigh;
    if (fabs(delta/(fabs(srlow)+cplow)) > 0.001) {
        writelog("\n\n**** WARNING ****\nFor species "+name+
                 ", discontinuity in s/R detected at Tmid = "
                 +fp2str(tmid)+"\n");
        writelog("\tValue computed using low-temperature polynomial:  "
                 +fp2str(srlow)+".\n");
        writelog("\tValue computed using high-temperature polynomial: "
                 +fp2str(srhigh)+".\n");
    }

    // Adjust coefficients to eliminate any discontinuity
    chigh[1] += 0.5 * delta;
    clow[1] -= 0.5 * delta;

    AssertThrowMsg(std::abs(entropy_R(tmid, clow) -
                            entropy_R(tmid, chigh)) < 1e-12,
                   "NasaThermo::ensureContinuity", "S/R does not match");
}

}
