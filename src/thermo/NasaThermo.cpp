/*!
 * @file NasaThermo.cpp Implementation of class Cantera::NasaThermo
 */
#include "NasaThermo.h"

#include "cantera/base/utilities.h"
#include "cantera/numerics/DenseMatrix.h"
#include "cantera/numerics/ctlapack.h"

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

    if (!m_allow_discontinuities) {
        doublereal maxError = checkContinuity(name, tmid, &clow[0], &chigh[0]);
        if (maxError > 1e-6) {
            fixDiscontinuities(tlow, tmid, thigh, &clow[0], &chigh[0]);
            AssertThrowMsg(checkContinuity(name, tmid, &clow[0], &chigh[0]) < 1e-12,
                   "NasaThermo::install", "Polynomials still not continuous");
        }
    }

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
    warn_deprecated("NasaThermo::reportParams");
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

doublereal NasaThermo::cp_R(double t, const doublereal* c)
{
    return poly4(t, c+2);
}

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

doublereal NasaThermo::checkContinuity(const std::string& name, double tmid,
                                       doublereal* clow, doublereal* chigh)
{
    // heat capacity
    doublereal cplow = cp_R(tmid, clow);
    doublereal cphigh = cp_R(tmid, chigh);
    doublereal delta = cplow - cphigh;
    doublereal maxError = abs(delta);
    if (fabs(delta/(fabs(cplow)+1.0E-4)) > 0.001) {
        writelog("\n\n**** WARNING ****\nFor species "+name+
                 ", discontinuity in cp/R detected at Tmid = "
                 +fp2str(tmid)+"\n");
        writelog("\tValue computed using low-temperature polynomial:  "
                 +fp2str(cplow)+".\n");
        writelog("\tValue computed using high-temperature polynomial: "
                 +fp2str(cphigh)+".\n");
    }

    // enthalpy
    doublereal hrtlow = enthalpy_RT(tmid, clow);
    doublereal hrthigh = enthalpy_RT(tmid, chigh);
    delta = hrtlow - hrthigh;
    maxError = std::max(std::abs(delta), maxError);
    if (fabs(delta/(fabs(hrtlow)+cplow*tmid)) > 0.001) {
        writelog("\n\n**** WARNING ****\nFor species "+name+
                 ", discontinuity in h/RT detected at Tmid = "
                 +fp2str(tmid)+"\n");
        writelog("\tValue computed using low-temperature polynomial:  "
                 +fp2str(hrtlow)+".\n");
        writelog("\tValue computed using high-temperature polynomial: "
                 +fp2str(hrthigh)+".\n");
    }

    // entropy
    doublereal srlow = entropy_R(tmid, clow);
    doublereal srhigh = entropy_R(tmid, chigh);
    delta = srlow - srhigh;
    maxError = std::max(std::abs(delta), maxError);
    if (fabs(delta/(fabs(srlow)+cplow)) > 0.001) {
        writelog("\n\n**** WARNING ****\nFor species "+name+
                 ", discontinuity in s/R detected at Tmid = "
                 +fp2str(tmid)+"\n");
        writelog("\tValue computed using low-temperature polynomial:  "
                 +fp2str(srlow)+".\n");
        writelog("\tValue computed using high-temperature polynomial: "
                 +fp2str(srhigh)+".\n");
    }

    return maxError;
}

void NasaThermo::fixDiscontinuities(doublereal Tlow, doublereal Tmid,
                                    doublereal Thigh, doublereal* clow,
                                    doublereal* chigh)
{
    // The thermodynamic parameters can be written in terms nondimensionalized
    // coefficients A[i] and the nondimensional temperature t = T/Tmid as:
    //
    //     C_low(t) = A[0] + A[i] * t**i
    //     H_low(t) = A[0] + A[i] / (i+1) * t**i + A[5] / t
    //     S_low(t) = A[0]*ln(t) + A[i] / i * t**i + A[6]
    //
    // where the implicit sum is over the range 1 <= i <= 4 and the
    // nondimensional coefficients are related to the dimensional coefficients
    // a[i] by:
    //
    //     A[0] = a[0]
    //     A[i] = Tmid**i * a[i], 1 <= i <= 4
    //     A[5] = a[5] / Tmid
    //     A[6] = a[6] + a[0] * ln(Tmid)
    //
    // and corresponding relationships hold for the high-temperature
    // polynomial coefficients B[i]. This nondimensionalization is necessary
    // in order for the resulting matrix to be well-conditioned.
    //
    // The requirement that C_low(1) = C_high(1) is satisfied by:
    //
    //     B[0] = A[0] + (A[i] - B[i])
    //     C_high(t) = A[0] + (A[i] + B[i] * t**i - 1)
    //
    // The requirement that H_low(1) = H_high(1) is satisfied by:
    //
    //     B[5] = A[5] + (i / (i+1) * (B[i] - A[i]))
    //     H_high(t) = A[0] + A[5] / t + (1 - i / (i+1) / t) * A[i] +
    //                 (t**i / (i+1) - 1 + i / (i+1) / t) * B[i]
    //
    // The requirement that S_low(1) = S_high(1) is satisfied by:
    //
    //    B[6] = A[6] + (A[i] - B[i]) / i
    //    S_high(t) = A[0] * ln(t) + A[6] + (ln(t) + 1 / i) * A[i] +
    //                (-ln(t) + t**i / i - 1 / i) * B[i]

    // Formulate a linear least squares problem for the nondimensionalized
    // coefficients. In the system of equations M*x = b:
    // - each row of M consists of the factors in one of the above equations
    //   for C_low, H_high, etc. evaluated at some temperature between Tlow
    //   and Thigh
    // - x is a vector of the 11 independent coefficients (A[0] through A[6]
    //   and B[1] through B[4])
    // - B is a vector of the corresponding value of C, H, or S computed using
    //   the original polynomial.

    const size_t nTemps = 12;
    const size_t nCols = 11; // number of independent coefficients
    const size_t nRows = 3*nTemps; // Evaluate C, H, and S at each temperature
    DenseMatrix M(nRows, nCols, 0.0);
    vector_fp b(nRows);
    doublereal sqrtDeltaT = sqrt(Thigh) - sqrt(Tlow);
    vector_fp tpow(5);
    for (size_t j = 0; j < nTemps; j++) {
        double T = pow(sqrt(Tlow) + sqrtDeltaT * j / (nTemps - 1.0), 2);
        double t = T / Tmid; // non-dimensionalized temperature
        double logt = std::log(t);
        size_t n = 3 * j; // row index
        for (int i = 1; i <= 4; i++) {
            tpow[i] = pow(t, i);
        }

        // row n: Cp/R
        // row n+1: H/RT
        // row n+2: S/R
        // columns 0 through 6 are for the low-T coefficients
        // columns 7 through 10 are for the independent high-T coefficients
        M(n, 0) = 1.0;
        M(n+1,0) = 1.0;
        M(n+2,0) = logt;
        M(n+1,5) = 1.0 / t;
        M(n+2,6) = 1.0;
        if (t <= 1.0) {
            for (int i = 1; i <= 4; i++) {
                M(n,i) = tpow[i];
                M(n+1,i) = tpow[i] / (i+1);
                M(n+2,i) = tpow[i] / i;
            }
            b[n] = cp_R(T, clow);
            b[n+1] = enthalpy_RT(T, clow);
            b[n+2] = entropy_R(T, clow);
        } else {
            for (int i = 1; i <= 4; i++) {
                M(n,i) = 1.0;
                M(n,i+6) = tpow[i] - 1.0;
                M(n+1,i) = 1 - i / ((i + 1.0) * t);
                M(n+1,i+6) = -1 + tpow[i] / (i+1) + i / ((i+1) * t);
                M(n+2,i) = logt + 1.0 / i;
                M(n+2,i+6) = -logt + (tpow[i] - 1.0) / i;
            }
            b[n] = cp_R(T, chigh);
            b[n+1] = enthalpy_RT(T, chigh);
            b[n+2] = entropy_R(T, chigh);
        }
    }

    // Solve the least squares problem
    vector_fp sigma(nRows);
    size_t rank;
    int info;
    vector_fp work(1);
    int lwork = -1;
    // First get the desired size of the work array
    ct_dgelss(nRows, nCols, 1, &M(0,0), nRows, &b[0], nRows,
              &sigma[0], -1, rank, &work[0], lwork, info);
    work.resize(work[0]);
    lwork = work[0];
    ct_dgelss(nRows, nCols, 1, &M(0,0), nRows, &b[0], nRows,
              &sigma[0], -1, rank, &work[0], lwork, info);

    AssertTrace(info == 0);
    AssertTrace(rank == nCols);
    AssertTrace(sigma[0] / sigma[10] < 1e20); // condition number

    // Compute the full set of nondimensionalized coefficients
    // (dgelss returns the solution of M*x = b in b).

    // Note that clow and chigh store the coefficients in the order:
    // clow = [a[5], a[6], a[0], a[1], a[2], a[3], a[4]]
    clow[2] = chigh[2] = b[0];
    clow[0] = chigh[0] = b[5];
    clow[1] = chigh[1] = b[6];
    for (int i = 1; i <= 4; i++) {
        clow[2+i] = b[i];
        chigh[2+i] = b[6+i];
        chigh[2] += clow[2+i] - chigh[2+i];
        chigh[0] += i / (i + 1.0) * (chigh[2+i] - clow[2+i]);
        chigh[1] += (clow[2+i] - chigh[2+i]) / i;
    }

    // redimensionalize
    for (int i = 1; i <= 4; i++) {
        clow[2+i] /= pow(Tmid, i);
        chigh[2+i] /= pow(Tmid, i);
    }
    clow[0] *= Tmid;
    chigh[0] *= Tmid;
    clow[1] -= clow[2] * std::log(Tmid);
    chigh[1] -= chigh[2] * std::log(Tmid);
}

}
