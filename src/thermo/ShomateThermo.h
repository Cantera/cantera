/**
 * @file ShomateThermo.h
 *   Header for the 2 regions Shomate polynomial
 *   for multiple species in a phase, derived from the
 *   \link Cantera::SpeciesThermo SpeciesThermo\endlink base class (see \ref mgrsrefcalc and
 *   \link Cantera::ShomateThermo ShomateThermo\endlink).
 */
// Copyright 2001  California Institute of Technology

#ifndef CT_SHOMATETHERMO_H
#define CT_SHOMATETHERMO_H

#include "cantera/thermo/SpeciesThermoMgr.h"
#include "cantera/thermo/ShomatePoly.h"
#include "cantera/base/global.h"
#include "cantera/base/utilities.h"

namespace Cantera
{
//! A species thermodynamic property manager for the Shomate polynomial parameterization.
/*!
 *  This is the parameterization used
 *  in the NIST Chemistry WebBook (http://webbook.nist.gov/chemistry)
 * The parameterization assumes there are two temperature regions
 * each with its own Shomate polynomial representation, for each
 * species in the phase.
 *
 * \f[
 * \tilde{c}_p^0(T) = A + B t + C t^2 + D t^3 + \frac{E}{t^2}
 * \f]
 * \f[
 * \tilde{h}^0(T) = A t + \frac{B t^2}{2} + \frac{C t^3}{3}
 + \frac{D t^4}{4}  - \frac{E}{t}  + F.
 * \f]
 * \f[
 * \tilde{s}^0(T) = A\ln t + B t + \frac{C t^2}{2}
 + \frac{D t^3}{3} - \frac{E}{2t^2}  + G.
 * \f]
 *
 * In the above expressions, the thermodynamic polynomials are expressed
 * in dimensional units, but the temperature,\f$ t \f$, is divided by 1000. The
 * following dimensions are assumed in the above expressions:
 *
 *    - \f$ \tilde{c}_p^0(T)\f$ = Heat Capacity (J/gmol*K)
 *    - \f$ \tilde{h}^0(T) \f$ = standard Enthalpy (kJ/gmol)
 *    - \f$ \tilde{s}^0(T) \f$= standard Entropy (J/gmol*K)
 *    - \f$ t \f$= temperature (K) / 1000.
 *
 *  Note, the polynomial data (i.e., A, ... , G) is entered in dimensional form.
 *
 *  This is in contrast to the NASA database polynomials which are entered in
 *  nondimensional form (i.e., NASA parameterizes C_p/R, while Shomate
 *  parameterizes C_p assuming units of J/gmol*K - and kJ/gmol*K for H).
 *  Note, also that the H - H_298.15 equation has units of kJ/gmol, because of
 *  the implicit integration of (t = T 1000), which provides a
 *  multiplier of 1000 to the Enthalpy equation.
 *
 * @deprecated To be removed after Cantera 2.2. Use GeneralSpeciesThermo instead.
 * @ingroup mgrsrefcalc
 */
class ShomateThermo : public SpeciesThermo
{
public:
    //! Initialized to the type of parameterization
    /*!
     * Note, this value is used in some template functions
     */
    const int ID;

    //! constructor
    ShomateThermo() :
        ID(SHOMATE),
        m_tlow_max(0.0),
        m_thigh_min(1.e30),
        m_p0(-1.0),
        m_ngroups(0) {
        warn_deprecated("class ShomateThermo", "To be removed after "
            "Cantera 2.2. Use GeneralSpeciesThermo instead.");
        m_t.resize(7);
    }

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    ShomateThermo(const ShomateThermo& right) :
        ID(SHOMATE),
        m_tlow_max(0.0),
        m_thigh_min(1.e30),
        m_p0(-1.0),
        m_ngroups(0) {
        *this = right;
    }

    //! Assignment Operator
    /*!
     * @param right Object to be copied
     */
    ShomateThermo& operator=(const ShomateThermo& right) {
        if (&right == this) {
            return *this;
        }

        SpeciesThermo::operator=(right);
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

        return *this;
    }

    virtual SpeciesThermo* duplMyselfAsSpeciesThermo() const {
        return new ShomateThermo(*this);
    }

    //! Install a new species thermodynamic property
    //! parameterization for one species using Shomate polynomials
    /*!
     *  Two temperature regions are assumed.
     *
     * @param name      Name of the species
     * @param index     Species index
     * @param type      int flag specifying the type of parameterization to be
     *                  installed.
     * @param c         Vector of coefficients for the parameterization.
     *                  There are 15 coefficients for the 2-zone Shomate polynomial.
     *                  The first coefficient is the value of Tmid. The next 7
     *                  coefficients are the low temperature range Shomate coefficients.
     *                  The last 7 are the high temperature range Shomate coefficients.
     *
     * @param minTemp  minimum temperature for which this parameterization
     *                 is valid.
     * @param maxTemp  maximum temperature for which this parameterization
     *                 is valid.
     * @param refPressure standard-state pressure for this
     *                    parameterization.
     *
     * @see ShomatePoly
     * @see ShomatePoly2
     */
    virtual void install(const std::string& name, size_t index, int type,
                         const doublereal* c,
                         doublereal minTemp, doublereal maxTemp,
                         doublereal refPressure) {
        if (type != SHOMATE) {
            throw CanteraError("ShomateThermo::install",
                               "Incompatible thermo parameterization: Got " +
                               int2str(type) + " but " + int2str(SHOMATE) +
                               " was expected.");
        }
        int imid = int(c[0]);       // midpoint temp converted to integer
        int igrp = m_index[imid];   // has this value been seen before?
        if (igrp == 0) {            // if not, prepare new group
            std::vector<ShomatePoly> v;
            m_high.push_back(v);
            m_low.push_back(v);
            m_tmid.push_back(c[0]);
            m_index[imid] = igrp = static_cast<int>(m_high.size());
            m_ngroups++;
        }
        m_group_map[index] = igrp;
        m_posInGroup_map[index] = (int) m_low[igrp-1].size();
        doublereal tlow  = minTemp;
        doublereal tmid  = c[0];
        doublereal thigh = maxTemp;

        const doublereal* clow = c + 1;
        const doublereal* chigh = c + 8;
        m_high[igrp-1].push_back(ShomatePoly(index, tmid, thigh,
                                             refPressure, chigh));
        m_low[igrp-1].push_back(ShomatePoly(index, tlow, tmid,
                                            refPressure, clow));
        m_tlow_max = std::max(m_tlow_max, tlow);
        m_thigh_min = std::min(m_thigh_min, thigh);
        if (m_tlow.size() < index + 1) {
            m_tlow.resize(index + 1,  tlow);
            m_thigh.resize(index + 1, thigh);
        }
        m_tlow[index] = tlow;
        m_thigh[index] = thigh;

        if (m_p0 < 0.0) {
            m_p0 = refPressure;
        } else if (fabs(m_p0 - refPressure) > 0.1) {
            std::string logmsg =  " ERROR ShomateThermo: New Species, " + name
                                  +  ", has a different reference pressure, "
                                  + fp2str(refPressure) + ", than existing reference pressure, " + fp2str(m_p0) + "\n";
            writelog(logmsg);
            logmsg = "                  This is now a fatal error\n";
            writelog(logmsg);
            throw CanteraError("install()", "Species have different reference pressures");
        }
        m_p0 = refPressure;
        markInstalled(index);
    }

    virtual void install_STIT(size_t index,
                              shared_ptr<SpeciesThermoInterpType> stit_ptr) {
        throw CanteraError("install_STIT", "not implemented");
    }

    //! Like update(), but only updates the single species k.
    /*!
     * @param k       species index
     * @param t       Temperature (Kelvin)
     * @param cp_R    Vector of Dimensionless heat capacities. (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies. (length m_kk).
     * @param s_R     Vector of Dimensionless entropies. (length m_kk).
     */
    virtual void update_one(size_t k, doublereal t, doublereal* cp_R,
                            doublereal* h_RT, doublereal* s_R) const {
        doublereal tt = 1.e-3*t;
        m_t[0] = tt;
        m_t[1] = tt*tt;
        m_t[2] = m_t[1]*tt;
        m_t[3] = 1.0/m_t[1];
        m_t[4] = log(tt);
        m_t[5] = 1.0/GasConstant;
        m_t[6] = 1.0/(GasConstant * t);

        size_t grp = getValue(m_group_map, k);
        size_t pos = getValue(m_posInGroup_map, k);
        const std::vector<ShomatePoly> &mlg = m_low[grp-1];
        const ShomatePoly* nlow = &(mlg[pos]);

        if (t < nlow->maxTemp()) {
            nlow->updateProperties(&m_t[0], cp_R, h_RT, s_R);
        } else {
            const std::vector<ShomatePoly> &mhg = m_high[grp-1];
            const ShomatePoly* nhigh = &(mhg[pos]);
            nhigh->updateProperties(&m_t[0], cp_R, h_RT, s_R);
        }
    }

    virtual void update(doublereal t, doublereal* cp_R,
                        doublereal* h_RT, doublereal* s_R) const {
        doublereal tt = 1.e-3*t;
        m_t[0] = tt;
        m_t[1] = tt*tt;
        m_t[2] = m_t[1]*tt;
        m_t[3] = 1.0/m_t[1];
        m_t[4] = log(tt);
        m_t[5] = 1.0/GasConstant;
        m_t[6] = 1.0/(GasConstant * t);

        std::vector<ShomatePoly>::const_iterator _begin, _end;
        for (int i = 0; i != m_ngroups; i++) {
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

    virtual doublereal minTemp(size_t k=npos) const {
        if (k == npos) {
            return m_tlow_max;
        } else {
            return m_tlow[k];
        }
    }

    virtual doublereal maxTemp(size_t k=npos) const {
        if (k == npos) {
            return m_thigh_min;
        } else {
            return m_thigh[k];
        }
    }

    virtual doublereal refPressure(size_t k=npos) const {
        return m_p0;
    }

    virtual int reportType(size_t index) const {
        return SHOMATE;
    }

    virtual void reportParams(size_t index, int& type,
                              doublereal* const c,
                              doublereal& minTemp,
                              doublereal& maxTemp,
                              doublereal& refPressure) const {
        type = reportType(index);
        if (type == SHOMATE) {
            size_t grp = getValue(m_group_map, index);
            size_t pos = getValue(m_posInGroup_map, index);
            int itype = SHOMATE;
            const std::vector<ShomatePoly> &mlg = m_low[grp-1];
            const std::vector<ShomatePoly> &mhg = m_high[grp-1];
            const ShomatePoly* lowPoly  = &(mlg[pos]);
            const ShomatePoly* highPoly = &(mhg[pos]);
            doublereal tmid = lowPoly->maxTemp();
            c[0] = tmid;
            size_t n;
            double ttemp;
            lowPoly->reportParameters(n, itype, minTemp, ttemp, refPressure,
                                      c + 1);
            if (n != index) {
                throw CanteraError("ShomateThermo::reportParams",
                                   "Index mismatch in low-T polynomial");
            }
            if (itype != SHOMATE && itype != SHOMATE1) {
                throw CanteraError("ShomateThermo::reportParams",
                                   "Thermo type mismatch in low-T polynomial");
            }
            highPoly->reportParameters(n, itype,  ttemp, maxTemp,
                                       refPressure, c + 8);
            if (n != index) {
                throw CanteraError("ShomateThermo::reportParams",
                                   "Index mismatch in high-T polynomial");
            }
            if (itype != SHOMATE && itype != SHOMATE1) {
                throw CanteraError("ShomateThermo::reportParams",
                                   "Thermo type mismatch in high-T polynomial");
            }
        } else {
            throw CanteraError("ShomateThermo::reportParams", "Thermo type mismatch");
        }
    }

    virtual doublereal reportOneHf298(const size_t k) const {
        size_t grp = getValue(m_group_map, k);
        size_t pos = getValue(m_posInGroup_map, k);
        const ShomatePoly& nlow = m_low[grp-1][pos];

        if (nlow.maxTemp() > 298.15) {
            return nlow.reportHf298();
        } else {
            const ShomatePoly& nhigh = m_high[grp-1][pos];
            return nhigh.reportHf298();
        }
    }

    virtual void modifyOneHf298(const size_t k, const doublereal Hf298New) {

        size_t grp = m_group_map[k];
        size_t pos = m_posInGroup_map[k];
        ShomatePoly& nlow = m_low[grp-1][pos];
        ShomatePoly& nhigh = m_high[grp-1][pos];

        double hnow = reportOneHf298(k);
        double delH =  Hf298New - hnow;
        if (nlow.maxTemp() > 298.15) {
            nlow.modifyOneHf298(k, Hf298New);
            double h = nhigh.reportHf298(0);
            double hnew = h + delH;
            nhigh.modifyOneHf298(k, hnew);
        } else {
            nhigh.modifyOneHf298(k, Hf298New);
            double h = nlow.reportHf298(0);
            double hnew = h + delH;
            nlow.modifyOneHf298(k, hnew);
        }

    }

protected:
    //! Vector of vector of NasaPoly1's for the high temp region.
    /*!
     * This is the high temp region representation. The first Length is equal
     * to the number of groups. The second vector is equal to the number of
     * species in that particular group.
     */
    std::vector<std::vector<ShomatePoly> > m_high;

    //! Vector of vector of NasaPoly1's for the low temp region.
    /*!
     * This is the low temp region representation. The first Length is equal
     * to the number of groups. The second vector is equal to the number of
     * species in that particular group.
     */
    std::vector<std::vector<ShomatePoly> > m_low;

    //! Map between the midpoint temperature, as an int, to the group number
    /*!
     * Length is equal to the number of groups. Only used in the setup.
     */
    std::map<int, int>              m_index;

    //! Vector of log temperature limits
    /*!
     * Length is equal to the number of groups.
     */
    vector_fp                  m_tmid;

    //! Maximum value of the low temperature limit
    doublereal                 m_tlow_max;

    //! Minimum value of the high temperature limit
    doublereal                 m_thigh_min;

    //! Vector of low temperature limits (species index)
    /*!
     * Length is equal to number of species
     */
    vector_fp                  m_tlow;

    //! Vector of low temperature limits (species index)
    /*!
     * Length is equal to number of species
     */
    vector_fp                  m_thigh;

    //! Reference pressure (Pa)
    /*!
     * all species must have the same reference pressure.
     */
    doublereal                 m_p0;

    //! number of groups
    int                        m_ngroups;

    //! Vector of temperature polynomials
    mutable vector_fp          m_t;

    /*!
     * This map takes as its index, the species index in the phase.
     * It returns the group index, where the temperature polynomials
     * for that species are stored. group indices start at 1,
     * so a decrement is always performed to access vectors.
     */
    std::map<size_t, size_t> m_group_map;

    /*!
     * This map takes as its index, the species index in the phase.
     * It returns the position index within the group, where the
     * temperature polynomials for that species are stored.
     */
    std::map<size_t, size_t> m_posInGroup_map;
};

}

#endif
