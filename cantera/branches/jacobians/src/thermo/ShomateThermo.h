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
#include "ShomatePoly.h"
#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/Cantera.h"

namespace Cantera {

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
 *  Note, the polynomial data (i.e., A, ... , G) is entered in dimensional
 *        form.
 *
 *  This is in contrast to the NASA database polynomials which are entered in
 *  nondimensional form (i.e., NASA parameterizes C_p/R, while Shomate
 *  parameterizes C_p assuming units of J/gmol*K - and kJ/gmol*K for H).
 *  Note, also that the H - H_298.15 equation has units of kJ/gmol, because of
 *  the implicit integration of (t = T 1000), which provides a
 *  multiplier of 1000 to the Enthalpy equation.
 *
 * @ingroup mgrsrefcalc
 */
template<typename ValAndDerivType>
class ShomateThermo: public SpeciesThermo<ValAndDerivType>
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
            m_ngroups(0)
    {
        m_t.resize(7);
    }

    //! destructor
    virtual ~ShomateThermo()
    {
    }

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    ShomateThermo(const ShomateThermo<ValAndDerivType>& right) :
            ID(SHOMATE),
            m_tlow_max(0.0),
            m_thigh_min(1.e30),
            m_p0(-1.0),
            m_ngroups(0)
    {
        operator=(right);
    }

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    template<typename ValAndDerivType2>
    ShomateThermo(const ShomateThermo<ValAndDerivType2>& right) :
            ID(SHOMATE),
            m_tlow_max(0.0),
            m_thigh_min(1.e30),
            m_p0(-1.0),
            m_ngroups(0)
    {
        operator=(right);
    }

    //=============================================================================================================================
    //! Assignment Operator
    /*!
     * @param right Object to be copied
     */
    ShomateThermo& operator=(const ShomateThermo<ValAndDerivType>& right);

    template<typename ValAndDerivType2>
    ShomateThermo& operator=(const ShomateThermo<ValAndDerivType2>& right);

    //=============================================================================================================================
    //! Duplication routine for objects which inherit from
    //! %SpeciesThermo
    /*!
     *  This virtual routine can be used to duplicate %SpeciesThermo  objects
     *  inherited from %SpeciesThermo even if the application only has
     *  a pointer to %SpeciesThermo to work with.
     *  ->commented out because we first need to add copy constructors
     *   and assignment operators to all of the derived classes.
     */
    virtual SpeciesThermo<ValAndDerivType>* duplMyselfAsSpeciesThermo() const
    {
        ShomateThermo<ValAndDerivType>* st = new ShomateThermo<ValAndDerivType>(*this);
        return (SpeciesThermo<ValAndDerivType> *) st;
    }

    //! Duplication routine for objects which inherit from %SpeciesThermo
    /*!
     *  This virtual routine can be used to duplicate %SpeciesThermo  objects
     *  inherited from %SpeciesThermo even if the application only has
     *  a pointer to %SpeciesThermo to work with.
     *
     *  This routine returns a doublereal templated version of SpeciesThermo no matter
     *  what templated version the underlying class is.
     *
     *  @return Duplicated <double> version of the SpeciesThermo
     */
    virtual SpeciesThermo<doublereal>* duplMyselfAsSpeciesThermoDouble() const
    {
        ShomateThermo<doublereal>* st = new ShomateThermo<doublereal>(*this);
        return (SpeciesThermo<doublereal> *) st;
    }

    //! Install a new species thermodynamic property
    //! parameterization for one species using Shomate polynomials
    //!
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
    virtual void install(const std::string& name, size_t index, int type, const doublereal* c, doublereal minTemp,
                         doublereal maxTemp, doublereal refPressure)
    {
        int imid = int(c[0]); // midpoint temp converted to integer
        int igrp = m_index[imid]; // has this value been seen before?
        if (igrp == 0) { // if not, prepare new group
            std::vector<ShomatePoly<ValAndDerivType> > v;
            m_high.push_back(v);
            m_low.push_back(v);
            m_tmid.push_back(c[0]);
            m_index[imid] = igrp = static_cast<int>(m_high.size());
            m_ngroups++;
        }
        m_group_map[index] = igrp;
        m_posInGroup_map[index] = (int) m_low[igrp - 1].size();
        doublereal tlow = minTemp;
        doublereal tmid = c[0];
        doublereal thigh = maxTemp;

        const doublereal* clow = c + 1;
        const doublereal* chigh = c + 8;
        m_high[igrp - 1].push_back(ShomatePoly<ValAndDerivType>(index, tmid, thigh, refPressure, chigh));
        m_low[igrp - 1].push_back(ShomatePoly<ValAndDerivType>(index, tlow, tmid, refPressure, clow));
        if (tlow > m_tlow_max) {
            m_tlow_max = tlow;
        }
        if (thigh < m_thigh_min) {
            m_thigh_min = thigh;
        }

        if (m_tlow.size() < index + 1) {
            m_tlow.resize(index + 1, tlow);
            m_thigh.resize(index + 1, thigh);
        }
        m_tlow[index] = tlow;
        m_thigh[index] = thigh;

        if (m_p0 < 0.0) {
            m_p0 = refPressure;
        } else if (fabs(m_p0 - refPressure) > 0.1) {
            std::string logmsg = " ERROR ShomateThermo: New Species, " + name + ", has a different reference pressure, "
                    + fp2str(refPressure) + ", than existing reference pressure, " + fp2str(m_p0) + "\n";
            writelog(logmsg);
            logmsg = "                  This is now a fatal error\n";
            writelog(logmsg);
            throw CanteraError("install()", "Species have different reference pressures");
        }
        m_p0 = refPressure;

    }

    //! Install a new species thermodynamic property
    //! parameterization for one species.
    /*!
     * @param stit_ptr Pointer to the SpeciesThermoInterpType object
     *          This will set up the thermo for one species
     */
    virtual void install_STIT(SpeciesThermoInterpType<ValAndDerivType> * stit_ptr)
    {
        throw CanteraError("install_STIT", "not implemented");
    }

    //! Like update(), but only updates the single species k.
    /*!
     * @param k       species index
     * @param t       Temperature (Kelvin)
     * @param cp_R    Vector of Dimensionless heat capacities.
     *                (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies.
     *                (length m_kk).
     * @param s_R     Vector of Dimensionless entropies.
     *                (length m_kk).
     */
    virtual void update_one(size_t k, doublereal t, ValAndDerivType* cp_R, ValAndDerivType* h_RT, ValAndDerivType* s_R) const
    {

        doublereal tt = 1.e-3 * t;
        m_t[0] = tt;
        m_t[1] = tt * tt;
        m_t[2] = m_t[1] * tt;
        m_t[3] = 1.0 / m_t[1];
        m_t[4] = log(tt);
        m_t[5] = 1.0 / GasConstant;
        m_t[6] = 1.0 / (GasConstant * t);

        size_t grp = m_group_map[k];
        size_t pos = m_posInGroup_map[k];
        const std::vector<ShomatePoly<ValAndDerivType> > &mlg = m_low[grp - 1];
        const ShomatePoly<ValAndDerivType>* nlow = & (mlg[pos]);

        doublereal tmid = nlow->maxTemp();
        if (t < tmid) {
            nlow->updateProperties(&m_t[0], cp_R, h_RT, s_R);
        } else {
            const std::vector<ShomatePoly<ValAndDerivType> > &mhg = m_high[grp - 1];
            const ShomatePoly<ValAndDerivType>* nhigh = & (mhg[pos]);
            nhigh->updateProperties(&m_t[0], cp_R, h_RT, s_R);
        }
    }

    //! Compute the reference-state properties for all species.
    /*!
     * Given temperature T in K, this method updates the values of
     * the non-dimensional heat capacity at constant pressure,
     * enthalpy, and entropy, at the reference pressure, Pref
     * of each of the standard states.
     *
     * @param t       Temperature (Kelvin)
     * @param cp_R    Vector of Dimensionless heat capacities.
     *                (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies.
     *                (length m_kk).
     * @param s_R     Vector of Dimensionless entropies.
     *                (length m_kk).
     */
    virtual void update(doublereal t, ValAndDerivType* cp_R, ValAndDerivType* h_RT, ValAndDerivType* s_R) const
    {
        int i;

        doublereal tt = 1.e-3 * t;
        m_t[0] = tt;
        m_t[1] = tt * tt;
        m_t[2] = m_t[1] * tt;
        m_t[3] = 1.0 / m_t[1];
        m_t[4] = log(tt);
        m_t[5] = 1.0 / GasConstant;
        m_t[6] = 1.0 / (GasConstant * t);

        typename std::vector<ShomatePoly<ValAndDerivType> >::const_iterator _begin, _end;
        for (i = 0; i != m_ngroups; i++) {
            if (t > m_tmid[i]) {
                _begin = m_high[i].begin();
                _end = m_high[i].end();
            } else {
                _begin = m_low[i].begin();
                _end = m_low[i].end();
            }
            for (; _begin != _end; ++_begin) {
                _begin->updateProperties(&m_t[0], cp_R, h_RT, s_R);
            }
        }
    }

    //! Minimum temperature.
    /*!
     * If no argument is supplied, this
     * method returns the minimum temperature for which \e all
     * parameterizations are valid. If an integer index k is
     * supplied, then the value returned is the minimum
     * temperature for species k in the phase.
     *
     * @param k    Species index
     */
    virtual doublereal minTemp(size_t k = npos) const
    {
        if (k == npos) {
            return m_tlow_max;
        } else {
            return m_tlow[k];
        }
    }

    //! Maximum temperature.
    /*!
     * If no argument is supplied, this
     * method returns the maximum temperature for which \e all
     * parameterizations are valid. If an integer index k is
     * supplied, then the value returned is the maximum
     * temperature for parameterization k.
     *
     * @param k  species index
     */
    virtual doublereal maxTemp(size_t k = npos) const
    {
        if (k == npos) {
            return m_thigh_min;
        } else {
            return m_thigh[k];
        }
    }

    //! The reference-state pressure for species k.
    /*!
     *
     * returns the reference state pressure in Pascals for
     * species k. If k is left out of the argument list,
     * it returns the reference state pressure for the first
     * species.
     * Note that some SpeciesThermo implementations, such
     * as those for ideal gases, require that all species
     * in the same phase have the same reference state pressures.
     *
     * @param k species index
     */
    virtual doublereal refPressure(size_t k = npos) const
    {
        return m_p0;
    }

    //! This utility function reports the type of parameterization
    //! used for the species with index number index.
    /*!
     *
     * @param index  Species index
     */
    virtual int reportType(size_t index) const
    {
        return SHOMATE;
    }

    /*!
     * This utility function reports back the type of
     * parameterization and all of the parameters for the
     * species, index.
     *
     * @param index     Species index
     * @param type      Integer type of the standard type
     * @param c         Vector of coefficients used to set the
     *                  parameters for the standard state.
     *
     * @param minTemp   output - Minimum temperature
     * @param maxTemp   output - Maximum temperature
     * @param refPressure output - reference pressure (Pa).
     */
    virtual void reportParams(size_t index, int& type, doublereal* const c, doublereal& minTemp, doublereal& maxTemp,
                              doublereal& refPressure) const
    {
        type = reportType(index);
        if (type == SHOMATE) {
            size_t grp = m_group_map[index];
            size_t pos = m_posInGroup_map[index];
            int itype = SHOMATE;
            const std::vector<ShomatePoly<ValAndDerivType> > &mlg = m_low[grp - 1];
            const std::vector<ShomatePoly<ValAndDerivType> > &mhg = m_high[grp - 1];
            const ShomatePoly<ValAndDerivType> * lowPoly = & (mlg[pos]);
            const ShomatePoly<ValAndDerivType> * highPoly = & (mhg[pos]);
            doublereal tmid = lowPoly->maxTemp();
            c[0] = tmid;
            size_t n;
            double ttemp;
            lowPoly->reportParameters(n, itype, minTemp, ttemp, refPressure, c + 1);
            if (n != index) {
                throw CanteraError("  ", "confused");
            }
            if (itype != SHOMATE && itype != SHOMATE1) {
                throw CanteraError("  ", "confused");
            }
            highPoly->reportParameters(n, itype, ttemp, maxTemp, refPressure, c + 8);
            if (n != index) {
                throw CanteraError("  ", "confused");
            }
            if (itype != SHOMATE && itype != SHOMATE1) {
                throw CanteraError("  ", "confused");
            }
        } else {
            throw CanteraError(" ", "confused");
        }
    }

#ifdef H298MODIFY_CAPABILITY

    virtual doublereal reportOneHf298(size_t k) const
    {
        doublereal h;
        doublereal t = 298.15;

        int grp = m_group_map[k];
        int pos = m_posInGroup_map[k];
        const std::vector<ShomatePoly<ValAndDerivType> > &mlg = m_low[grp - 1];
        const ShomatePoly<ValAndDerivType>* nlow = & (mlg[pos]);

        doublereal tmid = nlow->maxTemp();
        if (t <= tmid) {
            h = nlow->reportHf298();
        } else {
            const std::vector<ShomatePoly<ValAndDerivType> > &mhg = m_high[grp - 1];
            const ShomatePoly<ValAndDerivType>* nhigh = & (mhg[pos]);
            h = nhigh->reportHf298();
        }
        return h;
    }

    virtual void modifyOneHf298(const size_t k, const doublereal Hf298New)
    {

        int grp = m_group_map[k];
        int pos = m_posInGroup_map[k];
        std::vector<ShomatePoly<ValAndDerivType> > &mlg = m_low[grp - 1];
        ShomatePoly<ValAndDerivType>* nlow = & (mlg[pos]);
        std::vector<ShomatePoly<ValAndDerivType> > &mhg = m_high[grp - 1];
        ShomatePoly<ValAndDerivType>* nhigh = & (mhg[pos]);
        doublereal tmid = nlow->maxTemp();

        double hnow = reportOneHf298(k);
        double delH = Hf298New - hnow;
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
protected:

    //! Vector of vector of NasaPoly1's for the high temp region.
    /*!
     * This is the high temp region representation.
     * The first Length is equal to the number of groups.
     * The second vector is equal to the number of species
     * in that particular group.
     */
    std::vector<std::vector<ShomatePoly<ValAndDerivType> > > m_high;

    //! Vector of vector of NasaPoly1's for the low temp region.
    /*!
     * This is the low temp region representation.
     * The first Length is equal to the number of groups.
     * The second vector is equal to the number of species
     * in that particular group.
     */
    std::vector<std::vector<ShomatePoly<ValAndDerivType> > > m_low;

    //! Map between the midpoint temperature, as an int, to the group number
    /*!
     * Length is equal to the number of groups. Only used in the setup.
     */
    std::map<int, int> m_index;

    //! Vector of log temperature limits
    /*!
     * Length is equal to the number of groups.
     */
    vector_fp m_tmid;

    //! Maximum value of the low temperature limit
    doublereal m_tlow_max;

    //! Minimum value of the high temperature limit
    doublereal m_thigh_min;

    //! Vector of low temperature limits (species index)
    /*!
     * Length is equal to number of species
     */
    vector_fp m_tlow;

    //! Vector of low temperature limits (species index)
    /*!
     * Length is equal to number of species
     */
    vector_fp m_thigh;

    //! Reference pressure (Pa)
    /*!
     * all species must have the same reference pressure.
     */
    doublereal m_p0;

    //! number of groups
    int m_ngroups;

    //! Vector of temperature polynomials
    mutable std::vector<ValAndDerivType> m_t;

    /*!
     * This map takes as its index, the species index in the phase.
     * It returns the group index, where the temperature polynomials
     * for that species are stored. group indices start at 1,
     * so a decrement is always performed to access vectors.
     */
    mutable std::map<size_t, size_t> m_group_map;

    /*!
     * This map takes as its index, the species index in the phase.
     * It returns the position index within the group, where the
     * temperature polynomials for that species are stored.
     */
    mutable std::map<size_t, size_t> m_posInGroup_map;

    friend class ShomateThermo<doublereal> ;
    friend class ShomateThermo<doubleFAD> ;
};

}

#endif
