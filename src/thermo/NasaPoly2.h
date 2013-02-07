/**
 *  @file NasaPoly2.h
 *  Header for a single-species standard state object derived
 *  from \link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType\endlink  based
 *  on the NASA temperature polynomial form applied to two temperature regions
 *  (see \ref spthermo and class \link Cantera::NasaPoly2 NasaPoly2\endlink).
 *
 *  Two zoned Nasa polynomial parameterization
 */
// Copyright 2001  California Institute of Technology


#ifndef CT_NASAPOLY2_H
#define CT_NASAPOLY2_H

#include "cantera/thermo/SpeciesThermoInterpType.h"

namespace Cantera
{

/**
 *
 *
 * The NASA polynomial parameterization for two temperature ranges.
 * This parameterization expresses the heat capacity as a
 * fourth-order polynomial. Note that this is the form used in the
 * 1971 NASA equilibrium program and by the Chemkin software
 * package, but differs from the form used in the more recent NASA
 * equilibrium program.
 *
 * Seven coefficients \f$(a_0,\dots,a_6)\f$ are used to represent
 * \f$ c_p^0(T)\f$, \f$ h^0(T)\f$, and \f$ s^0(T) \f$ as
 * polynomials in \f$ T \f$ :
 * \f[
 * \frac{c_p(T)}{R} = a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4
 * \f]
 * \f[
 * \frac{h^0(T)}{RT} = a_0 + \frac{a_1}{2} T + \frac{a_2}{3} T^2
 * + \frac{a_3}{4} T^3 + \frac{a_4}{5} T^4  + \frac{a_5}{T}.
 * \f]
 * \f[
 * \frac{s^0(T)}{R} = a_0\ln T + a_1 T + \frac{a_2}{2} T^2
 + \frac{a_3}{3} T^3 + \frac{a_4}{4} T^4  + a_6.
 * \f]
 *
 * This class is designed specifically for use by the class
 * GeneralSpeciesThermo.
 *
 * @ingroup spthermo
 */
class NasaPoly2 : public SpeciesThermoInterpType
{

public:

    //! Empty constructor
    NasaPoly2()
        : m_lowT(0.0),
          m_midT(0.0),
          m_highT(0.0),
          m_Pref(0.0),
          m_index(0),
          m_coeff(15, 0.0) {
    }

    //! Full Constructor
    /*!
     * @param n         Species index
     * @param tlow      output - Minimum temperature
     * @param thigh     output - Maximum temperature
     * @param pref      output - reference pressure (Pa).
     * @param coeffs    Vector of coefficients used to set the
     *                  parameters for the standard state.
     */
    NasaPoly2(size_t n, doublereal tlow, doublereal thigh, doublereal pref,
              const doublereal* coeffs) :
        m_lowT(tlow),
        m_highT(thigh),
        m_Pref(pref),
        mnp_low(n, tlow, coeffs[0], pref, coeffs +1),
        mnp_high(n, tlow, thigh, pref, coeffs + 8),
        m_index(n),
        m_coeff(15, 0.0) {

        std::copy(coeffs, coeffs + 15, m_coeff.begin());
        m_midT = coeffs[0];
    }

    //! Copy Constructor
    /*!
     * @param b object to be copied.
     */
    NasaPoly2(const NasaPoly2& b) :
        m_lowT(b.m_lowT),
        m_midT(b.m_midT),
        m_highT(b.m_highT),
        m_Pref(b.m_Pref),
        mnp_low(b.mnp_low),
        mnp_high(b.mnp_high),
        m_index(b.m_index),
        m_coeff(b.m_coeff) {
    }

    //! Assignment operator
    /*!
     * @param b object to be copied.
     */
    NasaPoly2& operator=(const NasaPoly2& b) {
        if (&b != this) {
            m_lowT   = b.m_lowT;
            m_midT   = b.m_midT;
            m_highT  = b.m_highT;
            m_Pref   = b.m_Pref;
            m_index  = b.m_index;
            m_coeff  = b.m_coeff;
            mnp_low  = b.mnp_low;
            mnp_high = b.mnp_high;
        }
        return *this;
    }

    //! destructor
    virtual ~NasaPoly2() { }

    //! duplicator
    virtual SpeciesThermoInterpType*
    duplMyselfAsSpeciesThermoInterpType() const {
        NasaPoly2* np = new NasaPoly2(*this);
        return (SpeciesThermoInterpType*) np;
    }

    //! Returns the minimum temperature that the thermo
    //! parameterization is valid
    doublereal minTemp() const     {
        return m_lowT;
    }

    //! Returns the maximum temperature that the thermo
    //! parameterization is valid
    doublereal maxTemp() const     {
        return m_highT;
    }

    //! Returns the reference pressure (Pa)
    doublereal refPressure() const {
        return m_Pref;
    }

    //! Returns an integer representing the type of parameterization
    virtual int reportType() const {
        return NASA2;
    }

    //! Returns an integer representing the species index
    virtual size_t speciesIndex() const {
        return m_index;
    }


    //! Update the properties for this species, given a temperature polynomial
    /*!
     * This method is called with a pointer to an array containing the functions of
     * temperature needed by this  parameterization, and three pointers to arrays where the
     * computed property values should be written. This method updates only one value in
     * each array.
     *
     * Temperature Polynomial:
     *  tt[0] = t;
     *  tt[1] = t*t;
     *  tt[2] = m_t[1]*t;
     *  tt[3] = m_t[2]*t;
     *  tt[4] = 1.0/t;
     *  tt[5] = std::log(t);
     *
     * @param tt  vector of temperature polynomials
     * @param cp_R    Vector of Dimensionless heat capacities.
     *                (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies.
     *                (length m_kk).
     * @param s_R     Vector of Dimensionless entropies.
     *                (length m_kk).
     */
    void updateProperties(const doublereal* tt,
                          doublereal* cp_R, doublereal* h_RT, doublereal* s_R) const {

        double T = tt[0];
        if (T <= m_midT) {
            mnp_low.updateProperties(tt, cp_R, h_RT, s_R);
        } else {
            mnp_high.updateProperties(tt, cp_R, h_RT, s_R);
        }
    }

    //! Compute the reference-state property of one species
    /*!
     * Given temperature T in K, this method updates the values of
     * the non-dimensional heat capacity at constant pressure,
     * enthalpy, and entropy, at the reference pressure, Pref
     * of one of the species. The species index is used
     * to reference into the cp_R, h_RT, and s_R arrays.
     *
     * @param temp    Temperature (Kelvin)
     * @param cp_R    Vector of Dimensionless heat capacities.
     *                (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies.
     *                (length m_kk).
     * @param s_R     Vector of Dimensionless entropies.
     *                (length m_kk).
     */
    void updatePropertiesTemp(const doublereal temp,
                              doublereal* cp_R,
                              doublereal* h_RT,
                              doublereal* s_R) const {
        if (temp <= m_midT) {
            mnp_low.updatePropertiesTemp(temp, cp_R, h_RT, s_R);
        } else {
            mnp_high.updatePropertiesTemp(temp, cp_R, h_RT, s_R);
        }
    }

    //!This utility function reports back the type of
    //! parameterization and all of the parameters for the
    //! species, index.
    /*!
     * All parameters are output variables
     *
     * @param n         Species index
     * @param type      Integer type of the standard type
     * @param tlow      output - Minimum temperature
     * @param thigh     output - Maximum temperature
     * @param pref      output - reference pressure (Pa).
     * @param coeffs    Vector of coefficients used to set the
     *                  parameters for the standard state.
     */
    void reportParameters(size_t& n, int& type,
                          doublereal& tlow, doublereal& thigh,
                          doublereal& pref,
                          doublereal* const coeffs) const {
        n = m_index;
        type = NASA2;
        tlow = m_lowT;
        thigh = m_highT;
        pref = m_Pref;
        for (int i = 0; i < 15; i++) {
            coeffs[i] = m_coeff[i];
        }
    }

#ifdef H298MODIFY_CAPABILITY

    doublereal reportHf298(doublereal* const h298 = 0) const {
        double h;
        if (298.15 <= m_midT) {
            h = mnp_low.reportHf298(0);
        } else {
            h = mnp_high.reportHf298(0);
        }
        if (h298) {
            h298[m_index] = h;
        }
        return h;
    }

    void modifyOneHf298(const size_t& k, const doublereal Hf298New) {
        if (k != m_index) {
            return;
        }

        doublereal h298now = reportHf298(0);
        doublereal delH = Hf298New - h298now;
        double h = mnp_low.reportHf298(0);
        double hnew = h + delH;
        mnp_low.modifyOneHf298(k, hnew);
        h  = mnp_high.reportHf298(0);
        hnew = h + delH;
        mnp_high.modifyOneHf298(k, hnew);
    }

#endif

protected:
    //!  lowest valid temperature
    doublereal m_lowT;
    //! Midrange temperature
    doublereal m_midT;
    //! Highest valid temperature
    doublereal m_highT;
    //! Reference state pressure
    doublereal m_Pref;
    //! NasaPoly1 object for the low temperature region.
    NasaPoly1 mnp_low;
    //! NasaPoly1 object for the high temperature region.
    NasaPoly1 mnp_high;
    //! species index
    size_t m_index;
    //! array of polynomial coefficients
    vector_fp m_coeff;

};

}
#endif




