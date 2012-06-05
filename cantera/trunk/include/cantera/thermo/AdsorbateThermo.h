/**
 *  @file AdsorbateThermo.h
 *
 *  Header for a single-species standard
 *  state object derived from \link Cantera::SpeciesThermoInterpType
 *  SpeciesThermoInterpType\endlink based on the expressions for the
 *  thermo properties of a species with several vibrational models.
 *
 */
// Copyright 2007  California Institute of Technology


#ifndef CT_ADSORBATE_H
#define CT_ADSORBATE_H

#include "SpeciesThermoInterpType.h"

namespace Cantera
{

/**
 * An adsorbed surface species.
 *
 * This class is designed specifically for use by the class
 * GeneralSpeciesThermo. It implements a model for the
 * thermodynamic properties of a molecule that can be modeled as a
 * set of independent quantum harmonic oscillators.
 *
 * @ingroup spthermo
 */
class Adsorbate : public SpeciesThermoInterpType
{

public:

    //! Empty constructor
    Adsorbate()
        : m_lowT(0.0),
          m_highT(0.0),
          m_index(0),
          m_nFreqs(0) {
    }


    //! Full Constructor
    /*!
     * @param n         Species index
     * @param tlow      output - Minimum temperature
     * @param thigh     output - Maximum temperature
     * @param pref      output - reference pressure (Pa).
     */
    Adsorbate(size_t n, doublereal tlow, doublereal thigh, doublereal pref,
              const doublereal* coeffs) : m_lowT(tlow),
        m_highT(thigh),
        m_index(n) {
        m_nFreqs = int(coeffs[0]);
        m_be = coeffs[1];
        m_freq.resize(m_nFreqs);
        std::copy(coeffs+2, coeffs + 2 + m_nFreqs, m_freq.begin());
    }

    /// Copy Constructor
    Adsorbate(const Adsorbate& b) :
        m_lowT(b.m_lowT),
        m_highT(b.m_highT),
        m_Pref(b.m_Pref),
        m_index(b.m_index),
        m_be(b.m_be) {
        m_nFreqs = b.m_nFreqs;
        std::copy(b.m_freq.begin(), b.m_freq.begin() + m_nFreqs,
                  m_freq.begin());
    }

    //! destructor
    virtual ~Adsorbate() {}

    //! duplicator
    virtual SpeciesThermoInterpType*
    duplMyselfAsSpeciesThermoInterpType() const {
        Adsorbate* np = new Adsorbate(*this);
        return (SpeciesThermoInterpType*) np;
    }

    virtual void install(std::string name, size_t index, int type,
                         const doublereal* c, doublereal minTemp, doublereal maxTemp,
                         doublereal refPressure) {
        m_be = c[1];
        m_nFreqs = int(c[0]);
        for (size_t n = 0; n < m_nFreqs; n++) {
            m_freq[n] = c[n+2];
        }
        m_index = index;

        m_lowT  = minTemp;
        m_highT = maxTemp;
        m_Pref = refPressure;
    }


    //! Returns the minimum temperature that the thermo
    //! parameterization is valid
    virtual doublereal minTemp() const     {
        return m_lowT;
    }

    //! Returns the maximum temperature that the thermo
    //! parameterization is valid
    virtual doublereal maxTemp() const     {
        return m_highT;
    }

    //! Returns the reference pressure (Pa)
    virtual doublereal refPressure() const {
        return OneAtm;
    }

    //! Returns an integer representing the type of parameterization
    virtual int reportType() const {
        return ADSORBATE;
    }

    //! Returns an integer representing the species index
    virtual size_t speciesIndex() const {
        return m_index;
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
        h_RT[m_index] = _energy_RT(temp);
        cp_R[m_index] = (temp*h_RT[m_index]
                         - (temp-0.01)*_energy_RT(temp-0.01))/0.01;
        s_R[m_index] = h_RT[m_index] - _free_energy_RT(temp);
    }

    //! This utility function reports back the type of
    /*! parameterization and all of the parameters for the
     * species, index.
     *
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
        type = ADSORBATE;
        tlow = m_lowT;
        thigh = m_highT;
        pref = m_Pref;
        coeffs[0] = static_cast<double>(m_nFreqs);
        coeffs[1] = m_be;
        for (size_t i = 2; i < m_nFreqs+2; i++) {
            coeffs[i] = m_freq[i-2];
        }
    }

protected:
    //!  lowest valid temperature
    doublereal m_lowT;
    //! Highest valid temperatre
    doublereal m_highT;
    //! Reference state pressure
    doublereal m_Pref;
    //! species index
    size_t m_index;
    size_t m_nFreqs;
    //! array of vib frequencies
    vector_fp m_freq;
    doublereal m_be;


    doublereal _energy_RT(double T) const {
        doublereal x, hnu_kt, hnu, sum = 0.0;
        doublereal kt = T*Boltzmann;
        for (size_t i = 0; i < m_nFreqs; i++) {
            hnu = Planck * m_freq[i];
            hnu_kt = hnu/kt;
            x = exp(-hnu_kt);
            sum += hnu_kt * x/(1.0 - x);
        }
        return sum + m_be/(GasConstant*T);
    }

    doublereal _free_energy_RT(double T) const {
        doublereal x, hnu_kt, sum = 0.0;
        doublereal kt = T*Boltzmann;
        for (size_t i = 0; i < m_nFreqs; i++) {
            hnu_kt = Planck * m_freq[i] / kt;
            x = exp(-hnu_kt);
            sum += log(1.0 - x);
        }
        return sum + m_be/(GasConstant*T);
    }

    doublereal _entropy_R(double T) const {
        return _energy_RT(T) - _free_energy_RT(T);
    }

};

}
#endif




