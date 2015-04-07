/**
 *  @file AdsorbateThermo.h
 *
 *  Header for a single-species standard
 *  state object derived from \link Cantera::SpeciesThermoInterpType
 *  SpeciesThermoInterpType\endlink based on the expressions for the
 *  thermo properties of a species with several vibrational models.
 */
// Copyright 2007  California Institute of Technology

#ifndef CT_ADSORBATE_H
#define CT_ADSORBATE_H

#include "SpeciesThermoInterpType.h"
#include "cantera/base/global.h"

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
    Adsorbate() :
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
              const doublereal* coeffs)
        : SpeciesThermoInterpType(n, tlow, thigh, pref)
    {
        m_nFreqs = int(coeffs[0]);
        m_be = coeffs[1];
        m_freq.resize(m_nFreqs);
        std::copy(coeffs+2, coeffs + 2 + m_nFreqs, m_freq.begin());
    }

    /// Copy Constructor
    Adsorbate(const Adsorbate& b) :
        m_be(b.m_be) {
        m_nFreqs = b.m_nFreqs;
        std::copy(b.m_freq.begin(), b.m_freq.begin() + m_nFreqs,
                  m_freq.begin());
    }

    virtual SpeciesThermoInterpType*
    duplMyselfAsSpeciesThermoInterpType() const {
        Adsorbate* np = new Adsorbate(*this);
        return (SpeciesThermoInterpType*) np;
    }

    //! @deprecated Not a member of the base class
    virtual void install(const std::string& name, size_t index, int type,
                         const doublereal* c, doublereal minTemp_, doublereal maxTemp_,
                         doublereal refPressure_) {
        warn_deprecated("Adsorbate::install", "Not a member of the base class.");
        m_be = c[1];
        m_nFreqs = int(c[0]);
        for (size_t n = 0; n < m_nFreqs; n++) {
            m_freq[n] = c[n+2];
        }
        m_index = index;

        m_lowT  = minTemp_;
        m_highT = maxTemp_;
        m_Pref = refPressure_;
    }

    virtual int reportType() const {
        return ADSORBATE;
    }

    void updatePropertiesTemp(const doublereal temp,
                              doublereal* cp_R,
                              doublereal* h_RT,
                              doublereal* s_R) const {
        h_RT[m_index] = _energy_RT(temp);
        cp_R[m_index] = (temp*h_RT[m_index]
                         - (temp-0.01)*_energy_RT(temp-0.01))/0.01;
        s_R[m_index] = h_RT[m_index] - _free_energy_RT(temp);
    }

    //! @deprecated
    void reportParameters(size_t& n, int& type,
                          doublereal& tlow, doublereal& thigh,
                          doublereal& pref,
                          doublereal* const coeffs) const {
        warn_deprecated("AdsorbateThermo::reportParameters");
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
