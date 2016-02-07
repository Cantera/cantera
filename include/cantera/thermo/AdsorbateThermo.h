/**
 * @file AdsorbateThermo.h
 *
 * Header for a single-species standard state object derived from \link
 * Cantera::SpeciesThermoInterpType SpeciesThermoInterpType\endlink based on the
 * expressions for the thermo properties of a species with several vibrational
 * models.
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
 * GeneralSpeciesThermo. It implements a model for the thermodynamic properties
 * of a molecule that can be modeled as a set of independent quantum harmonic
 * oscillators.
 *
 * @ingroup spthermo
 */
class Adsorbate : public SpeciesThermoInterpType
{
public:
    //! Empty constructor
    Adsorbate() {}

    //! Full Constructor
    /*!
     * @param tlow      output - Minimum temperature
     * @param thigh     output - Maximum temperature
     * @param pref      output - reference pressure (Pa).
     */
    Adsorbate(double tlow, double thigh, double pref, const double* coeffs)
        : SpeciesThermoInterpType(tlow, thigh, pref)
    {
        m_freq.resize(int(coeffs[0]));
        m_be = coeffs[1];
        std::copy(coeffs+2, coeffs + 2 + m_freq.size(), m_freq.begin());
    }

    virtual SpeciesThermoInterpType*
    duplMyselfAsSpeciesThermoInterpType() const {
        return new Adsorbate(*this);
    }

    virtual int reportType() const {
        return ADSORBATE;
    }

    void updatePropertiesTemp(const doublereal temp,
                              doublereal* cp_R,
                              doublereal* h_RT,
                              doublereal* s_R) const {
        *h_RT = _energy_RT(temp);
        *cp_R = (temp**h_RT - (temp-0.01)*_energy_RT(temp-0.01))/0.01;
        *s_R = *h_RT - _free_energy_RT(temp);
    }

    void reportParameters(size_t& n, int& type,
                          doublereal& tlow, doublereal& thigh,
                          doublereal& pref,
                          doublereal* const coeffs) const {
        n = 0;
        type = ADSORBATE;
        tlow = m_lowT;
        thigh = m_highT;
        pref = m_Pref;
        coeffs[0] = static_cast<double>(m_freq.size());
        coeffs[1] = m_be;
        for (size_t i = 2; i < m_freq.size()+2; i++) {
            coeffs[i] = m_freq[i-2];
        }
    }

protected:
    //! array of vib frequencies
    vector_fp m_freq;
    doublereal m_be;

    doublereal _energy_RT(double T) const {
        doublereal x, hnu_kt, hnu, sum = 0.0;
        doublereal kt = T*Boltzmann;
        for (size_t i = 0; i < m_freq.size(); i++) {
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
        for (size_t i = 0; i < m_freq.size(); i++) {
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
