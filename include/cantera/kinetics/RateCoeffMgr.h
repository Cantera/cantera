/**
 *  @file RateCoeffMgr.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_RATECOEFF_MGR_H
#define CT_RATECOEFF_MGR_H

#include "RxnRates.h"

namespace Cantera
{

/**
 * This rate coefficient manager supports one parameterization of
 * the rate constant of any type.
 */
template<class R>
class Rate1
{
public:
    Rate1() {}
    virtual ~Rate1() {}

    /**
     * Install a rate coefficient calculator.
     * @param rxnNumber the reaction number
     * @param rate rate coefficient specification for the reaction
     */
    void install(size_t rxnNumber, const R& rate) {
        m_rxn.push_back(rxnNumber);
        m_rates.push_back(rate);
        m_indices[rxnNumber] = m_rxn.size() - 1;
    }

    //! Replace an existing rate coefficient calculator
    void replace(size_t rxnNumber, const R& rate) {
        size_t i = m_indices[rxnNumber];
        m_rates[i] = rate;
    }

    /**
     * Update the concentration-dependent parts of the rate coefficient, if any.
     * Used by class SurfaceArrhenius to compute coverage-dependent *
     * modifications to the Arrhenius parameters. The array c should contain
     * whatever data the particular rate coefficient class needs to update its
     * rates.  Note that this method does not return anything. To get the
     * updated rates, method update must be called after the call to update_C.
     */
    void update_C(const doublereal* c) {
        for (size_t i = 0; i != m_rates.size(); i++) {
            m_rates[i].update_C(c);
        }
    }

    /**
     * Write the rate coefficients into array values. Each calculator writes one
     * entry in values, at the location specified by the reaction number when it
     * was installed. Note that nothing will be done for reactions that have
     * constant rates. The array values should be preloaded with the constant
     * rate coefficients.
     */
    void update(doublereal T, doublereal logT, doublereal* values) {
        doublereal recipT = 1.0/T;
        for (size_t i = 0; i != m_rates.size(); i++) {
            values[m_rxn[i]] = m_rates[i].updateRC(logT, recipT);
        }
    }

    size_t nReactions() const {
        return m_rates.size();
    }

    //! Return effective preexponent for the specified reaction.
    /*!
     *  Returns effective preexponent, accounting for surface coverage
     *  dependencies. Used in InterfaceKinetics.
     *
     *  @param irxn Reaction number in the kinetics mechanism
     *  @return Effective preexponent
     */
    double effectivePreExponentialFactor(size_t irxn) {
        return m_rates[irxn].preExponentialFactor();
    }

    //! Return effective activation energy for the specified reaction.
    /*!
     *  Returns effective activation energy, accounting for surface coverage
     *  dependencies. Used in InterfaceKinetics.
     *
     *  @param irxn Reaction number in the kinetics mechanism
     *  @return Effective activation energy divided by the gas constant
     */
    double effectiveActivationEnergy_R(size_t irxn) {
        return m_rates[irxn].activationEnergy_R();
    }

    //! Return effective temperature exponent for the specified  reaction.
    /*!
     *  Returns effective temperature exponent, accounting for surface coverage
     *  dependencies. Used in InterfaceKinetics. Current parameterization in
     *  SurfaceArrhenius does not change this parameter with the change in
     *  surface coverages.
     *
     *  @param irxn Reaction number in the kinetics mechanism
     *  @return Effective temperature exponent
     */
    double effectiveTemperatureExponent(size_t irxn) {
        return m_rates[irxn].temperatureExponent();
    }

protected:
    std::vector<R> m_rates;
    std::vector<size_t> m_rxn;

    //! map reaction number to index in m_rxn / m_rates
    std::map<size_t, size_t> m_indices;
};

}

#endif
