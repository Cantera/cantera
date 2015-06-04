/**
 *  @file RateCoeffMgr.h
 */
// Copyright 2001  California Institute of Technology


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
     * @param rdata rate coefficient specification for the reaction
     */
    size_t install(size_t rxnNumber, const ReactionData& rdata) {
        /*
        * Check to see if the current reaction rate type is the same as the
        * type of this class. If not, throw an error condition.
        */
        if (rdata.rateCoeffType != R::type())
            throw CanteraError("Rate1::install",
                               "incorrect rate coefficient type: "+int2str(rdata.rateCoeffType) + ". Was Expecting type: "+ int2str(R::type()));

        // Install a rate calculator and return the index of the calculator.
        m_rxn.push_back(rxnNumber);
        m_rates.push_back(R(rdata));
        m_indices[rxnNumber] = m_rxn.size() - 1;
        return m_rates.size() - 1;
    }

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
     * Update the concentration-dependent parts of the rate
     * coefficient, if any. Used by class SurfaceArrhenius to
     * compute coverage-dependent * modifications to the Arrhenius
     * parameters. The array c should contain whatever data the
     * particular rate coefficient class needs to update its
     * rates.  Note that this method does not return anything. To
     * get the updated rates, method update must be called after
     * the call to update_C.
     */
    void update_C(const doublereal* c) {
        for (size_t i = 0; i != m_rates.size(); i++) {
            m_rates[i].update_C(c);
        }
    }

    /**
     * Write the rate coefficients into array values. Each
     * calculator writes one entry in values, at the location
     * specified by the reaction number when it was
     * installed. Note that nothing will be done for reactions
     * that have constant rates. The array values should be
     * preloaded with the constant rate coefficients.
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

protected:
    std::vector<R>             m_rates;
    std::vector<size_t>           m_rxn;

    //! map reaction number to index in m_rxn / m_rates
    std::map<size_t, size_t> m_indices;
};

}

#endif
