/**
 *  @file RateCoeffMgr.h
 */
// Copyright 2001  California Institute of Technology


#ifndef CT_RATECOEFF_MGR_H
#define CT_RATECOEFF_MGR_H

#include "cantera/base/utilities.h"
#include "RxnRates.h"

#include "cantera/base/ct_defs.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/stringUtils.h"

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
     * @param useAux flag to indicate that auxiliary rate information from
     *        rdata should be used.
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
        return m_rates.size() - 1;
    }

    /**
     * Return a reference to the nth rate coefficient calculator.
     * Note that this is not the same as the calculator for
     * reaction n, since reactions with constant rate coefficients
     * do not have a calculator.
     */
    const R& rateCoeff(int loc) const {
        return m_rates[loc];
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
        typename std::vector<R>::iterator b = m_rates.begin();
        typename std::vector<R>::iterator e = m_rates.end();
        int i = 0;
        for (; b != e; ++b, ++i) {
            b->update_C(c);
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
        typename std::vector<R>::const_iterator b = m_rates.begin();
        typename std::vector<R>::const_iterator e = m_rates.end();
        doublereal recipT = 1.0/T;
        int i = 0;
        for (; b != e; ++b, ++i) {
            // values[m_rxn[i]] = exp(b->update(logT, recipT));
            values[m_rxn[i]] = b->updateRC(logT, recipT);
        }
    }

    void writeUpdate(std::ostream& output1, std::string key) {
        output1 << key;
    }

    size_t nReactions() const {
        return m_rates.size();
    }

protected:
    std::vector<R>             m_rates;
    std::vector<size_t>           m_rxn;
    vector_fp              m_const; //!< @deprecated not used
};



/**
 * This rate coefficient manager supports two parameterizations of
 * any type.
 */
template<class R1, class R2>
class Rate2
{
public:

    Rate2() {}
    virtual ~Rate2() {}

    int install(size_t rxnNumber, int rateType, size_t m,
                const doublereal* c) {
        if (rateType == R1::type()) {
            return m_r1.install(rxnNumber, rateType, m, c);
        } else if (rateType == R2::type()) {
            return m_r2.install(rxnNumber, rateType, m, c);
        } else
            throw CanteraError("Rate2::install",
                               "unknown rate coefficient type");
        return -1;
    }

    void update(doublereal T, doublereal logT,
                doublereal* values) {
        m_r1.update(T, logT, values);
        m_r2.update(T, logT, values);
    }

protected:

    Rate1<R1> m_r1;
    Rate1<R2> m_r2;
};

}

#endif
