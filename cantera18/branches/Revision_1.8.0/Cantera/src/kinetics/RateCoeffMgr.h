/**
 *  @file RateCoeffMgr.h
 */

/*
 * $Author: dggoodwin $
 * $Revision: 1.1 $
 * $Date: 2007/05/04 14:27:23 $
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_RATECOEFF_MGR_H
#define CT_RATECOEFF_MGR_H

#include "utilities.h"
#include "RxnRates.h"

#include "ct_defs.h"
#include "ctexceptions.h"


namespace Cantera {

    /**
     * This rate coefficient manager supports one parameterization of
     * the rate constant of any type.
     */
    template<class R>
    class Rate1 {

    public:

        Rate1(){}
        virtual ~Rate1(){}

        /**
         * Install a rate coefficient calculator.
         * @param rxnNumber the reaction number
         * @param rateType  the rate type
         * @param m length of coefficient array
         * @param coefficients
         */
        int install( int rxnNumber,  int rateType, int m, 
            const doublereal* c ) {
            /*
	     * Check to see if the current reaction rate type
	     * is the same as the type of this class. If not,
	     * throw an error condition. 
	     */
            if (rateType != R::type()) 
                throw CanteraError("Rate1::install",
                    "incorrect rate coefficient type: "+int2str(rateType));

            // if any coefficient other than the first is non-zero, or
            // if alwaysComputeRate() is true, install a rate
            // calculator and return the index of the calculator.
            for (int i = 1; i < m; i++) {
                if (c[i] != 0.0 || R::alwaysComputeRate() ) {
                    m_rxn.push_back(rxnNumber);
                    m_rates.push_back(R(m, c));
                    return static_cast<int>(m_rates.size()) - 1;
                }
            }
            return -1;
        }
  
        /**
         * Return a reference to the nth rate coefficient calculator.
         * Note that this is not the same as the calculator for
         * reaction n, since reactions with constant rate coefficients
         * do not have a calculator.
         */  
        const R& rateCoeff(int loc) const { return m_rates[loc]; }

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
            TYPENAME_KEYWORD std::vector<R>::iterator b = m_rates.begin();
            TYPENAME_KEYWORD std::vector<R>::iterator e = m_rates.end();
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
            TYPENAME_KEYWORD std::vector<R>::const_iterator b = m_rates.begin();
            TYPENAME_KEYWORD std::vector<R>::const_iterator e = m_rates.end();
            doublereal recipT = 1.0/T;
            int i = 0;
            for (; b != e; ++b, ++i) {
	      // values[m_rxn[i]] = exp(b->update(logT, recipT));
	      values[m_rxn[i]] = b->updateRC(logT, recipT);
            }
        }

        void writeUpdate(std::ostream & output1, std::string key) {
	    output1 << key;
        }

    protected:
        std::vector<R>             m_rates;
        std::vector<int>           m_rxn;
        array_fp              m_const; // not used
    };



    /**
     * This rate coefficient manager supports two parameterizations of
     * any type.
     */
    template<class R1, class R2>
    class Rate2 {
    public:

        Rate2(){}
        virtual ~Rate2(){}

        int install( int rxnNumber,  int rateType, int m,
            const doublereal* c) {
            if (rateType == R1::type())
                return m_r1.install(rxnNumber, rateType, m, c);
            else if (rateType == R2::type())
                return m_r2.install(rxnNumber, rateType, m, c);
            else
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
