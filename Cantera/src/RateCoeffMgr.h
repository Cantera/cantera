/**
 *  @file RateCoeffMgr.h
 *
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_RATECOEFF_MGR_H
#define CT_RATECOEFF_MGR_H

#include "utilities.h"
#include "RxnRates.h"

#ifdef HAVE_INTEL_MKL
#include "mkl_vml.h"
#endif

#include "ct_defs.h"
#include "ctexceptions.h"

namespace Cantera {

    // exception class
    class UnknownRateCoefficient : public CanteraError {
    public:
        UnknownRateCoefficient() {
            m_msg += "Unknown rate coefficient type.";
        }
    };

    /**
     * Virtual base class for rate coefficient managers. 
     */
//     class RateCoeffMgr {
//     public:
//         virtual int install( int rxnNumber,  int rateType, 
//             const vector_fp& c )=0;
//         virtual void update(doublereal T, doublereal logT, vector_fp& values)=0;
//         virtual void writeUpdate(ostream& s, string name) {}
//     };



    /**
     * This rate coefficient manager supports one parameterization of
     * any type.
     */
    template<class R>
    class Rate1 {

    public:

        Rate1(){}
        virtual ~Rate1(){}

        int install( int rxnNumber,  int rateType, int m, 
            const doublereal* c ) {
            if (rateType != R::type()) 
                throw UnknownRateCoefficient();
        
            //int m = c.size();

            // if any coefficient other than the first is non-zero,
            // install a rate calculator and return the index of the 
            // calculator.

            for (int i = 1; i < m; i++) {
                if (c[i] != 0.0) {
                    m_rxn.push_back(rxnNumber);
                    m_rates.push_back(R(c));
                    return m_rates.size() - 1;
                }
            }
            return -1;
        }
    
        const R& rateCoeff(int loc) const {
            return m_rates[loc];
        }

        void update_C(const doublereal* c) {
            TYPENAME_KEYWORD vector<R>::iterator b = m_rates.begin();
            TYPENAME_KEYWORD vector<R>::iterator e = m_rates.end();
            int i = 0;
            for (; b != e; ++b, ++i) {
                b->update_C(c);
            }
        }

        void update(doublereal T, doublereal logT, doublereal* values) {
            TYPENAME_KEYWORD vector<R>::const_iterator b = m_rates.begin();
            TYPENAME_KEYWORD vector<R>::const_iterator e = m_rates.end();
            doublereal recipT = 1.0/T;
            int i = 0;
            for (; b != e; ++b, ++i) {
                values[m_rxn[i]] = exp(b->update(logT, recipT));
            }
        }

        void update_dT(doublereal T, doublereal logT, doublereal dT,  
            doublereal* values) {
            TYPENAME_KEYWORD vector<R>::const_iterator b = m_rates.begin();
            TYPENAME_KEYWORD vector<R>::const_iterator e = m_rates.end();
            doublereal recipT = 1.0/T;
            int i = 0;
            for (; b != e; ++b, ++i) {
                values[m_rxn[i]] *= (1.0 + dT*b->update_dT(logT, recipT));
            }
        }

        virtual void writeUpdate(ostream& s, string name) {
            int nrates = m_rates.size();
            for (int i = 0; i < nrates; i++) {
                s << "    " << name << "[" << m_rxn[i] << "] = ";
                m_rates[i].writeUpdateRHS(s);
            }
        }

    protected:
        vector<R>             m_rates;
        vector<int>           m_rxn;
        array_fp              m_const;
    };



    /**
     * This rate coefficient manager supports two parameterizations of
     * any type.
     */
//     template<class R1, class R2>
//     class Rate2 : public RateCoeffMgr {
//     public:

//         Rate2(){}
//         virtual ~Rate2(){}

//         virtual int install( int rxnNumber,  int rateType, const vector_fp& c ) {
//             if (rateType == R1::type())
//                 return m_r1.install(rxnNumber, rateType, c);
//             else if (rateType == R2::type())
//                 return m_r2.install(rxnNumber, rateType, c);
//             else
//                 throw UnknownRateCoefficient();
//             return -1;
//         }
    
//         virtual void update(doublereal T, doublereal logT, vector_fp& values) {
//             m_r1.update(T, logT, values);
//             m_r2.update(T, logT, values);
//         }
//     protected:
//         Rate1<R1> m_r1;
//         Rate1<R2> m_r2;
//     };
    
}

#endif
