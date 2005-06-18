/**
 *  @file RxnRates.h
 *
 */

/* $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_RXNRATES_H
#define CT_RXNRATES_H

#include "reaction_defs.h"
#include "ctexceptions.h"

namespace Cantera {
    

    class Arrhenius {

    public:
        static int type(){ return ARRHENIUS; }        
        Arrhenius() : m_b (0.0), m_E (0.0) {}
        Arrhenius( int csize, const doublereal* c )
            : m_b (c[1]), m_E (c[2]) { m_logA = log(c[0]);}
        Arrhenius( doublereal A, doublereal b, doublereal E)
            : m_b (b), m_E (E) { m_logA = log(A);}

        void update_C(const doublereal* c) {}
        
        doublereal update(doublereal logT, doublereal recipT) const {
            return m_logA + m_b*logT - m_E*recipT;
        }

        /// no longer used
        //doublereal update_dT(doublereal logT, doublereal recipT) const {
        //    return recipT*(m_b + m_E*recipT);
        //}
        
        void writeUpdateRHS(ostream& s) const {
            s << " exp(" << m_logA;
            if (m_b != 0.0) s << " + " << m_b << " * tlog"; 
            if (m_E != 0.0) s << " - " << m_E << " * rt";
            s << ");" << endl;
        }

        doublereal activationEnergy_R() const {
            return m_E;
        }

        static bool alwaysComputeRate() { return false;}

    protected:
        doublereal m_logA, m_b, m_E;
    };


    class ArrheniusSum {

    public:
        static int type(){ return ARRHENIUS_SUM; }        
        ArrheniusSum() : m_nterms(0) {}
        ArrheniusSum( int csize, const doublereal* c ) {
            m_nterms = 0;
            addArrheniusTerm(c[0], c[1], c[2]);
        }
        void addArrheniusTerm(doublereal A, doublereal b, doublereal E) {
            m_terms.push_back(Arrhenius(A, b, E));
            m_nterms++;
        }
            
        void update_C(const doublereal* c) {}
        
        doublereal update(doublereal logT, doublereal recipT) const {
            int n;
            doublereal f, fexp = 0.0;
            for (n = 0; n < m_nterms; n++) {
                f = m_terms[n].update(logT, recipT);
                fexp += exp(f);
            }
            return log(fexp);
        }

        //        doublereal update_dT(doublereal logT, doublereal recipT) const {
        //    throw CanteraError("ArrheniusSum::update_dT","not implemented.");
        //}
        
        void writeUpdateRHS(ostream& s) const {
            ;
        }

        //doublereal activationEnergy_R() const {
        //    return m_E;
        //}

        static bool alwaysComputeRate() { return true;}

    protected:
        vector<Arrhenius> m_terms;
        int m_nterms;
    };



    /**
     * An Arrhenius rate with coverage-dependent terms.
     */
    class SurfaceArrhenius {

    public:
        static int type(){ return ARRHENIUS; }        
        SurfaceArrhenius() : m_b (0.0), m_E (0.0),
                             m_acov(0.0), m_ecov(0.0), 
                             m_mcov(0.0), m_ncov(0), m_nmcov(0)
            {}
        SurfaceArrhenius( int csize, const doublereal* c )
            : m_b (c[1]), m_E (c[2]), 
              m_acov(0.0), m_ecov(0.0), m_mcov(0.0), m_ncov(0), m_nmcov(0) 
            { m_logA = log(c[0]);
              if (csize >= 7) {
                  for (int n = 3; n < csize-3; n += 4) {
                      addCoverageDependence(int(c[n]), 
                          c[n+1], c[n+2], c[n+3]);
                  }
              }
            }

        void addCoverageDependence(int k, doublereal a, 
            doublereal m, doublereal e) {
            m_ncov++;
            m_sp.push_back(k);
            m_ac.push_back(a);
            m_ec.push_back(e);
            if (m != 0.0) {
                m_msp.push_back(k);
                m_mc.push_back(m);
                m_nmcov++;
            }
        }
            
        void update_C(const doublereal* theta) {
            m_acov = 0.0;
            m_ecov = 0.0;
            m_mcov = 0.0;
            int n, k;
            doublereal th;
            for (n = 0; n < m_ncov; n++) {
                k = m_sp[n];
                m_acov += m_ac[n] * theta[k];
                m_ecov += m_ec[n] * theta[k];
            }
            for (n = 0; n < m_nmcov; n++) {
                k = m_msp[n];
                // changed n to k, dgg 1/22/04
                th = fmaxx(theta[k], Tiny);
                //                th = fmaxx(theta[n], Tiny);
                m_mcov += m_mc[n]*log(th);
            }
        }
        
        doublereal update(doublereal logT, doublereal recipT) const {
            return m_logA + m_acov + m_b*logT 
                - (m_E + m_ecov)*recipT + m_mcov;
        }

        doublereal activationEnergy_R() const {
            return m_E + m_ecov;
        }

        static bool alwaysComputeRate() { return true;}

    protected:
        doublereal m_logA, m_b, m_E;
        doublereal m_acov, m_ecov, m_mcov;
        vector_int m_sp, m_msp;
        vector_fp m_ac, m_ec, m_mc;
        int m_ncov, m_nmcov; 
    };


#ifdef INCL_TST

    class TST {

    public:
        static int type(){ return TSTRATE; }        
        TST() {}
        TST( const vector_fp& c ) {
            m_b.resize(10);
            copy(c.begin(), c.begin() + 10, m_b.begin());
            m_k = int(c[10]);
        }
        
        void update_C(const vector_fp& c) {
            doublereal ck = c[m_k];
            delta_s0 = m_b[0] + m_b[1]*ck + m_b[2]*ck*ck;
            delta_e0 = m_b[5] + m_b[6]*ck + m_b[7]*ck*ck;
        }

        doublereal update(doublereal logT, doublereal recipT) const {
            doublereal delta_s = delta_s0*(1.0 + m_b[3]*logT + m_b[4]*recipT);
            doublereal delta_E = delta_e0*(1.0 + m_b[8]*logT + m_b[9]*recipT);
            return logBoltz_Planck + logT + delta_s -  delta_E*recipT;
        }
        
        void writeUpdateRHS(ostream& s) const {}

    protected:
        doublereal delta_s0, delta_e0;
        int m_k;
        vector_fp m_b;
    };

#endif
    
}


//     class LandauTeller {

//     public:
//         static int type(){ return LANDAUTELLER; }        
//         LandauTeller(){}
//         LandauTeller( const vector_fp& c ) : m_c(c) { m_c[0] = log(c[0]); }
        
//         doublereal update(doublereal logT, doublereal recipT) const {
//             return m_c[0] + m_c[1]*tt[1] - m_c[2]*tt[2] 
//                 + m_c[3]*tt[3] + m_c[4]*tt[4];
//         }
        
//         //void writeUpdateRHS(ostream& s) const {
//         //     s << exp(m_logA); 
//         //    s << " * exp(";
//         //    if (m_b != 0.0) s << m_b << " * tlog"; 
//         //    if (m_E != 0.0) s << " - " << m_E << " * rt";
//         //    if (m_E != 0.0) s << " - " << m_E << " * rt";
//         //        s << ");" << endl;
//         //    }
//         //}

//     protected:
//         doublereal m_logA, m_b, m_E;
//     };
    
//}


#endif


