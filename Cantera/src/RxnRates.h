/**
 *  @file RxnRates.h
 *
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_RXNRATES_H
#define CT_RXNRATES_H

#include "reaction_defs.h"

namespace Cantera {
    

    class Arrhenius {

    public:
        static int type(){ return ARRHENIUS; }        
        Arrhenius() : m_b (0.0), m_E (0.0) {}
        Arrhenius( const doublereal* c )
            : m_b (c[1]), m_E (c[2]) { m_logA = log(c[0]);}

        void update_C(const doublereal* c) {}
        
        doublereal update(doublereal logT, doublereal recipT) const {
            return m_logA + m_b*logT - m_E*recipT;
        }

        doublereal update_dT(doublereal logT, doublereal recipT) const {
            return recipT*(m_b + m_E*recipT);
        }
        
        void writeUpdateRHS(ostream& s) const {
            s << " exp(" << m_logA;
            if (m_b != 0.0) s << " + " << m_b << " * tlog"; 
            if (m_E != 0.0) s << " - " << m_E << " * rt";
            s << ");" << endl;
        }

        doublereal activationEnergy_R() const {
            return m_E;
        }

    protected:
        doublereal m_logA, m_b, m_E;
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


