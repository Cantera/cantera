/**
 *
 * @file GasKineticsWriter.h
 *
 * $Author: hkmoffa $
 * $Revision: 1.3 $
 * $Date: 2008/12/17 17:09:37 $
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_GASKINETICSWRITER_H
#define CT_GASKINETICSWRITER_H

#define WRITE_UPDATE

#include <fstream>
#include <map>

#include "mix_defs.h"
#include "Kinetics.h"

#include "utilities.h"
#include "StoichManager.h"
#include "ThirdBodyMgr.h"
#include "FalloffMgr.h"
#include "RateCoeffMgr.h"
#include "Phase.h"

#include <cmath>
#include <cstdlib>

namespace Cantera {

    // forward references
    class Enhanced3BConc;
    class ReactionData;

    //! Class to write a hard-coded version of a mechanism.
    /*!
     * @ingroup kineticsmgr
     */
    class GasKineticsWriter {

    public:
        
        /// Default constructor.
        GasKineticsWriter();

        /// Destructor.
        virtual ~GasKineticsWriter(){}

        void init(int nsp);
        doublereal reactantStoichCoeff(int k, int i) const {
            return m_rrxn[k][i];
        }

        doublereal productStoichCoeff(int k, int i) const {
            return m_prxn[k][i];
        }
        void writeUpdateROP(){}

		void writeGetNetProductionRates(std::ostream& s, int nsp, int nrxns) {
            int i, k;
			s << "void get_wdot(const double* rop, double* wdot) {" << std::endl;
            for (k = 0; k < nsp; k++) {
                s << "  wdot[" << k << "] = ";
                doublereal net;
                bool empty = true;
                for (i = 0; i < nrxns; i++) {
                    net = productStoichCoeff(k,i) - reactantStoichCoeff(k,i);
                    if (net > 0.0) {
                        empty = false;
                        if (net == 1.0) 
                            s << " + rop[" << i << "]";
                        else 
                            s << " + " << net << "*rop[" << i << "]";
                    }
                    else if (net < 0.0) {
                        empty = false;
                        if (net == -1.0) 
                            s << " - rop[" << i << "]";
                        else 
                            s << " - " << -net << "*rop[" << i << "]";
                    }
                }
                if (empty) s << "0.0";
                s << ";" << endl;
            }
            s << "}" << endl;
        }

 
		void writeUpdateKc(std::ostream& s, int nsp, int nrxns) {
            int i, k, n, nn, ir;
            s << "void update_kc(const double* a, "
                "double exp_c0, double* rkc) {" << endl;
            for (i = 0; i != m_nrev; i++) {
                //if (isReversible(i)) {
                ir = m_revindex[i];
                s << "  rkc[" << ir << "] = ";
                bool empty = true;
                for (k = 0; k < nsp; k++) {
                    n = int(productStoichCoeff(k,ir));
                    for (nn = 0; nn != n; nn++) {
                        if (!empty) s << "*";
                        s << "a[" << k << "]";
                        empty = false;
                    }
                }
                if (m_dn[i] < 0.0) {
                    n = -m_dn[i];
                    for (nn = 0; nn < n; nn++) s << "*exp_c0";
                }
                s << "/(";
                empty = true;
                for (k = 0; k < nsp; k++) {
                    n = int(reactantStoichCoeff(k,ir));
                    for (nn = 0; nn < n; nn++) {
                        if (!empty) s << "*";
                        s << "a[" << k << "]";
                        empty = false;
                    }
                }
                if (m_dn[i] > 0.0) {
                    n = m_dn[i];
                    for (nn = 0; nn != n; nn++) s << "*exp_c0";
                }
                s << ");" << endl;
            }
            s << "}" << endl;
        }

		void writeEvalRopnet(std::ostream& s) {
            int i;
            s << "void eval_ropnet(const double* c, "
                "const double* rf, const double* rkc, double* r) {" << endl;
            for (i = 0; i < m_ii; i++) {
                s << "  r[" << i << "] = rf[" << i << "] * (" 
                  << m_reactantWriter.mult(i);
                if (isReversible(i)) {
                    s << " - rkc[" << i << "] * " 
                      << m_revProductWriter.mult(i);
                }
                s << ");" << endl;
            } 
            s << "}" << endl;
        }



		void writeUpdateRates(std::ostream& s) {
            s << "void update_rates(double t, double tlog, double* rf) {" << endl;
            s << "  double rt = 1.0/t;" << endl;
            m_rates.writeUpdate(s, "rf");
            s << "}" << endl;
        }

        ///  Add a reaction to the mechanism. 
        void addReaction(const ReactionData& r);
        
    protected:

        int m_kk, m_ii, m_nfall, m_nrev, m_nirrev;

        vector_int m_fallindx;
        
        Rate1<Arrhenius>                    m_falloff_low_rates;
        Rate1<Arrhenius>                    m_falloff_high_rates;        
        Rate1<Arrhenius>                    m_rates;        
        
		std::vector<int> m_irrev;
#ifdef INCL_STOICH_WRITER
        StoichWriter                       m_reactantWriter;
        StoichWriter                       m_revProductWriter;
        StoichWriter                       m_irrevProductWriter;
#endif
		mutable std::vector<std::map<int, doublereal> >     m_rrxn;
		mutable std::vector<std::map<int, doublereal> >     m_prxn;

        vector_int m_dn;
        vector_int m_revindex;

    private:

        int reactionNumber(){ return m_ii;}
        void addElementaryReaction(const ReactionData& r);
        void addThreeBodyReaction(const ReactionData& r);
        void addFalloffReaction(const ReactionData& r);
        
        void installReagents(const vector_int& r,
            const vector_int& p, bool reversible);

        virtual bool isReversible(int i) {
            if (find(m_revindex.begin(), m_revindex.end(), i) 
                < m_revindex.end()) return true;
            else return false;
        }
        bool m_finalized;
    };
}

#endif
