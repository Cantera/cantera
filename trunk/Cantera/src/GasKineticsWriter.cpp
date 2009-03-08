/**
 *  @file GasKineticsWriter.cpp 
 *
 */

// Copyright 2001  California Institute of Technology


// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ReactionData.h"
#include "GasKineticsWriter.h"

#include "StoichManager.h"
#include "Enhanced3BConc.h"
#include "ThirdBodyMgr.h"
#include "RateCoeffMgr.h"

//#include "ThermoPhase.h"

#include <iostream>
using namespace std;

namespace Cantera {    

    /**
     * Construct an empty reaction mechanism.
     */    
    GasKineticsWriter::
    GasKineticsWriter() : m_kk(0), m_ii(0), m_nfall(0), m_nrev(0), m_nirrev(0), 
                    m_finalized(false) {}

     void GasKineticsWriter::
     addReaction(const ReactionData& r) {

         if (r.reactionType == ELEMENTARY_RXN)      addElementaryReaction(r);
         else if (r.reactionType == THREE_BODY_RXN) addThreeBodyReaction(r);
         else if (r.reactionType == FALLOFF_RXN)    addFalloffReaction(r);

         // operations common to all reaction types
         installReagents( r.reactants, r.products, r.reversible );
         m_ii++;
     }


    void GasKineticsWriter::
    addFalloffReaction(const ReactionData& r) {

        // install high and low rate coeff calculators
        m_falloff_high_rates.install( m_nfall,
            r.rateCoeffType, r.rateCoeffParameters.size(),
            r.rateCoeffParameters.begin() );         
        m_falloff_low_rates.install( m_nfall, 
            r.rateCoeffType, r.auxRateCoeffParameters.size(), 
            r.auxRateCoeffParameters.begin() );
                
        // add this reaction number to the list of 
        // falloff reactions
        m_fallindx.push_back( reactionNumber() );

        // increment the falloff reaction counter
        ++m_nfall;
    }


    void GasKineticsWriter::
    addElementaryReaction(const ReactionData& r) {

        int iloc;
        // install rate coeff calculator
        iloc = m_rates.install( reactionNumber(),
            r.rateCoeffType, r.rateCoeffParameters.size(),
            r.rateCoeffParameters.begin() );
    }


    void GasKineticsWriter::
    addThreeBodyReaction(const ReactionData& r) {
            
        int iloc;
        // install rate coeff calculator
        iloc = m_rates.install( reactionNumber(),
            r.rateCoeffType, r.rateCoeffParameters.size(),
            r.rateCoeffParameters.begin() );
    }

        
    void GasKineticsWriter::installReagents(const vector_int& r,
        const vector_int& p, bool reversible) {
            
         int nr = r.size();
         int rnum = reactionNumber();
         int i;
         for (i = 0; i < nr; i++) {
             m_rrxn[r[i]][rnum] += 1.0;
         }

        m_reactantWriter.add( reactionNumber(), r);

         int np = p.size();

         for (i = 0; i < np; i++) {
             m_prxn[p[i]][rnum] += 1.0;
         }

        if (reversible) {
            m_revProductWriter.add(reactionNumber(), p);
            m_dn.push_back(np - nr);
            m_revindex.push_back(reactionNumber());
            m_nrev++;
        }
        else {
            m_irrevProductWriter.add(reactionNumber(), p);
            m_irrev.push_back( reactionNumber() );
            m_nirrev++;
        }        
    }

    void GasKineticsWriter::init(int nsp) {
        m_rrxn.resize(nsp);
        m_prxn.resize(nsp);
    }

}








