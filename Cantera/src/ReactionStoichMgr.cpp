//------------------------------------------------
///
///  @file ReactionStoichMgr.cpp
///
///
//------------------------------------------------

// $Author$
// $Revision$
// $Date$

// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif


#include "ReactionStoichMgr.h"
#include "StoichManager.h"
#include "ctexceptions.h"
#include "diagnostics.h"
#include "ReactionData.h"

namespace Cantera {

  // create stoichiometry managers for the reactants of all reactions,
  // for the products of the reversible reactions, and for the
  // products of the irreversible reactions.
  ReactionStoichMgr::
  ReactionStoichMgr() {
    m_reactants = new StoichManagerN;
    m_revproducts = new StoichManagerN;
    m_irrevproducts = new StoichManagerN;
    m_global = new StoichManagerN;
    
  }

  // delete the three stoichiometry managers
  ReactionStoichMgr::~ReactionStoichMgr() {
    delete m_reactants;
    delete m_revproducts;
    delete m_irrevproducts;
    delete m_global;
  }


  void ReactionStoichMgr::
  add(int rxn, const vector_int& reactants, const vector_int& products,
      bool reversible) {

    m_reactants->add(rxn, reactants);

    if (reversible) 
      m_revproducts->add(rxn, products);
    else
      m_irrevproducts->add(rxn, products);
  }


  void ReactionStoichMgr::
  add(int rxn, const ReactionData& r) {

        vector_int rk;
        doublereal frac;
        int n, ns, m, nr = r.reactants.size();
        for (n = 0; n < nr; n++) {
            ns = int(r.rstoich[n]);
            frac = r.rstoich[n] - 1.0*int(r.rstoich[n]);
            if (frac != 0.0) {
                throw CanteraError("ReactionStoichMgr::add",
                    "fractional reactant stoichiometric coefficient not allowed!");
            }
            for (m = 0; m < ns; m++) {
                rk.push_back(r.reactants[n]);
            }
        }

        m_reactants->add( rxn, rk);
        
        vector_int pk;
        bool isfrac = false;
        int np = r.products.size();
        for (n = 0; n < np; n++) {
            ns = int(r.pstoich[n]);
            frac = r.pstoich[n] - 1.0*int(r.pstoich[n]);
            if (frac != 0.0) isfrac = true;
            for (m = 0; m < ns; m++) {
                pk.push_back(r.products[n]);
            }
        }

        if (isfrac) {            
            if (r.reversible) {
                throw CanteraError("ReactionStoichMgr::add",
                    "fractional product stoichiometric coefficients only allowed "
                    "\nfor irreversible reactions");
            }
            else {            
                cout << "adding fractional reaction" << endl;
                m_irrevproducts->add(rxn, r.products, r.order, r.pstoich);
            }
        }
        else if (r.reversible) {
            m_revproducts->add(rxn, pk);
        }
        else {   
            m_irrevproducts->add(rxn, pk);
        }

        // if the reaction is a global one, then add it to the global
        // stoichiometry manager. The reaction was added to
        // m_reactants above, but this will compute the rate using
        // stoichiometric coefficients. The global stoich manager then
        // corrects this, using orders that are the difference between
        // the true order and the reactant stoichiometric coefficient.
        if (r.global) {
            vector_fp delta_order(nr,0.0);
            for (n = 0; n < nr; n++) {
                delta_order[n] = r.order[n] - r.rstoich[n];
            }
            m_global->add(rxn, r.reactants, delta_order);
        }
  }

  void ReactionStoichMgr::
  getCreationRates(int nsp, const doublereal* ropf, 
      const doublereal* ropr, doublereal* c) {
      // zero out the output array
    fill(c, c + nsp, 0.0);
    
    m_revproducts->incrementSpecies(ropf, c);
    m_irrevproducts->incrementSpecies(ropf, c);

    // the reverse direction creates reactant species
    m_reactants->incrementSpecies(ropr, c);
  }

    void ReactionStoichMgr::
    getDestructionRates(int nsp, const doublereal* ropf, 
        const doublereal* ropr, doublereal* d) {
        fill(d, d + nsp, 0.0);
        m_revproducts->incrementSpecies(ropr, d);
        m_reactants->incrementSpecies(ropf, d);
    }

    void ReactionStoichMgr::
    getNetProductionRates(int nsp, const doublereal* ropnet, doublereal* w) {
        fill(w, w + nsp, 0.0);
        m_revproducts->incrementSpecies(ropnet, w);
        m_irrevproducts->incrementSpecies(ropnet, w);
        m_reactants->decrementSpecies(ropnet, w);
    }

    void ReactionStoichMgr::
    getReactionDelta(int nr, const doublereal* g, doublereal* dg) {
        fill(dg, dg + nr, 0.0);
        m_revproducts->incrementReactions(g, dg);
        m_irrevproducts->incrementReactions(g, dg);
        m_reactants->decrementReactions(g, dg);
    }

    void ReactionStoichMgr::
    getRevReactionDelta(int nr, const doublereal* g, doublereal* dg) {
        fill(dg, dg + nr, 0.0);
        m_revproducts->incrementReactions(g, dg);
        m_reactants->decrementReactions(g, dg);
    }

    void ReactionStoichMgr::
    multiplyReactants(const doublereal* c, doublereal* r) {
        m_reactants->multiply(c, r);
        m_global->power(c, r);
    }

    void ReactionStoichMgr::
    multiplyRevProducts(const doublereal* c, doublereal* r) {
        m_revproducts->multiply(c, r);
    }
}
