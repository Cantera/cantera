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
        int n, ns, m, nr = r.reactants.size();
        for (n = 0; n < nr; n++) {
            ns = r.rstoich[n];
            //   m_rrxn[r.reactants[n]][rnum] = ns;
            for (m = 0; m < ns; m++) {
                rk.push_back(r.reactants[n]);
            }
        }
        
        vector_int pk;
        int np = r.products.size();
        for (n = 0; n < np; n++) {
            ns = r.pstoich[n];
            //  m_prxn[r.products[n]][rnum] = ns;
            for (m = 0; m < ns; m++) {
                pk.push_back(r.products[n]);
            }
        }

        m_reactants->add( rxn, rk);

        if (r.reversible) {
            m_revproducts->add(rxn, pk);
        }
        else {
            m_irrevproducts->add(rxn, pk);
        }

        if (r.global) {
            vector_fp delta_order(nr,0.0);
            for (n = 0; n < nr; n++) {
                delta_order[n] = r.order[n] - r.rstoich[n];
                //cout << "rxn stoich " << r.reactants[n] << "   " << r.order[n] << "  " << delta_order[n] << endl;
            }
            m_global->add(rxn, r.reactants, delta_order);
        }
  }


//   void ReactionStoichMgr::
//   add(int rxn, const vector_int& reactants, const vector_int& products,
//       bool reversible, const vector_fp& fwdOrder) {

// #ifdef DIAGNOSE_RXNSTOICHMGR
//             printf("ReactionStoichMgr::add    adding reaction number %d\n", rxn);
// #endif

//     // add the reactants with the specified forward order
//     m_reactants->add(rxn, reactants, fwdOrder);

//     // check whether any orders are not equal to 1.0
//     // if so, add an entry to the global stoich manager 
//     // with orders decremented by one
//     int nr = reactants.size();
//     int n;
//     bool global = false;
//     for (n = 0; n < nr; n++) {
//         if (fwdOrder[n] != 1.0) {

// #ifdef DIAGNOSE_RXNSTOICHMGR
//             printf(".... global reaction: fwdOrder[%d] = %f.\n", n, fwdOrder[n]);
// #endif
//             global = true;
//             if (reversible) {
//                 throw CanteraError("ReactionStoichMgr::add",
//                     "reversible global reactions not allowed");
//             }
//         }
//     }
//     if (global) {
//         vector_fp fwdOrder_minus_one;
//         for (n = 0; n < nr; n++) 
//             fwdOrder_minus_one[n] = fwdOrder[n] - 1.0;
//         m_global->add(rxn, reactants, fwdOrder_minus_one);
//     }


//     // depending on whether the reversible flag is set or not, add the
//     // products either to the reversible or irreversible product
//     // stoichiometry manager.
//     if (reversible) 
//       m_revproducts->add(rxn, products);
//     else
//       m_irrevproducts->add(rxn, products);
//   }


  void ReactionStoichMgr::
  getCreationRates(int nsp, const doublereal* ropf, const doublereal* ropr, doublereal* c) {
    fill(c, c + nsp, 0.0);
    m_revproducts->incrementSpecies(ropf, c);
    m_irrevproducts->incrementSpecies(ropf, c);
    m_reactants->incrementSpecies(ropr, c);
  }

    void ReactionStoichMgr::
    getDestructionRates(int nsp, const doublereal* ropf, const doublereal* ropr, doublereal* d) {
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
