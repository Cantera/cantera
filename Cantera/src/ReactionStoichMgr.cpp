// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif


#include "ReactionStoichMgr.h"
#include "StoichManager.h"

namespace Cantera {

    ReactionStoichMgr::
    ReactionStoichMgr() {
        m_reactants = new StoichManagerN;
        m_revproducts = new StoichManagerN;
        m_irrevproducts = new StoichManagerN;
    }

    ReactionStoichMgr::~ReactionStoichMgr() {
        delete m_reactants;
        delete m_revproducts;
        delete m_irrevproducts;
    }

    void ReactionStoichMgr::
    add(int rxn, const vector_int& reactants, const vector_int& products,
        bool reversible) {
        vector_fp forder(reactants.size(), 1.0);
        add(rxn, reactants, products, reversible, forder);
    }

    void ReactionStoichMgr::
    add(int rxn, const vector_int& reactants, const vector_int& products,
        bool reversible, const vector_fp& fwdOrder) {
        m_reactants->add(rxn, reactants, fwdOrder);
        if (reversible) 
            m_revproducts->add(rxn, products);
        else
            m_irrevproducts->add(rxn, products);
    }

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
    }

    void ReactionStoichMgr::
    multiplyRevProducts(const doublereal* c, doublereal* r) {
        m_revproducts->multiply(c, r);
    }
}
