#ifndef CT_RXN_STOICH
#define CT_RXN_STOICH

#include "ct_defs.h"

namespace Cantera {

    class StoichManagerN;

    class ReactionStoichMgr {

    public:

        ReactionStoichMgr();
        virtual ~ReactionStoichMgr();

        void add(int rxn, const vector_int& reactants, const vector_int& products,
            bool reversible, const vector_fp& fwdOrder); 
        void add(int rxn, const vector_int& reactants, const vector_int& products,
            bool reversible); 

        void getCreationRates(int nsp, const doublereal* ropf, const doublereal* ropr, doublereal* c);
        void getDestructionRates(int nsp, const doublereal* ropf, const doublereal* ropr, doublereal* d);
        void getNetProductionRates(int nsp, const doublereal* ropnet, doublereal* w);
        void getReactionDelta(int nr, const doublereal* g, doublereal* dg);
        void getRevReactionDelta(int nr, const doublereal* g, doublereal* dg);

        void multiplyReactants(const doublereal* c, doublereal* r);
        void multiplyRevProducts(const doublereal* c, doublereal* r);

    protected:

        StoichManagerN*  m_reactants;
        StoichManagerN*  m_revproducts;
        StoichManagerN*  m_irrevproducts;
    };
}

#endif
