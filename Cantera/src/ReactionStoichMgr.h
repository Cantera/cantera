/**
 * @file ReactionStoichMgr.h
 *
 * $Author$
 * $Revision$
 * $Date$
 */

#ifndef CT_RXN_STOICH
#define CT_RXN_STOICH


#include "ct_defs.h"

namespace Cantera {

    class StoichManagerN;

    /**
     * This class handles calculations involving reaction stoichiometry.
     */
    class ReactionStoichMgr {

    public:

        ReactionStoichMgr();
        virtual ~ReactionStoichMgr();

        /**
         * Add a reaction with specified, possibly non-integral, reaction orders. 
         * @param rxn Reaction number
         * @param reactants vector of integer reactant indices
         * @param reactants vector of integer product indices
         * @param reversible true if the reaction is reversible, false otherwise
         * @param fwdOrder reaction orders for the reactants. This vector must
         * be the same length as 'reactants,' and the reaction orders are for the
         * species with index in the corresponding location in 'reactants.'
         */
        void add(int rxn, const vector_int& reactants, const vector_int& products,
            bool reversible, const vector_fp& fwdOrder);

        /**
         * Add a reaction with mass-action kinetics.
         * @param rxn Reaction number
         * @param reactants vector of integer reactant indices
         * @param reactants vector of integer product indices
         * @param reversible true if the reaction is reversible, false otherwise
         */
        void add(int rxn, const vector_int& reactants, const vector_int& products,
            bool reversible); 

        /**
         * Given the arrays of the forward and reverse rates of progress for all reactions,
         * compute the species creation rates and return them in array c.
         */
        void getCreationRates(int nsp, const doublereal* ropf, const doublereal* ropr, doublereal* c);

        /**
         * Given the arrays of the forward and reverse rates of progress for all reactions,
         * compute the species destruction rates and return them in array d.
         */
        void getDestructionRates(int nsp, const doublereal* ropf, const doublereal* ropr, doublereal* d);

        /**
         * Given the array of the net rates of progress for all reactions,
         * compute the species net production rates and return them in array w.
         */
        void getNetProductionRates(int nsp, const doublereal* ropnet, doublereal* w);

        /**
         * Given an array of species properties 'g', return in array 'dg' the change in this quantity 
         * in the reactions. Array 'g' must have a length at least as great
         * as the number of species, and array 'dg' must have a length
         * as great as the total number of reactions.  
         */
        void getReactionDelta(int nr, const doublereal* g, doublereal* dg);

        /**
         * Given an array of species properties 'g', return in array
         * 'dg' the change in this quantity in the reversible
         * reactions. Array 'g' must have a length at least as great
         * as the number of species, and array 'dg' must have a length
         * as great as the total number of reactions.  This method
         * only computes 'dg' for the reversible reactions, and the
         * entries of 'dg' for the irreversible reactions are
         * unaltered. This is primarily designed for use in
         * calculating reveerse rate coefficients from thermochemistry
         * for reversible reactions.
         */
        void getRevReactionDelta(int nr, const doublereal* g, doublereal* dg);

        /** 
         * Given an array of concentrations C, multiply the entries in array R by 
         * the concentration products for the reactants:
         * \f[
         *  R_i = R_i * \prod_k C_k^{o_{k,i}}
         * \f]
         * Here \f$ o_{k,i} \f$ is the reaction order of species k in reaction i.
         */
        void multiplyReactants(const doublereal* C, doublereal* R);

       /** 
         * Given an array of concentrations C, multiply the entries in array R by 
         * the concentration products for the products:
         * \f[
         *  R_i = R_i * \prod_k C_k^{\nu^{(p)}_{k,i}}
         * \f]
         * Here \f$ \nu^{(p)}_{k,i} \f$ is the product stoichiometric coefficient
         * of species k in reaction i.
         */
        void multiplyRevProducts(const doublereal* c, doublereal* r);

    protected:

        StoichManagerN*  m_reactants;
        StoichManagerN*  m_revproducts;
        StoichManagerN*  m_irrevproducts;
    };
}

#endif
