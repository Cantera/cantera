/**
 * @file ReactionStoichMgr.h
 *
 * Header file declaring class ReactionStoichMgr.
 */
#ifndef CT_RXN_STOICH
#define CT_RXN_STOICH

#include "cantera/kinetics/StoichManager.h"

namespace Cantera
{

class ReactionData;
class Reaction;

/**
 * Reaction mechanism stoichiometry manager. This is an internal class used
 * by kinetics manager classes, and is not meant for direct use in
 * user programs.
 *
 * Class ReactionStoichMgr handles the calculation of quantities involving
 * the stoichiometry of a set of reactions. The reactions may have integer
 * or non-integer stoichiometric coefficients. Specifically, its methods compute
 * - species creation rates
 * - species destruction rates
 * - species net production rates
 * - the change in molar species properties in the reactions
 * - concentration products
 *
 * To use this class, method add() is first used to add each reaction.
 * Once all reactions have been added, the methods that compute various
 * quantities may be called.
 *
 * The nomenclature used below to document the methods is as follows.
 *
 * - \f$ N_r \f$
 *     Integer reactant stoichiometric coefficient matrix. The (k,i)
 *     element of this matrix is the stoichiometric coefficient of
 *     species \e k as a reactant in reaction \e i.
 * - \f$ N_p \f$
 *     Integer product stoichiometric coefficient matrix. The (k,i)
 *     element of this matrix is the stoichiometric coefficient of
 *     species \e k as a product in reaction \e i.
 * - \f$ Q_{\rm fwd} \f$
 *     Vector of length I of forward rates of progress.
 * - \f$ Q_{\rm rev} \f$
 *     Vector of length I of reverse rates of progress.
 * - \f$ C \f$
 *     Vector of K species creation rates.
 * - \f$ D \f$
 *     Vector of K species destruction rates.
 * - \f$ W = C - D \f$
 *     Vector of K species net production rates.
 * @deprecated Unused; Functionality merged into class Kinetics. To be removed
 *     after Cantera 2.2.
 */
class ReactionStoichMgr
{
public:
    /// Constructor.
    ReactionStoichMgr();

    /// Destructor.
    virtual ~ReactionStoichMgr() {}

    ReactionStoichMgr(const ReactionStoichMgr& right);

    ReactionStoichMgr& operator=(const ReactionStoichMgr& right);

    
    //! Add a reaction with mass-action kinetics. 
    /*!Vectors
     * 'reactants' and 'products' contain the integer species
     * indices of the reactants and products, respectively.  Note
     * that if more than one molecule of a given species is
     * involved in the reaction, then its index is repeated.
     *
     * For example, suppose a reaction mechanism involves the
     * species N2, O2, O, N, NO. N2 is assigned index number 0, O2
     * number 1, and so on through NO with number 4.  Then the
     * representation of the following reactions is as shown here.
     *
     * - N + O = NO
     *   - reactants: (3, 2)
     *   - products:  (4)
     *
     * - O + O = O2
     *   - reactants: (2, 2)   [ note repeated index ]
     *   - products:  (1)
     *
     * @param rxn           Reaction number. This number will be used as the index
     *                      into the rate of progress vector in the methods below.
     * @param reactants     Vector of integer reactant indices
     * @param products      Vector of integer product indices
     * @param reversible    True if the reaction is reversible, false otherwise
     */
    virtual void add(size_t rxn, const std::vector<size_t>& reactants,
                     const std::vector<size_t>& products, bool reversible);

    /**
     * Add a reaction with specified, possibly non-integral, reaction orders.
     * @param rxn Reaction number
     * @param r Data structure containing reactant and product vectors, etc.
     */
    virtual void add(size_t rxn, const ReactionData& r);

    /**
     * Species creation rates. Given the arrays of the forward and reverse
     * rates of progress for all reactions, compute the species creation
     * rates, given by
     * \f[
     *  C = N_p Q_f  + N_r Q_r.
     * \f]
     */
    virtual void getCreationRates(size_t nSpecies,
                                  const doublereal* fwdRatesOfProgress,
                                  const doublereal* revRatesOfProgress,
                                  doublereal* creationRates);

    /**
     * Species destruction rates. Given the arrays of the forward and reverse
     * rates of progress for all reactions, compute the species destruction
     * rates, given by
     * \f[
     *  D = N_r Q_f  + N_p Q_r,
     * \f]
     * Note that the stoichiometric coefficient matrices are very sparse, integer
     * matrices.
     */
    virtual void getDestructionRates(size_t nSpecies,
                                     const doublereal* fwdRatesOfProgress,
                                     const doublereal* revRatesOfProgress,
                                     doublereal* destructionRates);

    /**
     * Species net production rates. Given the array of the net rates of
     * progress for all reactions, compute the species net production rates,
     * given by
     * \f[
     *  W = (N_r - N_p) Q_{\rm net},
     * \f]
     */
    virtual void getNetProductionRates(size_t nsp, const doublereal* ropnet, doublereal* w);

    //! Calculates the change of a molar species property in a reaction.
    /*!
     * Given an array of species properties 'g', return in array 'dg' the
     * change in this quantity in the reactions. Array 'g' must have a length
     * at least as great as the number of species, and array 'dg' must have a
     * length as great as the total number of reactions.
     *  \f[
     *      \delta g_i = \sum_k{\nu_{i,k} g_k   }
     *  \f]
     *
     * @param nReactions  Number of reactions
     * @param g           Molar property of the species.
     *                    An example would be the partial molar enthalpy
     *                    Length is equal to number of kinetic species
     * @param dg          Calculated property change of the reaction.
     *                    An example would be the delta change in enthalpy,
     *                    i.e., the enthalpy of reaction.
     */
    virtual void getReactionDelta(size_t nReactions,
                                  const doublereal* g,
                                  doublereal* dg);

    /**
     * Given an array of species properties 'g', return in array 'dg' the
     * change in this quantity in the reversible reactions. Array 'g' must
     * have a length at least as great as the number of species, and array
     * 'dg' must have a length as great as the total number of reactions.
     * This method only computes 'dg' for the reversible reactions, and the
     * entries of 'dg' for the irreversible reactions are unaltered. This is
     * primarily designed for use in calculating reverse rate coefficients
     * from thermochemistry for reversible reactions.
     */
    virtual void getRevReactionDelta(size_t nr, const doublereal* g, doublereal* dg);

    /**
     * Given an array of concentrations C, multiply the entries in array R by
     * the concentration products for the reactants.
     * \f[
     *  R_i = R_i * \prod_k C_k^{o_{k,i}}
     * \f]
     *
     * Here \f$ o_{k,i} \f$ is the reaction order of species k in reaction i.
     */
    virtual void multiplyReactants(const doublereal* C, doublereal* R);

    /**
     * Given an array of concentrations C, multiply the entries in array R by
     * the concentration products for the products.
     * \f[
     *  R_i = R_i * \prod_k C_k^{\nu^{(p)}_{k,i}}
     * \f]
     * Here \f$ \nu^{(p)}_{k,i} \f$ is the product stoichiometric coefficient
     * of species k in reaction i.
     */
    virtual void multiplyRevProducts(const doublereal* c, doublereal* r);

    //! @deprecated To be removed after Cantera 2.2
    virtual void write(const std::string& filename);

protected:
    //! @deprecated To be removed after Cantera 2.2
    void writeCreationRates(std::ostream& f);
    //! @deprecated To be removed after Cantera 2.2
    void writeDestructionRates(std::ostream& f);
    //! @deprecated To be removed after Cantera 2.2
    void writeNetProductionRates(std::ostream& f);
    //! @deprecated To be removed after Cantera 2.2
    void writeMultiplyReactants(std::ostream& f);
    //! @deprecated To be removed after Cantera 2.2
    void writeMultiplyRevProducts(std::ostream& f);
    StoichManagerN m_reactants;
    StoichManagerN m_revproducts;
    StoichManagerN m_irrevproducts;
    vector_fp m_dummy;
};
}

#endif
