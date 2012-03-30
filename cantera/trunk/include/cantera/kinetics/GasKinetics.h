/**
 * @file GasKinetics.h
 *
 * @ingroup chemkinetics
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_GASKINETICS_H
#define CT_GASKINETICS_H

#include <fstream>
#include <map>

#include "cantera/thermo/mix_defs.h"
#include "Kinetics.h"

#include "cantera/base/utilities.h"

#include "ReactionStoichMgr.h"
#include "ThirdBodyMgr.h"
#include "FalloffMgr.h"
#include "RateCoeffMgr.h"

#include <cmath>
#include <cstdlib>

void get_wdot(const doublereal* rop, doublereal* wdot);

namespace Cantera
{

// forward references
class Enhanced3BConc;
class ReactionData;

/**
 * Kinetics manager for elementary gas-phase chemistry. This
 * kinetics manager implements standard mass-action reaction rate
 * expressions for low-density gases.
 * @ingroup kinetics
 */
class GasKinetics : public Kinetics
{

public:

    /**
     * @name Constructors and General Information
     */
    //@{

    //! Constructor.
    /*!
     *  @param thermo  Pointer to the gas ThermoPhase (optional)
     */
    GasKinetics(thermo_t* thermo = 0);


    //!Copy Constructor for the %GasKinetics object.
    /*!
     * Currently, this is not fully implemented. If called it will
     * throw an exception.
     *
     * @param right  object to be copied
     */
    GasKinetics(const GasKinetics& right);

    //! Destructor.
    virtual ~GasKinetics();

    //! Assignment operator
    /*!
     *  This is NOT a virtual function.
     *
     * @param right    Reference to %GasKinetics object to be copied into the
     *                 current one.
     */
    GasKinetics& operator=(const GasKinetics& right);

    //! Duplication routine for objects which inherit from Kinetics
    /*!
     *  This virtual routine can be used to duplicate %Kinetics objects
     *  inherited from %Kinetics even if the application only has
     *  a pointer to %Kinetics to work with.
     *
     *  These routines are basically wrappers around the derived copy  constructor.
     *
     * @param  tpVector Vector of shallow pointers to ThermoPhase objects. this is the
     *                  m_thermo vector within this object
     */
    virtual Kinetics* duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const;


    //! Identifies the subclass of the Kinetics manager type.
    /*!
     * These are listed in mix_defs.h.
     * @deprecated use type() instead
     */
    DEPRECATED(virtual int ID() const) {
        return cGasKinetics;
    }

    //!  Identifies the kinetics manager type.
    /*!
     *   Each class derived from Kinetics should overload this method to
     *   return a unique integer. Standard values are defined in file
     *   mix_defs.h.
     */
    virtual int type() const {
        return cGasKinetics;
    }

    virtual doublereal reactantStoichCoeff(size_t k, size_t i) const {
        return m_rrxn[k][i];
    }

    virtual doublereal productStoichCoeff(size_t k, size_t i) const {
        return m_prxn[k][i];
    }

    //@}
    /**
     * @name Reaction Rates Of Progress
     */
    //@{
    /**
     * Forward rates of progress.
     * Return the forward rates of progress in array fwdROP, which
     * must be dimensioned at least as large as the total number
     * of reactions.
     */
    virtual void getFwdRatesOfProgress(doublereal* fwdROP) {
        updateROP();
        std::copy(m_ropf.begin(), m_ropf.end(), fwdROP);
    }

    /**
     * Reverse rates of progress.
     * Return the reverse rates of progress in array revROP, which
     * must be dimensioned at least as large as the total number
     * of reactions.
     */
    virtual void getRevRatesOfProgress(doublereal* revROP) {
        updateROP();
        std::copy(m_ropr.begin(), m_ropr.end(), revROP);
    }

    /**
     * Net rates of progress.  Return the net (forward - reverse)
     * rates of progress in array netROP, which must be
     * dimensioned at least as large as the total number of
     * reactions.
     */
    virtual void getNetRatesOfProgress(doublereal* netROP) {
        updateROP();
        std::copy(m_ropnet.begin(), m_ropnet.end(), netROP);
    }


    /**
     * Equilibrium constants. Return the equilibrium constants of
     * the reactions in concentration units in array kc, which
     * must be dimensioned at least as large as the total number
     * of reactions.
     */
    virtual void getEquilibriumConstants(doublereal* kc);

    /**
     * Return the array of values for the reaction gibbs free energy
     * change.
     * These values depend on the species concentrations.
     *
     *  units = J kmol-1
     */
    virtual void getDeltaGibbs(doublereal* deltaG);

    /**
     * Return the array of values for the reaction enthalpy change.
     * These values depend upon the species concentrations.
     *
     *  units = J kmol-1
     */
    virtual void getDeltaEnthalpy(doublereal* deltaH);

    /**
     * Return the array of values for the reactions change in
     * entropy.
     * These values depend upon the concentration
     * of the solution.
     *
     *  units = J kmol-1 Kelvin-1
     */
    virtual void getDeltaEntropy(doublereal* deltaS);

    /**
     * Return the array of values for the reaction
     * standard state Gibbs free energy change.
     * These values do not depend on the species
     * concentrations.
     *
     *  units = J kmol-1
     */
    virtual void getDeltaSSGibbs(doublereal* deltaG);

    /**
     * Return the array of values for the change in the
     * standard state enthalpies of reaction.
     * These values do not depend upon the concentration
     * of the solution.
     *
     *  units = J kmol-1
     */
    virtual void getDeltaSSEnthalpy(doublereal* deltaH);

    /**
     * Return the array of values for the change in the
     * standard state entropies for each reaction.
     * These values do not depend upon the concentration
     * of the solution.
     *
     *  units = J kmol-1 Kelvin-1
     */
    virtual void getDeltaSSEntropy(doublereal* deltaS);

    //@}
    /**
     * @name Species Production Rates
     */
    //@{

    //! Return the species net production rates
    /*!
     * Species net production rates [kmol/m^3/s]. Return the species
     * net production rates (creation - destruction) in array
     * wdot, which must be dimensioned at least as large as the
     * total number of species.
     *
     *  @param net  Array of species production rates.
     *             units kmol m-3 s-1
     */
    virtual void getNetProductionRates(doublereal* net);

    //! Return the species creation rates
    /*!
     * Species creation rates [kmol/m^3]. Return the species
     * creation rates in array cdot, which must be
     * dimensioned at least as large as the total number of
     * species.
     *
     *  @param cdot  Array of species creation rates.
     *              units kmol m-3 s-1
     */
    virtual void getCreationRates(doublereal* cdot);

    //! Return a vector of the species destruction rates
    /*!
     * Species destruction rates [kmol/m^3]. Return the species
     * destruction rates in array ddot, which must be
     * dimensioned at least as large as the total number of
     * species.
     *
     *
     *  @param ddot  Array of species destruction rates.
     *               units kmol m-3 s-1
     *
     */
    virtual void getDestructionRates(doublereal* ddot);

    //@}
    /**
     * @name Reaction Mechanism Informational Query Routines
     */
    //@{

    /**
     * Flag specifying the type of reaction. The legal values and
     * their meaning are specific to the particular kinetics
     * manager.
     */
    virtual int reactionType(size_t i) const {
        return m_index[i].first;
    }

    virtual std::string reactionString(size_t i) const {
        return m_rxneqn[i];
    }

    /**
     * True if reaction i has been declared to be reversible. If
     * isReversible(i) is false, then the reverse rate of progress
     * for reaction i is always zero.
     */
    virtual bool isReversible(size_t i) {
        if (std::find(m_revindex.begin(), m_revindex.end(), i)
                < m_revindex.end()) {
            return true;
        } else {
            return false;
        }
    }

    /**
     * Return the forward rate constants
     *
     * length is the number of reactions. units depends
     * on many issues.
     */
    virtual void getFwdRateConstants(doublereal* kfwd);

    /**
     * Return the reverse rate constants.
     *
     * length is the number of reactions. units depends
     * on many issues. Note, this routine will return rate constants
     * for irreversible reactions if the default for
     * doIrreversible is overridden.
     */
    virtual void getRevRateConstants(doublereal* krev,
                                     bool doIrreversible = false);

    //@}
    /**
     * @name Reaction Mechanism Setup Routines
     */
    //@{

    virtual void init();

    ///  Add a reaction to the mechanism.
    virtual void addReaction(ReactionData& r);

    virtual void finalize();
    virtual bool ready() const;

    virtual void update_T();
    virtual void update_C();

    void updateROP();


    const std::vector<grouplist_t>& reactantGroups(size_t i) {
        return m_rgroups[i];
    }
    const std::vector<grouplist_t>& productGroups(size_t i) {
        return m_pgroups[i];
    }


    void _update_rates_T();

    //! Update properties that depend on concentrations.
    //! Currently the enhanced collision partner concentrations are updated
    //! here, as well as the pressure-dependent portion of P-log and Chebyshev
    //! reactions.
    void _update_rates_C();

    //@}

protected:

    size_t m_nfall;

    std::vector<size_t> m_fallindx;

    Rate1<Arrhenius>                    m_falloff_low_rates;
    Rate1<Arrhenius>                    m_falloff_high_rates;
    Rate1<Arrhenius>                    m_rates;

    mutable std::map<size_t, std::pair<int, size_t> > m_index;

    FalloffMgr                          m_falloffn;

    ThirdBodyMgr<Enhanced3BConc>        m_3b_concm;
    ThirdBodyMgr<Enhanced3BConc>        m_falloff_concm;

    std::vector<size_t> m_irrev;

    Rate1<Plog> m_plog_rates;
    Rate1<ChebyshevRate> m_cheb_rates;

    ReactionStoichMgr*                   m_rxnstoich;

    std::vector<size_t> m_fwdOrder;

    size_t m_nirrev;
    size_t m_nrev;

    std::map<size_t, std::vector<grouplist_t> > m_rgroups;
    std::map<size_t, std::vector<grouplist_t> > m_pgroups;

    std::vector<int>                         m_rxntype;

    mutable std::vector<std::map<size_t, doublereal> > m_rrxn;
    mutable std::vector<std::map<size_t, doublereal> > m_prxn;

    /**
     * Difference between the input global reactants order
     * and the input global products order. Changed to a double
     * to account for the fact that we can have real-valued
     * stoichiometries.
     */
    vector_fp  m_dn;
    std::vector<size_t> m_revindex;

    std::vector<std::string> m_rxneqn;

    //! @name Reaction rate data
    //!@{
    doublereal m_logp_ref;
    doublereal m_logc_ref;
    doublereal m_logStandConc;
    vector_fp m_ropf;
    vector_fp m_ropr;
    vector_fp m_ropnet;
    vector_fp m_rfn_low;
    vector_fp m_rfn_high;
    bool m_ROP_ok;

    doublereal m_temp;
    vector_fp m_rfn;
    vector_fp falloff_work;
    vector_fp concm_3b_values;
    vector_fp concm_falloff_values;
    vector_fp m_rkcn;
    //!@}

    vector_fp m_conc;
    void processFalloffReactions();
    vector_fp m_grt;


private:

    size_t reactionNumber() {
        return m_ii;
    }
    std::vector<std::map<int, doublereal> > m_stoich;

    void addElementaryReaction(ReactionData& r);
    void addThreeBodyReaction(ReactionData& r);
    void addFalloffReaction(ReactionData& r);
    void addPlogReaction(ReactionData& r);
    void addChebyshevReaction(ReactionData& r);

    void installReagents(const ReactionData& r);

    void installGroups(size_t irxn, const std::vector<grouplist_t>& r,
                       const std::vector<grouplist_t>& p);
    void updateKc();

    void registerReaction(size_t rxnNumber, int type, size_t loc) {
        m_index[rxnNumber] = std::pair<int, size_t>(type, loc);
    }
    bool m_finalized;
};
}

#endif
