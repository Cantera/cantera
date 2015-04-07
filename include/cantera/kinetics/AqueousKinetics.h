/**
 * @file AqueousKinetics.h
 *
 * @ingroup chemkinetics
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_AQUEOUSKINETICS_H
#define CT_AQUEOUSKINETICS_H

#include "cantera/thermo/mix_defs.h"
#include "Kinetics.h"

#include "cantera/base/utilities.h"

#include "ReactionStoichMgr.h"
#include "ThirdBodyMgr.h"
#include "FalloffMgr.h"
#include "RateCoeffMgr.h"

namespace Cantera
{

// forward references

class ReactionData;

/**
 * Kinetics manager for elementary aqueous-phase chemistry. This
 * kinetics manager implements standard mass-action reaction rate
 * expressions for liquids
 *
 *
 *   Concentration
 *
 * @ingroup kinetics
 * @deprecated Not actually implemented
 */
class AqueousKinetics : public Kinetics
{

public:

    //! @name Constructors
    //! @{

    /// Constructor. Creates an empty reaction mechanism.
    AqueousKinetics(thermo_t* thermo = 0);

    AqueousKinetics(const AqueousKinetics& right);

    AqueousKinetics& operator=(const AqueousKinetics& right);

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
    //@}

    virtual int type() const {
        return cAqueousKinetics;
    }

    virtual doublereal reactantStoichCoeff(size_t k, size_t i) const {
        return m_rrxn[k][i];
    }

    virtual doublereal productStoichCoeff(size_t k, size_t i) const {
        return m_prxn[k][i];
    }

    //! @name Reaction Rates Of Progress
    //@{
    virtual void getFwdRatesOfProgress(doublereal* fwdROP) {
        updateROP();
        std::copy(m_ropf.begin(), m_ropf.end(), fwdROP);
    }

    virtual void getRevRatesOfProgress(doublereal* revROP) {
        updateROP();
        std::copy(m_ropr.begin(), m_ropr.end(), revROP);
    }

    virtual void getNetRatesOfProgress(doublereal* netROP) {
        updateROP();
        std::copy(m_ropnet.begin(), m_ropnet.end(), netROP);
    }

    virtual void getEquilibriumConstants(doublereal* kc);

    virtual void getDeltaGibbs(doublereal* deltaG);
    virtual void getDeltaEnthalpy(doublereal* deltaH);
    virtual void getDeltaEntropy(doublereal* deltaS);

    virtual void getDeltaSSGibbs(doublereal* deltaG);
    virtual void getDeltaSSEnthalpy(doublereal* deltaH);
    virtual void getDeltaSSEntropy(doublereal* deltaS);

    //! @}
    //! @name Species Production Rates
    //! @{

    virtual void getNetProductionRates(doublereal* net) {
        updateROP();
        m_rxnstoich.getNetProductionRates(m_kk, &m_ropnet[0], net);
    }

    virtual void getCreationRates(doublereal* cdot) {
        updateROP();
        m_rxnstoich.getCreationRates(m_kk, &m_ropf[0], &m_ropr[0], cdot);
    }

    virtual void getDestructionRates(doublereal* ddot) {
        updateROP();
        m_rxnstoich.getDestructionRates(m_kk, &m_ropf[0], &m_ropr[0], ddot);
    }

    //! @}
    //! @name Reaction Mechanism Informational Query Routines
    //! @{

    virtual int reactionType(size_t i) const {
        return m_index[i].first;
    }

    virtual std::string reactionString(size_t i) const {
        return m_rxneqn[i];
    }

    virtual bool isReversible(size_t i) {
        if (std::find(m_revindex.begin(), m_revindex.end(), i)
                < m_revindex.end()) {
            return true;
        } else {
            return false;
        }
    }

    virtual void getFwdRateConstants(doublereal* kfwd);
    virtual void getRevRateConstants(doublereal* krev,
                                     bool doIrreversible = false);

    //! @}
    //! @name Reaction Mechanism Setup Routines
    //! @{

    virtual void init();
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

    /*!
     * Update temperature-dependent portions of reaction rates and
     * falloff functions.
     */
    void _update_rates_T();

    /*!
     * Update properties that depend on concentrations. Currently only
     * the enhanced collision partner concentrations are updated here.
     */
    void _update_rates_C();

    //@}

protected:

    size_t m_nfall;

    Rate1<Arrhenius>                    m_rates;

    mutable std::map<size_t, std::pair<int, size_t> >   m_index;

    std::vector<size_t> m_irrev;

    ReactionStoichMgr m_rxnstoich;

    std::vector<size_t> m_fwdOrder;

    size_t m_nirrev;
    size_t m_nrev;

    std::map<size_t, std::vector<grouplist_t> > m_rgroups;
    std::map<size_t, std::vector<grouplist_t> > m_pgroups;

    std::vector<int>                         m_rxntype;

    mutable std::vector<std::map<size_t, doublereal> >     m_rrxn;
    mutable std::vector<std::map<size_t, doublereal> >     m_prxn;

    /**
     * Difference between the input global reactants order
     * and the input global products order. Changed to a double
     * to account for the fact that we can have real-valued
     * stoichiometries.
     */
    vector_fp  m_dn;
    std::vector<size_t> m_revindex;

    std::vector<std::string> m_rxneqn;

    vector_fp m_conc;
    vector_fp m_grt;

    //! @name Aqueous kinetics data
    //!@{
    vector_fp m_ropf;
    vector_fp m_ropr;
    vector_fp m_ropnet;
    bool m_ROP_ok;

    doublereal m_temp;
    vector_fp  m_rfn;

    vector_fp m_rkcn;
    //!@}

private:

    size_t reactionNumber() {
        return m_ii;
    }
    std::vector<std::map<int, doublereal> > m_stoich;

    void addElementaryReaction(ReactionData& r);

    void installReagents(const ReactionData& r);

    void installGroups(size_t irxn, const std::vector<grouplist_t>& r,
                       const std::vector<grouplist_t>& p);

    /**
     * Update the equilibrium constants in molar units.
     */
    void updateKc();

    void registerReaction(size_t rxnNumber, int type, size_t loc) {
        m_index[rxnNumber] = std::pair<int, size_t>(type, loc);
    }
    bool m_finalized;
};
}

#endif
