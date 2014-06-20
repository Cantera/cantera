/**
 * @file GasKinetics.h
 *
 * @ingroup chemkinetics
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_GASKINETICS_H
#define CT_GASKINETICS_H

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
    //! @name Constructors and General Information
    //! @{

    //! Constructor.
    /*!
     *  @param thermo  Pointer to the gas ThermoPhase (optional)
     */
    GasKinetics(thermo_t* thermo = 0);

    //! Copy Constructor
    GasKinetics(const GasKinetics& right);

    //! Assignment operator
    GasKinetics& operator=(const GasKinetics& right);

    virtual Kinetics* duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const;

    virtual int type() const {
        return cGasKinetics;
    }

    virtual doublereal reactantStoichCoeff(size_t k, size_t i) const {
        return m_rrxn[k][i];
    }

    virtual doublereal productStoichCoeff(size_t k, size_t i) const {
        return m_prxn[k][i];
    }

    //! @}
    //! @name Reaction Rates Of Progress
    //! @{

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

    virtual void getNetProductionRates(doublereal* net);
    virtual void getCreationRates(doublereal* cdot);
    virtual void getDestructionRates(doublereal* ddot);

    //! @}
    //! @name Reaction Mechanism Informational Query Routines
    //! @{

    virtual int reactionType(size_t i) const {
        return m_index[i].first;
    }

    virtual std::string reactionString(size_t i) const {
        return m_rxneqn[i];
    }

    virtual std::string reactantString(size_t i) const {
        return m_reactantStrings[i];
    }

    virtual std::string productString(size_t i) const {
        return m_productStrings[i];
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
    //@}

    void updateROP();

    const std::vector<grouplist_t>& reactantGroups(size_t i) {
        return m_rgroups[i];
    }
    const std::vector<grouplist_t>& productGroups(size_t i) {
        return m_pgroups[i];
    }

    //! Update temperature-dependent portions of reaction rates and falloff
    //! functions.
    virtual void update_rates_T();

    //! Update properties that depend on concentrations.
    //! Currently the enhanced collision partner concentrations are updated
    //! here, as well as the pressure-dependent portion of P-log and Chebyshev
    //! reactions.
    virtual void update_rates_C();

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

    ReactionStoichMgr m_rxnstoich;

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
    std::vector<std::string> m_reactantStrings;
    std::vector<std::string> m_productStrings;

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
    doublereal m_pres; //!< Last pressure at which rates were evaluated
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

    //! Update the equilibrium constants in molar units.
    void updateKc();

    void registerReaction(size_t rxnNumber, int type_, size_t loc) {
        m_index[rxnNumber] = std::pair<int, size_t>(type_, loc);
    }
    bool m_finalized;
};
}

#endif
