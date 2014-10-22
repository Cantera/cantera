/**
 * @file GasKinetics.h
 *
 * @ingroup chemkinetics
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_GASKINETICS_H
#define CT_GASKINETICS_H

#include "Kinetics.h"
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
        return getValue(m_rrxn[k], i, 0.0);
    }

    virtual doublereal productStoichCoeff(size_t k, size_t i) const {
        return getValue(m_prxn[k], i, 0.0);
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

    std::map<size_t, std::pair<int, size_t> > m_index;

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

    std::vector<std::map<size_t, doublereal> > m_rrxn;
    std::vector<std::map<size_t, doublereal> > m_prxn;

    /**
     * Difference between the input global reactants order
     * and the input global products order. Changed to a double
     * to account for the fact that we can have real-valued
     * stoichiometries.
     */
    vector_fp  m_dn;
    std::vector<size_t> m_revindex;

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

    virtual void installReagents(const ReactionData& r);

    //! Update the equilibrium constants in molar units.
    void updateKc();

    void registerReaction(size_t rxnNumber, int type_, size_t loc) {
        m_index[rxnNumber] = std::pair<int, size_t>(type_, loc);
    }
    bool m_finalized;
};
}

#endif
