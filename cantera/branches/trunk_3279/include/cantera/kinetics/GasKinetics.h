/**
 * @file GasKinetics.h
 *
 * @ingroup chemkinetics
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_GASKINETICS_H
#define CT_GASKINETICS_H

#include "BulkKinetics.h"
#include "ThirdBodyMgr.h"
#include "FalloffMgr.h"

namespace Cantera
{

/**
 * Kinetics manager for elementary gas-phase chemistry. This
 * kinetics manager implements standard mass-action reaction rate
 * expressions for low-density gases.
 * @ingroup kinetics
 */
class GasKinetics : public BulkKinetics
{
public:
    //! @name Constructors and General Information
    //! @{

    //! Constructor.
    /*!
     *  @param thermo  Pointer to the gas ThermoPhase (optional)
     */
    GasKinetics(thermo_t* thermo = 0);

    virtual Kinetics* duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const;

    virtual int type() const {
        return cGasKinetics;
    }

    //! @}
    //! @name Reaction Rates Of Progress
    //! @{

    virtual void getEquilibriumConstants(doublereal* kc);
    virtual void getFwdRateConstants(doublereal* kfwd);

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

    FalloffMgr                          m_falloffn;

    ThirdBodyMgr<Enhanced3BConc>        m_3b_concm;
    ThirdBodyMgr<Enhanced3BConc>        m_falloff_concm;

    Rate1<Plog> m_plog_rates;
    Rate1<ChebyshevRate> m_cheb_rates;

    //! @name Reaction rate data
    //!@{
    doublereal m_logp_ref;
    doublereal m_logc_ref;
    doublereal m_logStandConc;
    vector_fp m_rfn_low;
    vector_fp m_rfn_high;

    doublereal m_pres; //!< Last pressure at which rates were evaluated
    vector_fp falloff_work;
    vector_fp concm_3b_values;
    vector_fp concm_falloff_values;
    //!@}

    void processFalloffReactions();

    void addThreeBodyReaction(ReactionData& r);
    void addFalloffReaction(ReactionData& r);
    void addPlogReaction(ReactionData& r);
    void addChebyshevReaction(ReactionData& r);

    //! Update the equilibrium constants in molar units.
    void updateKc();

    bool m_finalized;
};
}

#endif
