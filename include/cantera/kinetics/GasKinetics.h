/**
 * @file GasKinetics.h
 * @ingroup chemkinetics
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_GASKINETICS_H
#define CT_GASKINETICS_H

#include "BulkKinetics.h"
#include "FalloffMgr.h"
#include "Reaction.h"

namespace Cantera
{

/**
 * Kinetics manager for elementary gas-phase chemistry. This kinetics manager
 * implements standard mass-action reaction rate expressions for low-density
 * gases.
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
    GasKinetics(ThermoPhase* thermo = 0);

    virtual std::string kineticsType() const {
        return "Gas";
    }

    //! @}
    //! @name Reaction Rates Of Progress
    //! @{

    virtual void getEquilibriumConstants(doublereal* kc);
    virtual void getFwdRateConstants(double* kfwd);

    //! @}
    //! @name Reaction Mechanism Setup Routines
    //! @{
    virtual void init();
    virtual bool addReaction(shared_ptr<Reaction> r);
    virtual void modifyReaction(size_t i, shared_ptr<Reaction> rNew);
    virtual void invalidateCache();
    //@}

    virtual void finalizeSetup();
    virtual void updateROP();

    virtual void getJacobianSettings(AnyMap& settings) const;
    virtual void setJacobianSettings(const AnyMap& settings);
    virtual Eigen::SparseMatrix<double> fwdRatesOfProgress_ddC();
    virtual Eigen::SparseMatrix<double> revRatesOfProgress_ddC();
    virtual Eigen::SparseMatrix<double> netRatesOfProgress_ddC();
    virtual Eigen::VectorXd fwdRatesOfProgress_ddT();
    virtual Eigen::VectorXd revRatesOfProgress_ddT();
    virtual Eigen::VectorXd netRatesOfProgress_ddT();

    //! Update temperature-dependent portions of reaction rates and falloff
    //! functions.
    virtual void update_rates_T();

    //! Update properties that depend on concentrations.
    //! Currently the enhanced collision partner concentrations are updated
    //! here, as well as the pressure-dependent portion of P-log and Chebyshev
    //! reactions.
    virtual void update_rates_C();

private:
    //! @name Internal service methods
    /*!
     * These methods are for @internal use, and seek to avoid code duplication
     * while evaluating terms used for rate constants, rates of progress, and
     * their derivatives. Frequently used methods are defined in the header
     * file and use raw data pointers to allow inlining and avoid overhead.
     */
    //!@{

    //! Calculate rate coefficients
    void processFwdRateCoefficients(double* ropf)
    {
        update_rates_C();
        update_rates_T();

        // copy rate coefficients into ropf
        copy(m_rfn.begin(), m_rfn.end(), ropf);

        if (m_falloff_high_rates.nReactions()) {
            processFalloffReactions(ropf);
        }

        // Scale the forward rate coefficient by the perturbation factor
        for (size_t i = 0; i < nReactions(); ++i) {
            ropf[i] *= m_perturb[i];
        }
    }

    //! Multiply rate with third-body collider concentrations
    void processThirdBodies(double* rop)
    {
        // multiply rop by enhanced 3b conc for all 3b rxns
        if (!concm_3b_values.empty()) {
            m_3b_concm.multiply(rop, concm_3b_values.data());
        }

        // reactions involving third body
        if (!concm_multi_values.empty()) {
            m_multi_concm.multiply(rop, concm_multi_values.data());
        }
    }

    //! Multiply rate with inverse equilibrium constant
    void processEquilibriumConstants(double* rop)
    {
        // For reverse rates computed from thermochemistry, multiply the forward
        // rate coefficients by the reciprocals of the equilibrium constants
        for (size_t i = 0; i < nReactions(); ++i) {
            rop[i] *= m_rkcn[i];
        }
    }

    //! Multiply by scaling for mole-fraction-based species derivatives
    void processDensityConversion(double* rop);

    //! Multiply rate with scaled temperature derivatives of the inverse
    //! equilibrium constant
    /*!
     *  This (scaled) derivative is handled by a finite difference.
     */
    void processEquilibriumConstants_ddTscaled(double* drkcn);

    //! Routine to calculate numerical temperature derivatives
    /*!
     *  @TODO  This is a 'work-around', as there is no consistent handling
     *      of reaction rates yet (FalloffReaction and BlowersMaselReaction
     *      still use 'legacy' implementations). Once the transition to the
     *      'new' ReactionRate framework is complete, individual numerical
     *      derivatives should be handled there (potentially as a
     *      fall-back option if no exact derivative is available).
     */
    Eigen::VectorXd ratesOfProgress_ddT(bool forward);
    //!@}

    //! Buffers for partial rop results with length nReactions()
    vector_fp m_rbuf0;
    vector_fp m_rbuf1;
    vector_fp m_rbuf2;

    //! Jacobian settings
    bool m_jac_exact_temperature_derivatives;
    bool m_jac_skip_third_bodies;
    bool m_jac_skip_falloff;
    double m_jac_atol_deltaT;

protected:
    //! Reaction index of each falloff reaction
    std::vector<size_t> m_fallindx;

    //! Reaction index of each legacy reaction (old framework)
    std::vector<size_t> m_legacy;

    //! Map of reaction index to falloff reaction index (i.e indices in
    //! #m_falloff_low_rates and #m_falloff_high_rates)
    std::map<size_t, size_t> m_rfallindx;

    //! Rate expressions for falloff reactions at the low-pressure limit
    Rate1<Arrhenius> m_falloff_low_rates;

    //! Rate expressions for falloff reactions at the high-pressure limit
    Rate1<Arrhenius> m_falloff_high_rates;

    FalloffMgr m_falloffn;

    ThirdBodyCalc m_3b_concm;
    ThirdBodyCalc m_falloff_concm;

    Rate1<Plog> m_plog_rates;
    Rate1<Chebyshev> m_cheb_rates;
    Rate1<BlowersMasel> m_blowersmasel_rates;

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

    void processFalloffReactions(double* ropf);

    // functions marked as deprecated below are only used for XML import and
    // transitional reaction types are marked as '-legacy'

    //! @deprecated To be removed after Cantera 2.6 (replaced by MultiRate approach)
    void addThreeBodyReaction(ThreeBodyReaction2& r);
    void addFalloffReaction(FalloffReaction& r);
    //! @deprecated To be removed after Cantera 2.6 (replaced by MultiRate approach)
    void addPlogReaction(PlogReaction2& r);
    //! @deprecated To be removed after Cantera 2.6 (replaced by MultiRate approach)
    void addChebyshevReaction(ChebyshevReaction2& r);
    void addBlowersMaselReaction(BlowersMaselReaction& r);

    //! @deprecated To be removed after Cantera 2.6 (replaced by MultiRate approach)
    void modifyThreeBodyReaction(size_t i, ThreeBodyReaction2& r);
    void modifyFalloffReaction(size_t i, FalloffReaction& r);
    //! @deprecated To be removed after Cantera 2.6 (replaced by MultiRate approach)
    void modifyPlogReaction(size_t i, PlogReaction2& r);
    //! @deprecated To be removed after Cantera 2.6 (replaced by MultiRate approach)
    void modifyChebyshevReaction(size_t i, ChebyshevReaction2& r);
    void modifyBlowersMaselReaction(size_t i, BlowersMaselReaction& r);

    //! Update the equilibrium constants in molar units.
    void updateKc();
};

}

#endif
