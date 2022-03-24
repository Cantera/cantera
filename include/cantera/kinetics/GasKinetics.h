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

    virtual void getThirdBodyConcentrations(double* concm);
    virtual const vector_fp& thirdBodyConcentrations() const {
        return m_concm;
    }

    //! @}
    //! @name Reaction Rates Of Progress
    //! @{

    virtual void getEquilibriumConstants(doublereal* kc);
    virtual void getFwdRateConstants(double* kfwd);

    //! @}
    //! @name Reaction Mechanism Setup Routines
    //! @{
    virtual bool addReaction(shared_ptr<Reaction> r, bool resize=true);
    virtual void modifyReaction(size_t i, shared_ptr<Reaction> rNew);
    virtual void invalidateCache();
    //! @}

    virtual void resizeReactions();
    void updateROP();

    virtual void getDerivativeSettings(AnyMap& settings) const;
    virtual void setDerivativeSettings(const AnyMap& settings);
    virtual void getFwdRateConstants_ddT(double* dkfwd);
    virtual void getFwdRatesOfProgress_ddT(double* drop);
    virtual void getRevRatesOfProgress_ddT(double* drop);
    virtual void getNetRatesOfProgress_ddT(double* drop);
    virtual void getFwdRateConstants_ddP(double* dkfwd);
    virtual void getFwdRatesOfProgress_ddP(double* drop);
    virtual void getRevRatesOfProgress_ddP(double* drop);
    virtual void getNetRatesOfProgress_ddP(double* drop);
    virtual void getFwdRateConstants_ddC(double* dkfwd);
    virtual void getFwdRatesOfProgress_ddC(double* drop);
    virtual void getRevRatesOfProgress_ddC(double* drop);
    virtual void getNetRatesOfProgress_ddC(double* drop);
    virtual Eigen::SparseMatrix<double> fwdRatesOfProgress_ddX();
    virtual Eigen::SparseMatrix<double> revRatesOfProgress_ddX();
    virtual Eigen::SparseMatrix<double> netRatesOfProgress_ddX();

    //! Update temperature-dependent portions of reaction rates and falloff
    //! functions.
    virtual void update_rates_T();

    //! Update properties that depend on concentrations.
    //! Currently the enhanced collision partner concentrations are updated
    //! here, as well as the pressure-dependent portion of P-log and Chebyshev
    //! reactions.
    virtual void update_rates_C();

protected:
    //! @name Internal service methods
    /*!
     * These methods are for @internal use, and seek to avoid code duplication
     * while evaluating terms used for rate constants, rates of progress, and
     * their derivatives.
     */
    //! @{

    //! Calculate rate coefficients
    void processFwdRateCoefficients(double* ropf);

    //! Multiply rate with third-body collider concentrations
    void processThirdBodies(double* rop);

    //! Multiply rate with inverse equilibrium constant
    void processEquilibriumConstants(double* rop);

    //! Multiply rate with scaled temperature derivatives of the inverse
    //! equilibrium constant
    /*!
     *  This (scaled) derivative is handled by a finite difference.
     */
    void processEquilibriumConstants_ddT(double* drkcn);

    //! Process temperature derivative
    //! @param in  rate expression used for the derivative calculation
    //! @param drop  pointer to output buffer
    void process_ddT(const vector_fp& in, double* drop);

    //! Process pressure derivative
    //! @param in  rate expression used for the derivative calculation
    //! @param drop  pointer to output buffer
    void process_ddP(const vector_fp& in, double* drop);

    //! Process concentration (molar density) derivative
    //! @param stoich  stoichiometry manager
    //! @param in  rate expression used for the derivative calculation
    //! @param drop  pointer to output buffer
    //! @param mass_action  boolean indicating whether law of mass action applies
    void process_ddC(StoichManagerN& stoich, const vector_fp& in,
                     double* drop, bool mass_action=true);

    //! Process mole fraction derivative
    //! @param stoich  stoichiometry manager
    //! @param in  rate expression used for the derivative calculation
    Eigen::SparseMatrix<double> process_ddX(StoichManagerN& stoich,
                                            const vector_fp& in);

    //! Helper function ensuring that all rate derivatives can be calculated
    //! @param name  method name used for error output
    //! @throw  CanteraError if legacy rates are present
    void assertDerivativesValid(const std::string& name);

    //! @}

    //! Reaction index of each falloff reaction
    std::vector<size_t> m_fallindx; //!< @deprecated (legacy only)

    //! Reaction index of each legacy reaction (old framework)
    std::vector<size_t> m_legacy;

    //! Map of reaction index to falloff reaction index (i.e indices in
    //! #m_falloff_low_rates and #m_falloff_high_rates)
    std::map<size_t, size_t> m_rfallindx; //!< @deprecated (legacy only)

    //! Rate expressions for falloff reactions at the low-pressure limit
    Rate1<Arrhenius2> m_falloff_low_rates; //!< @deprecated (legacy only)

    //! Rate expressions for falloff reactions at the high-pressure limit
    Rate1<Arrhenius2> m_falloff_high_rates; //!< @deprecated (legacy only)

    FalloffMgr m_falloffn; //!< @deprecated (legacy only)

    ThirdBodyCalc m_3b_concm; //!< @deprecated (legacy only)
    ThirdBodyCalc m_falloff_concm; //!< @deprecated (legacy only)

    Rate1<Plog> m_plog_rates; //!< @deprecated (legacy only)
    Rate1<ChebyshevRate> m_cheb_rates; //!< @deprecated (legacy only)

    //! @name Reaction rate data
    //!@{
    doublereal m_logStandConc;
    vector_fp m_rfn_low; //!< @deprecated (legacy only)
    vector_fp m_rfn_high; //!< @deprecated (legacy only)

    doublereal m_pres; //!< Last pressure at which rates were evaluated
    vector_fp falloff_work; //!< @deprecated (legacy only)
    vector_fp concm_3b_values; //!< @deprecated (legacy only)
    vector_fp concm_falloff_values; //!< @deprecated (legacy only)

    //!@}

    //! Buffers for partial rop results with length nReactions()
    vector_fp m_rbuf0;
    vector_fp m_rbuf1;
    vector_fp m_rbuf2;
    vector_fp m_sbuf0;
    vector_fp m_state;

    //! Derivative settings
    bool m_jac_skip_third_bodies;
    bool m_jac_skip_falloff;
    double m_jac_rtol_delta;

    // functions marked as deprecated below are only used for XML import and
    // transitional reaction types that are marked as '-legacy'

    //! @deprecated To be removed after Cantera 2.6 (replaced by MultiRate approach)
    void processFalloffReactions(double* ropf);

    //! @deprecated To be removed after Cantera 2.6 (replaced by MultiRate approach)
    void addThreeBodyReaction(ThreeBodyReaction2& r);
    //! @deprecated To be removed after Cantera 2.6 (replaced by MultiRate approach)
    void addFalloffReaction(FalloffReaction2& r);
    //! @deprecated To be removed after Cantera 2.6 (replaced by MultiRate approach)
    void addPlogReaction(PlogReaction2& r);
    //! @deprecated To be removed after Cantera 2.6 (replaced by MultiRate approach)
    void addChebyshevReaction(ChebyshevReaction2& r);

    //! @deprecated To be removed after Cantera 2.6 (replaced by MultiRate approach)
    void modifyThreeBodyReaction(size_t i, ThreeBodyReaction2& r);
    //! @deprecated To be removed after Cantera 2.6 (replaced by MultiRate approach)
    void modifyFalloffReaction(size_t i, FalloffReaction2& r);
    //! @deprecated To be removed after Cantera 2.6 (replaced by MultiRate approach)
    void modifyPlogReaction(size_t i, PlogReaction2& r);
    //! @deprecated To be removed after Cantera 2.6 (replaced by MultiRate approach)
    void modifyChebyshevReaction(size_t i, ChebyshevReaction2& r);

    //! Update the equilibrium constants in molar units.
    void updateKc();
};

}

#endif
