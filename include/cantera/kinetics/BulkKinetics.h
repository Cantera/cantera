/**
 * @file BulkKinetics.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_BULKKINETICS_H
#define CT_BULKKINETICS_H

#include "Kinetics.h"
#include "MultiRate.h"
#include "ThirdBodyCalc.h"

namespace Cantera
{

//! Specialization of Kinetics for chemistry in a single bulk phase
//! @ingroup kineticsmgr
class BulkKinetics : public Kinetics
{
public:
    //! @name Constructors and General Information
    //! @{
    BulkKinetics();

    string kineticsType() const override {
        return "bulk";
    }

    bool isReversible(size_t i) override;
    //! @}

    //! @name Reaction Mechanism Setup Routines
    //! @{
    bool addReaction(shared_ptr<Reaction> r, bool resize=true) override;
    void addThirdBody(shared_ptr<Reaction> r);
    void modifyReaction(size_t i, shared_ptr<Reaction> rNew) override;
    void resizeSpecies() override;
    void resizeReactions() override;
    void setMultiplier(size_t i, double f) override;
    void invalidateCache() override;
    //! @}

    //! @name Reaction rate constants, rates of progress, and thermodynamic properties
    //! @{
    void getFwdRateConstants(double* kfwd) override;
    void getEquilibriumConstants(double* kc) override;
    void getRevRateConstants(double* krev, bool doIrreversible=false) override;

    void getDeltaGibbs(double* deltaG) override;
    void getDeltaEnthalpy(double* deltaH) override;
    void getDeltaEntropy(double* deltaS) override;

    void getDeltaSSGibbs(double* deltaG) override;
    void getDeltaSSEnthalpy(double* deltaH) override;
    void getDeltaSSEntropy(double* deltaS) override;
    //! @}

    //! @name Derivatives of rate constants and rates of progress
    //! @{
    void getDerivativeSettings(AnyMap& settings) const override;
    void setDerivativeSettings(const AnyMap& settings) override;
    void getFwdRateConstants_ddT(double* dkfwd) override;
    void getFwdRatesOfProgress_ddT(double* drop) override;
    void getRevRatesOfProgress_ddT(double* drop) override;
    void getNetRatesOfProgress_ddT(double* drop) override;
    void getFwdRateConstants_ddP(double* dkfwd) override;
    void getFwdRatesOfProgress_ddP(double* drop) override;
    void getRevRatesOfProgress_ddP(double* drop) override;
    void getNetRatesOfProgress_ddP(double* drop) override;
    void getFwdRateConstants_ddC(double* dkfwd) override;
    void getFwdRatesOfProgress_ddC(double* drop) override;
    void getRevRatesOfProgress_ddC(double* drop) override;
    void getNetRatesOfProgress_ddC(double* drop) override;
    Eigen::SparseMatrix<double> fwdRatesOfProgress_ddX() override;
    Eigen::SparseMatrix<double> revRatesOfProgress_ddX() override;
    Eigen::SparseMatrix<double> netRatesOfProgress_ddX() override;
    Eigen::SparseMatrix<double> fwdRatesOfProgress_ddCi() override;
    Eigen::SparseMatrix<double> revRatesOfProgress_ddCi() override;
    Eigen::SparseMatrix<double> netRatesOfProgress_ddCi() override;
    //! @}

    //! @name Rate calculation intermediate methods
    //! @{

    void updateROP() override;

    void getThirdBodyConcentrations(double* concm) override;
    const vector<double>& thirdBodyConcentrations() const override {
        return m_concm;
    }

    //! @}

protected:
    //! @name Internal service methods
    //!
    //! @note These methods are for internal use, and seek to avoid code duplication
    //! while evaluating terms used for rate constants, rates of progress, and
    //! their derivatives.
    //! @{

    //! Multiply rate with third-body collider concentrations
    void processThirdBodies(double* rop);

    //! Multiply rate with inverse equilibrium constant
    void applyEquilibriumConstants(double* rop);

    //! Multiply rate with scaled temperature derivatives of the inverse
    //! equilibrium constant
    /*!
     *  This (scaled) derivative is handled by a finite difference.
     */
    void applyEquilibriumConstants_ddT(double* drkcn);

    //! Process temperature derivative
    //! @param in  rate expression used for the derivative calculation
    //! @param drop  pointer to output buffer
    void process_ddT(const vector<double>& in, double* drop);

    //! Process pressure derivative
    //! @param in  rate expression used for the derivative calculation
    //! @param drop  pointer to output buffer
    void process_ddP(const vector<double>& in, double* drop);

    //! Process concentration (molar density) derivative
    //! @param stoich  stoichiometry manager
    //! @param in  rate expression used for the derivative calculation
    //! @param drop  pointer to output buffer
    //! @param mass_action  boolean indicating whether law of mass action applies
    void process_ddC(StoichManagerN& stoich, const vector<double>& in,
                     double* drop, bool mass_action=true);

    //! Process derivatives
    //! @param stoich  stoichiometry manager
    //! @param in  rate expression used for the derivative calculation
    //! @param ddX true: w.r.t mole fractions false: w.r.t species concentrations
    //! @return a sparse matrix of derivative contributions for each reaction of
    //! dimensions nTotalReactions by nTotalSpecies
    Eigen::SparseMatrix<double> calculateCompositionDerivatives(
        StoichManagerN& stoich, const vector<double>& in, bool ddX=true);

    //! Helper function ensuring that all rate derivatives can be calculated
    //! @param name  method name used for error output
    //! @throw CanteraError if ideal gas assumption does not hold
    void assertDerivativesValid(const string& name);

    //! @}

    //! Vector of rate handlers
    vector<unique_ptr<MultiRateBase>> m_bulk_rates;
    map<string, size_t> m_bulk_types; //!< Mapping of rate handlers

    vector<size_t> m_revindex; //!< Indices of reversible reactions
    vector<size_t> m_irrev; //!< Indices of irreversible reactions

    //! Difference between the global reactants order and the global products
    //! order. Of type "double" to account for the fact that we can have real-
    //! valued stoichiometries.
    vector<double> m_dn;

    ThirdBodyCalc m_multi_concm; //!< used with MultiRate evaluator

    //! Third body concentrations
    vector<double> m_concm;

    //! Activity concentrations, as calculated by ThermoPhase::getActivityConcentrations
    vector<double> m_act_conc;

    //! Physical concentrations, as calculated by ThermoPhase::getConcentrations
    vector<double> m_phys_conc;

    //! Derivative settings
    bool m_jac_skip_third_bodies;
    bool m_jac_skip_falloff;
    double m_jac_rtol_delta;

    bool m_ROP_ok = false;

    //! Buffers for partial rop results with length nReactions()
    vector<double> m_rbuf0;
    vector<double> m_rbuf1;
    vector<double> m_rbuf2;
    vector<double> m_kf0; //!< Forward rate constants without perturbation
    vector<double> m_sbuf0;
    vector<double> m_state;
    vector<double> m_grt; //!< Standard chemical potentials for each species
};

}

#endif
