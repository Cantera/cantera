/**
 * @file InterfaceRate.h
 * Header for reaction rates that occur at interfaces.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_INTERFACERATE_H
#define CT_INTERFACERATE_H

#include "cantera/base/global.h"
#include "cantera/kinetics/BlowersMaselRate.h"
#include "MultiRate.h"

namespace Cantera
{

class AnyMap;

/**
 *  @defgroup surfaceGroup  Coverage-dependent rate parameterizations
 *
 *  This section describes the parameterizations used to describe rate
 *  parameterization that involve interfaces.
 *
 *  @ingroup chemkinetics
 */


//! Data container holding shared data for reaction rate specification with interfaces
/**
 * The data container InterfaceData holds precalculated data common to
 * InterfaceRate and StickingRate objects.
 *
 * The data container inherits from BlowersMaselData, where density is used to
 * hold the site density [kmol/m^2].
 */
struct InterfaceData : public BlowersMaselData
{
    InterfaceData();

    virtual bool update(const ThermoPhase& bulk, const Kinetics& kin) override;

    virtual void update(double T) override;

    virtual void update(double T, const vector_fp& values) override;

    using BlowersMaselData::update;

    virtual void perturbTemperature(double deltaT);

    virtual void resize(size_t nSpecies, size_t nReactions, size_t nPhases) override {
        coverages.resize(nSpecies, 0.);
        logCoverages.resize(nSpecies, 0.);
        partialMolarEnthalpies.resize(nSpecies, 0.);
        electricPotentials.resize(nPhases, 0.);
        standardChemPotentials.resize(nSpecies, 0.);
        standardConcentrations.resize(nSpecies, 0.);
        ready = true;
    }

    double sqrtT; //!< square root of temperature

    vector_fp coverages; //!< surface coverages
    vector_fp logCoverages; //!< logarithm of surface coverages
    vector_fp electricPotentials; //!< electric potentials of phases
    vector_fp standardChemPotentials; //!< standard state chemical potentials
    vector_fp standardConcentrations; //!< standard state concentrations
};


//! Base class for rate parameterizations that involve interfaces
/**
 * Rate expressions defined for interfaces may include coverage dependent terms,
 * where an example is given by [Kee, R. J., Coltrin, M. E., & Glarborg, P.(2003).
 * Chemically reacting flow: theory and practice. John Wiley & Sons. Eq 11.113].
 * Using Cantera nomenclature, this expression can be rewritten as
 *  \f[
 *      k_f = A T^b \exp \left( - \frac{E_a}{RT} \right)
 *          \prod_k 10^{a_k \theta_k} \theta_k^{m_k}
 *          \exp \left( \frac{- E_k \theta_k}{RT} \right)
 *  \f]
 * It is evident that this expression combines a regular modified Arrhenius rate
 * expression \f$ A T^b \exp \left( - \frac{E_a}{RT} \right) \f$ with coverage-related
 * terms, where the parameters \f$ (a_k, E_k, m_k) \f$ describe the dependency on the
 * surface coverage of species \f$ k, \theta_k \f$. The InterfaceRateBase class implements
 * terms related to coverage only, which allows for combinations with arbitrary rate
 * parameterizations (for example Arrhenius and BlowersMaselRate).
 */
class InterfaceRateBase
{
public:
    InterfaceRateBase();

    //! Perform object setup based on AnyMap node information
    //! @param node  AnyMap object containing reaction rate specification
    void setParameters(const AnyMap& node);

    //! Store parameters needed to reconstruct an identical object
    //! @param node  AnyMap object receiving reaction rate specification
    void getParameters(AnyMap& node) const;

    //! Set coverage dependencies based on AnyMap node information
    //! @param dependencies  Coverage dependencies
    //! @param units  Unit system
    void setCoverageDependencies(const AnyMap& dependencies,
                                 const UnitSystem& units=UnitSystem());

    //! Store parameters needed to reconstruct coverage dependencies
    //! @param dependencies  AnyMap receiving coverage information
    //! @param asVector  Optional boolean flag to override map output
    //! @todo  Remove vector version (which currently only serves testing purposes)
    void getCoverageDependencies(AnyMap& dependencies, bool asVector=false) const;

    //! Add a coverage dependency for species *sp*, with exponential dependence
    //! *a*, power-law exponent *m*, and activation energy dependence *e*,
    //! where *e* is in Kelvin, that is, energy divided by the molar gas constant.
    void addCoverageDependence(const std::string& sp, double a, double m, double e);

    //! Boolean indicating whether rate uses exchange current density formulation
    bool exchangeCurrentDensityFormulation() {
        return m_exchangeCurrentDensityFormulation;
    }

    //! Build rate-specific parameters based on Reaction and Kinetics context
    //! @param rxn  Reaction associated with rate parameterization
    //! @param kin  Kinetics object associated with rate parameterization
    void setContext(const Reaction& rxn, const Kinetics& kin);

    //! Set association with an ordered list of all species associated with a given
    //! `Kinetics` object.
    void setSpecies(const std::vector<std::string>& species);

    //! Update reaction rate parameters
    //! @param shared_data  data shared by all reactions of a given type
    void updateFromStruct(const InterfaceData& shared_data);

    //! Calculate modifications for the forward reaction rate for interfacial charge
    //! transfer reactions.
    /*!
     *  For reactions that transfer charge across a potential difference, the
     *  activation energies are modified by the potential difference. The correction
     *  factor is based on the net electric potential energy change
     *  \f[
     *   \Delta E_{p,j} = \sum_i E_{p,i} \nu_{i,j}
     *  \f]
     *  where potential energies are calculated as \f$ E_{p,i} = F \phi_i z_i \f$.
     *  Here, \f$ F \f$ is Faraday's constant, \f$ \phi_i \f$ is the electric potential
     *  of the species phase and \f$ z_i \f$ is the charge of the species.
     *
     *  When an electrode reaction rate is specified in terms of its exchange current
     *  density, the correction factor is adjusted to the standard reaction rate
     *  constant form and units. Specifically, this converts a reaction rate constant
     *  that was specified in units of A/m2 to kmol/m2/s.
     *
     *  @warning  The updated calculation of voltage corrections is an experimental
     *      part of the %Cantera API and may be changed or removed without notice.
     */
    double voltageCorrection() const {
        // Calculate reaction rate correction. Only modify those with a non-zero
        // activation energy.
        double correction = 1.;
        if (m_deltaPotential_RT != 0.) {
            // Comments preserved from previous implementation:
            // Below we decrease the activation energy below zero.
            // NOTE, there is some discussion about this point. Should we decrease the
            // activation energy below zero? I don't think this has been decided in any
            // definitive way. The treatment below is numerically more stable, however.
            correction = exp(-m_beta * m_deltaPotential_RT);
        }

        // Update correction if exchange current density formulation format is used.
        if (m_exchangeCurrentDensityFormulation) {
            // Comment preserved from previous implementation:
            // We need to have the straight chemical reaction rate constant to
            // come out of this calculation.
            double tmp = exp(-m_beta * m_deltaGibbs0_RT);
            tmp /= m_prodStandardConcentrations * Faraday;
            correction *= tmp;
        }
        return correction;
    }

    //! Boolean indicating whether rate uses electrochemistry
    /*!
     *  If this is true, the Butler-Volmer correction
     *  \f[
     *    f_{BV} = \exp ( - \beta * Delta E_{p,j} / R T )
     *  \f]
     *  is applied to the forward reaction rate, @see voltageCorrection.
     */
    bool usesElectrochemistry() {
        return m_chargeTransfer;
    }

    //! Return the charge transfer beta parameter
    double beta() const {
        if (m_chargeTransfer) {
            return m_beta;
        }
        return NAN;
    }

    //! Return site density [kmol/m^2]
    /*!
     *  @warning  This method is an experimental part of the %Cantera API and
     *      may be changed or removed without notice.
     */
    double siteDensity() const {
        return m_siteDensity;
    }

    //! Set site density [kmol/m^2]
    /*!
     *  @internal  This method is used for testing purposes only as the site density
     *      is a property of InterfaceKinetics and will be overwritten during an update
     *      of the thermodynamic state.
     *
     *  @warning  This method is an experimental part of the %Cantera API and
     *      may be changed or removed without notice.
     */
    void setSiteDensity(double siteDensity) {
        m_siteDensity = siteDensity;
    }

protected:
    double m_siteDensity; //!< Site density [kmol/m^2]
    double m_acov; //!< Coverage contribution to pre-exponential factor
    double m_ecov; //!< Coverage contribution to activation energy
    double m_mcov; //!< Coverage term in reaction rate
    bool m_chargeTransfer; //!< Boolean indicating use of electrochemistry
    bool m_exchangeCurrentDensityFormulation; //! Electrochemistry only
    double m_beta; //!< Forward value of apparent electrochemical transfer coefficient
    double m_deltaPotential_RT; //!< Normalized electric potential energy change
    double m_deltaGibbs0_RT; //!< Normalized standard state Gibbs free energy change
    double m_prodStandardConcentrations; //!< Products of standard concentrations
    std::map<size_t, size_t> m_indices; //!< Map holding indices of coverage species
    std::vector<std::string> m_cov; //!< Vector holding names of coverage species
    vector_fp m_ac; //!< Vector holding coverage-specific exponential dependence
    vector_fp m_ec; //!< Vector holding coverage-specific activation energy dependence
    vector_fp m_mc; //!< Vector holding coverage-specific power-law exponents

private:
    //! Pairs of species index and multiplers to calculate enthalpy change
    std::vector<std::pair<size_t, double>> m_stoichCoeffs;

    //! Pairs of phase index and net electric charges (same order as m_stoichCoeffs)
    std::vector<std::pair<size_t, double>> m_netCharges;
};


//! Base class for rate parameterizations that implement sticking coefficients
/**
 * The StickingCoverage class enhances Coverage to accommodate sticking coefficients.
 */
class StickingCoverage : public InterfaceRateBase
{
public:
    StickingCoverage();

    //! Perform object setup based on AnyMap node information
    //! @param node  Sticking coefficient parameters
    void setStickingParameters(const AnyMap& node);

    //! Store parameters needed to reconstruct an identical object
    //! @param node  Sticking coefficient parameters
    void getStickingParameters(AnyMap& node) const;

    //! Get flag indicating whether sticking rate uses the correction factor developed
    //! by Motz & Wise for reactions with high (near-unity) sticking coefficients.
    //! Defaults to 'false'.
    bool motzWiseCorrection() const {
        return m_motzWise;
    }

    //! Set flag for Motz & Wise correction factor
    void setMotzWiseCorrection(bool motz_wise) {
        m_motzWise = motz_wise;
        m_explicitMotzWise = true;
    }

    //! Get sticking species.
    std::string stickingSpecies() const {
        return m_stickingSpecies;
    }

    //! Set sticking species
    /*!
     *  For reactions with multiple non-surface species, the sticking species needs
     *  to be explicitly identified. Note that species have to be specified prior
     *  to adding a reaction to a Kinetics object.
     */
    void setStickingSpecies(const std::string& stickingSpecies) {
        m_stickingSpecies = stickingSpecies;
        m_explicitSpecies = true;
    }

    //! Get exponent applied to site density (sticking order)
    double stickingOrder() {
        return m_surfaceOrder;
    }

    //! Set exponent applied to site density (sticking order)
    /*!
     *  @internal  This method is used for testing purposes only as the value is
     *      determined automatically by setContext.
     *
     *  @warning  This method is an experimental part of the %Cantera API and
     *      may be changed or removed without notice.
     */
    void setStickingOrder(double order) {
        m_surfaceOrder = order;
    }

    //! Get the molecular weight of the sticking species
    double stickingWeight() {
        return GasConstant / (2 * Pi * m_multiplier * m_multiplier);
    }

    //! Set the molecular weight of the sticking species
    /*!
     *  @internal  This method is used for testing purposes only as the value is
     *      determined automatically by setContext.
     *
     *  @warning  This method is an experimental part of the %Cantera API and
     *      may be changed or removed without notice.
     */
    void setStickingWeight(double weight) {
        m_multiplier = sqrt(GasConstant / (2 * Pi * weight));
    }

    //! Build rate-specific parameters based on Reaction and Kinetics context
    //! @param rxn  Reaction associated with the sticking coefficient
    //! @param kin  Kinetics object associated with the sticking coefficient
    //! Parameters can be accessed using the method stickingSpecies, stickingOrder
    //! and stickingWeight.
    void setContext(const Reaction& rxn, const Kinetics& kin);

protected:
    bool m_motzWise; //!< boolean indicating whether Motz & Wise correction is used
    bool m_explicitMotzWise; //!< Correction cannot be overriden by default
    std::string m_stickingSpecies; //!< string identifying sticking species
    bool m_explicitSpecies; //!< Boolean flag
    double m_surfaceOrder; //!< exponent applied to site density term
    double m_multiplier; //!< multiplicative factor in rate expression
    double m_factor; //!< cached factor
};


//! A class template for interface reaction rate specifications
template <class RateType, class DataType>
class InterfaceRate : public RateType, public InterfaceRateBase
{
    CT_DEFINE_HAS_MEMBER(has_update, updateFromStruct)

public:
    InterfaceRate() = default;
    using RateType::RateType; // inherit constructors

    //! Constructor based on AnyMap content
    InterfaceRate(const AnyMap& node, const UnitStack& rate_units) {
        setParameters(node, rate_units);
    }
    explicit InterfaceRate(const AnyMap& node) {
        setParameters(node, {});
    }

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return unique_ptr<MultiRateBase>(
            new MultiRate<InterfaceRate<RateType, DataType>, DataType>);
    }

    //! Identifier of reaction rate type
    virtual const std::string type() const override {
        return "interface-" + RateType::type();
    }

    virtual void setParameters(
        const AnyMap& node, const UnitStack& rate_units) override
    {
        InterfaceRateBase::setParameters(node);
        RateType::setParameters(node, rate_units);
    }

    virtual void getParameters(AnyMap& node) const override {
        RateType::getParameters(node);
        node["type"] = type();
        InterfaceRateBase::getParameters(node);
    }

    virtual void setContext(const Reaction& rxn, const Kinetics& kin) override {
        RateType::setContext(rxn, kin);
        InterfaceRateBase::setContext(rxn, kin);
    }

    //! Update reaction rate parameters
    //! @param shared_data  data shared by all reactions of a given type
    void updateFromStruct(const DataType& shared_data) {
        _update(shared_data);
        InterfaceRateBase::updateFromStruct(shared_data);
    }

    //! Evaluate reaction rate
    //! @param shared_data  data shared by all reactions of a given type
    double evalFromStruct(const DataType& shared_data) const {
        double out = RateType::evalRate(shared_data.logT, shared_data.recipT) *
            std::exp(std::log(10.0) * m_acov - m_ecov * shared_data.recipT + m_mcov);
        if (m_chargeTransfer) {
            out *= voltageCorrection();
        }
        return out;
    }

    //! Evaluate derivative of reaction rate with respect to temperature
    //! divided by reaction rate
    //! @param shared_data  data shared by all reactions of a given type
    double ddTScaledFromStruct(const DataType& shared_data) const {
        throw NotImplementedError("InterfaceRate<>::ddTScaledFromStruct");
    }

    virtual double preExponentialFactor() const override {
        return RateType::preExponentialFactor() *
            std::exp(std::log(10.0) * m_acov + m_mcov);
    }

    virtual double activationEnergy() const override {
        return RateType::activationEnergy() + m_ecov * GasConstant;
    }

protected:
    //! Helper function to process updates for rate types that implement the
    //! `updateFromStruct` method.
    template <typename T=RateType,
        typename std::enable_if<has_update<T>::value, bool>::type = true>
    void _update(const DataType& shared_data) {
        T::updateFromStruct(shared_data);
    }

    //! Helper function for rate types that do not implement `updateFromStruct`.
    //! Does nothing, but exists to allow generic implementations of update().
    template <typename T=RateType,
        typename std::enable_if<!has_update<T>::value, bool>::type = true>
    void _update(const DataType& shared_data) {
    }
};

using InterfaceArrheniusRate = InterfaceRate<ArrheniusRate, InterfaceData>;
using InterfaceBlowersMaselRate = InterfaceRate<BlowersMaselRate, InterfaceData>;


//! A class template for interface sticking rate specifications
template <class RateType, class DataType>
class StickingRate : public RateType, public StickingCoverage
{
    CT_DEFINE_HAS_MEMBER(has_update, updateFromStruct)

public:
    StickingRate() = default;
    using RateType::RateType; // inherit constructors

    //! Constructor based on AnyMap content
    StickingRate(const AnyMap& node, const UnitStack& rate_units) {
        // sticking coefficients are dimensionless
        setParameters(node, Units(1.0));
    }
    explicit StickingRate(const AnyMap& node) {
        // sticking coefficients are dimensionless
        setParameters(node, Units(1.0));
    }

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return unique_ptr<MultiRateBase>(
            new MultiRate<StickingRate<RateType, DataType>, DataType>);
    }

    //! Identifier of reaction rate type
    virtual const std::string type() const override {
        return "sticking-" + RateType::type();
    }

    virtual void setParameters(
        const AnyMap& node, const UnitStack& rate_units) override
    {
        InterfaceRateBase::setParameters(node);
        RateType::m_negativeA_ok = node.getBool("negative-A", false);
        setStickingParameters(node);
        if (!node.hasKey("sticking-coefficient")) {
            RateType::setRateParameters(AnyValue(), node.units(), rate_units);
            return;
        }
        RateType::setRateParameters(
            node["sticking-coefficient"], node.units(), rate_units);
    }

    virtual void getParameters(AnyMap& node) const override {
        node["type"] = type();
        if (RateType::m_negativeA_ok) {
            node["negative-A"] = true;
        }
        AnyMap rateNode;
        RateType::getRateParameters(rateNode);
        getStickingParameters(node);
        if (!rateNode.empty()) {
            // RateType object is configured
            node["sticking-coefficient"] = std::move(rateNode);
        }
        InterfaceRateBase::getParameters(node);
    }

    virtual void setContext(const Reaction& rxn, const Kinetics& kin) override {
        RateType::setContext(rxn, kin);
        InterfaceRateBase::setContext(rxn, kin);
        StickingCoverage::setContext(rxn, kin);
    }

    virtual void validate(const std::string &equation, const Kinetics& kin) override {
        RateType::validate(equation, kin);
        fmt::memory_buffer err_reactions;
        double T[] = {200.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0};
        for (size_t i=0; i < 6; i++) {
            double k = RateType::evalRate(log(T[i]), 1 / T[i]);
            if (k > 1) {
                fmt_append(err_reactions,
                    "\n Sticking coefficient is greater than 1 for reaction '{}'\n"
                    " at T = {:.1f}\n", equation, T[i]);
            }
        }
        if (err_reactions.size()) {
            warn_user("StickingRate::validate", to_string(err_reactions));
        }
    }

    //! Update reaction rate parameters
    //! @param shared_data  data shared by all reactions of a given type
    void updateFromStruct(const DataType& shared_data) {
        _update(shared_data);
        InterfaceRateBase::updateFromStruct(shared_data);
        m_factor = pow(m_siteDensity, -m_surfaceOrder);
    }

    //! Evaluate reaction rate
    //! @param shared_data  data shared by all reactions of a given type
    double evalFromStruct(const DataType& shared_data) const {
        double out = RateType::evalRate(shared_data.logT, shared_data.recipT) *
            std::exp(std::log(10.0) * m_acov - m_ecov * shared_data.recipT + m_mcov);
        if (m_chargeTransfer) {
            // @todo  the physical interpretation of a 'sticking' charge transfer
            //      reaction remains to be resolved.
            out *= voltageCorrection();
        }
        if (m_motzWise) {
            out /= 1 - 0.5 * out;
        }
        return out * m_factor * shared_data.sqrtT * m_multiplier;
    }

    //! Evaluate derivative of reaction rate with respect to temperature
    //! divided by reaction rate
    //! @param shared_data  data shared by all reactions of a given type
    double ddTScaledFromStruct(const DataType& shared_data) const {
        throw NotImplementedError("StickingRate<>::ddTScaledFromStruct");
    }

    virtual double preExponentialFactor() const override {
        return RateType::preExponentialFactor() *
            std::exp(std::log(10.0) * m_acov + m_mcov);
    }

    virtual double activationEnergy() const override {
        return RateType::activationEnergy() + m_ecov * GasConstant;
    }

protected:
    //! Helper function to process updates for rate types that implement the
    //! `updateFromStruct` method.
    template <typename T=RateType,
        typename std::enable_if<has_update<T>::value, bool>::type = true>
    void _update(const DataType& shared_data) {
        T::updateFromStruct(shared_data);
    }

    //! Helper function for rate types that do not implement `updateFromStruct`.
    //! Does nothing, but exists to allow generic implementations of update().
    template <typename T=RateType,
        typename std::enable_if<!has_update<T>::value, bool>::type = true>
    void _update(const DataType& shared_data) {
    }
};

using StickingArrheniusRate = StickingRate<ArrheniusRate, InterfaceData>;
using StickingBlowersMaselRate = StickingRate<BlowersMaselRate, InterfaceData>;

}
#endif
