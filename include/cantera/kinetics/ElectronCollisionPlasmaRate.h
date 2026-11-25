//! @file ElectronCollisionPlasmaRate.h   Header for plasma reaction rates parameterized
//!     by electron collision cross section and electron energy distribution.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_ELECTRONCOLLISIONPLASMARATE_H
#define CT_ELECTRONCOLLISIONPLASMARATE_H

#include "cantera/kinetics/ReactionData.h"
#include "ReactionRate.h"
#include "MultiRate.h"

namespace Cantera
{

//! Data container holding shared data specific to ElectronCollisionPlasmaRate
/**
 * The data container `ElectronCollisionPlasmaData` holds precalculated data common to
 * all `ElectronCollisionPlasmaRate` objects.
 */
struct ElectronCollisionPlasmaData : public ReactionData
{
    ElectronCollisionPlasmaData();

    virtual bool update(const ThermoPhase& phase, const Kinetics& kin) override;
    using ReactionData::update;

    virtual void invalidateCache() override {
        ReactionData::invalidateCache();
        energyLevels.resize(0);
        distribution.resize(0);
        m_dist_number = -1;
    }

    vector<double> energyLevels; //!< electron energy levels
    vector<double> distribution; //!< electron energy distribution

    //! integer that is incremented when electron energy levels change
    int levelNumber = -1;

protected:
    //! integer that is incremented when electron energy distribution changes
    int m_dist_number = -1;
};


//! Electron collision plasma reaction rate type
/*!
 * The electron collision plasma reaction rate uses the electron collision
 * data and the electron energy distribution to calculate the reaction rate.
 * Hagelaar and Pitchford @cite hagelaar2005 define the reaction rate
 * coefficient (Eqn. 63) as,
 *
 *   @f[
 *        k =  \gamma \int_0^{\infty} \epsilon \sigma F_0 d\epsilon,
 *   @f]
 *
 * where @f$ \gamma = \sqrt{2/m_e} @f$ (Eqn.4 in @cite hagelaar2015),
 * @f$ m_e @f$ [kg] is the electron mass, @f$ \epsilon @f$ [J] is the electron
 * energy, @f$ \sigma @f$ [m²] is the reaction collision cross section,
 * @f$ F_0 @f$ [@f$ \mbox{J}^{-3/2} @f$] is the normalized electron energy
 * distribution function, and @f$ k @f$ has the unit of [m³/s]. The collision
 * process is treated as a bimolecular reaction and should have units of [m³/kmol/s].
 * Therefore the forward reaction coefficient becomes,
 *
 *   @f[
 *        k_f = \gamma N_A \int_0^{\infty} \epsilon \sigma F_0 d\epsilon,
 *   @f]
 *
 * where @f$ N_A @f$ [1/kmol] is the Avogadro's number. Since the unit of the
 * electron energy downloaded from LXCat is in [eV], the elementary charge, e,
 * can be taken out of the integral and combine with @f$ \gamma @f$ to get,
 *
 *   @f[
 *        k_f = \sqrt{\frac{2e}{m_e}} N_A \int_0^{\infty} \epsilon_V \sigma F_{0,V} d\epsilon_V.
 *   @f]
 *
 * In addition to the forward reaction coefficient, the reverse (super-elastic) reaction
 * coefficient can be written as,
 *
 *   @f[
 *        k_r = \sqrt{\frac{2e}{m_e}} N_A \int_0^{\infty} \epsilon_V \sigma_{super}
 *              F_{0,V} d\epsilon_V,
 *   @f]
 *
 * where @f$ \sigma_{super} @f$ is the super-elastic cross section, defined by the principle
 * of detailed balancing as:
 *
 *   @f[
 *        \sigma_{super}(\epsilon) = \frac{\epsilon + U}{\epsilon} \sigma(\epsilon + U),
 *   @f]
 *
 * where @f$ U @f$ is the collision threshold energy [eV], as returned
 * by #threshold().
 * 
 * Electron collision cross sections may be specified either as kinetic
 * reactions or as simple collision entries, and are in any case located
 * in the root-level `electron-collisions` section.
 * Entries without an equation contribute to the plasma collision model
 * without being registered as kinetic reactions.
 *
 * Reaction rate coefficients are integrated over a merged energy grid
 * containing the EEDF grid points and all cross-section breakpoints
 * within the EEDF domain. Piecewise-linear interpolation makes the
 * integrand cubic on each interval, where Simpson's rule is exact.
 *
 * Cross sections are held constant outside their tabulated range,
 * unless the phase enables the optional high-energy extrapolation.
 *
 * @ingroup otherRateGroup
 * @since New in %Cantera 3.1.
 */
class ElectronCollisionPlasmaRate : public ReactionRate
{
public:
    ElectronCollisionPlasmaRate() = default;
    //! Constructor from YAML input for ElectronCollisionPlasmaRate.
    /*!
    * This constructor is used to initialize an electron collision plasma rate
    * from an input YAML file. It extracts the energy levels, cross-sections,
    * and reaction metadata used in the rate coefficient calculation.
    *
    * @param node         The AnyMap node containing rate fields from YAML
    * @param rate_units   Units used for interpreting the rate fields
    */
    ElectronCollisionPlasmaRate(const AnyMap& node,
                                const UnitStack& rate_units={})
    {
        setParameters(node, rate_units);
    }

    virtual void setParameters(const AnyMap& node, const UnitStack& units) override;

    virtual void getParameters(AnyMap& node) const override;

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return unique_ptr<MultiRateBase>(
            new MultiRate<ElectronCollisionPlasmaRate, ElectronCollisionPlasmaData>);
    }

    virtual const std::string type() const override {
        return "electron-collision-plasma";
    }

    virtual void setContext(const Reaction& rxn, const Kinetics& kin) override;

    //! Evaluate reaction rate
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    double evalFromStruct(const ElectronCollisionPlasmaData& shared_data);

    //! Calculate the reverse rate coefficient for super-elastic collisions
    //! @param shared_data Data structure with energy levels and EEDF
    //! @param kf Forward rate coefficient (input, unused)
    //! @param kr Reverse rate coefficient (output, modified)
    void modifyRateConstants(const ElectronCollisionPlasmaData& shared_data,
                             double& kf, double& kr);

    //! Evaluate derivative of reaction rate with respect to temperature
    //! divided by reaction rate
    //! @param shared_data  data shared by all reactions of a given type
    double ddTScaledFromStruct(const ElectronCollisionPlasmaData& shared_data) const {
        throw NotImplementedError("ElectronCollisionPlasmaRate::ddTScaledFromStruct");
    }

    //! The kind of the process which will be one of the following:
    //! - `"effective"`: Effective momentum-transfer cross section containing
    //!   elastic and inelastic contributions
    //! - `"elastic"`: Elastic momentum-transfer collision
    //! - `"excitation"`: Electronic, vibrational or rotational excitation
    //! - `"ionization"`: Electron-impact ionization
    //! - `"attachment"`: Electron attachment
    //! @since New in Cantera 3.2.
    const string& kind() const {
        return m_kind;
    }

    //! Get the target species of the electron collision process.
    //! This is the name of the neutral or ionic species that the electron interacts with
    //! @since New in Cantera 3.2.
    const string& target() const {
        return m_target;
    }

    //! Get the product of the electron collision process.
    //! This is the name of the species or excited state of
    //! some species resulting from the process.
    //! @note this may not necessarily be represented by
    //! a distinct species in the mixture.
    //! @since New in Cantera 3.2.
    const string& product() const {
        return m_product;
    }

    //! Get the electron collision threshold [eV].
    //!
    //! If explicitly specified, the threshold is preserved, including if
    //! presenting of value of zero.
    //!
    //! Otherwise, for excitation, ionization, and attachment collisions,
    //! the threshold is inferred from the first positive cross-section
    //! value. If that value is preceded by a zero cross section, the energy
    //! of the preceding zero is used; otherwise, the energy of the first
    //! positive value is used.
    //!
    //! Excitation and ionization cross sections are normalized to have
    //! a zero value at the threshold, following the convention in BOLSIG +.
    //!
    //! @note This behavior may be subject to change.
    //!
    //! @since New in Cantera 3.2. Modified in 4.0.
    double threshold() const {
        return m_threshold;
    }

    //! The value of #m_energyLevels [eV]
    span<const double> energyLevels() const {
        return m_energyLevels;
    }

    //! The value of #m_crossSections [m2]
    span<const double> crossSections() const {
        return m_crossSections;
    }

    //! Return the local cross-section extrapolation preference.
    //! @since New in %Cantera 4.0.
    bool extrapolateCrossSectionsRequested() const {
        return m_extrapolateCrossSections;
    }

    //! Enable or disable high-energy cross-section extrapolation.
    //! The phase determines the effective setting.
    //! @since New in %Cantera 4.0.
    void setCrossSectionExtrapolation(bool enabled);

    //! The value of #m_crossSectionsInterpolated [m2]
    span<const double> crossSectionInterpolated() const {
        return m_crossSectionsInterpolated;
    }

    //! Update the value of #m_crossSectionsInterpolated [m2]
    void updateInterpolatedCrossSection(span<const double>);

    //! Return whether tabulated cross-section data have been assigned.
    //! @since New in %Cantera 4.0.
    bool hasCrossSectionData() const {
        return m_hasCrossSectionData;
    }

    //! Infer collision metadata directly from the reaction equation.
    //! Does not require a Kinetics object or a fully initialized phase.
    //! @since New in %Cantera 4.0.
    void setReactionContext(const Reaction& rxn, const string& electronName);

private:
    //! The name of the kind of electron collision
    string m_kind;

    //! The name of the target of electron collision
    string m_target;

    //! The product of electron collision
    string m_product;

    //! Original note attached to the collision cross section
    string m_note;

    //! The energy threshold of electron collision
    double m_threshold = 0.0;

    //! electron energy levels [eV]
    vector<double> m_energyLevels;

    //! Counter used to indicate when #m_energyLevels needs to be synced with the phase
    int m_levelNumber = -3;

    //! collision cross sections [m2] at #m_energyLevels
    vector<double> m_crossSections;

    //! collision cross sections [m2] after interpolation
    vector<double> m_crossSectionsInterpolated;

    //! Contains all EEDF energy levels and all cross-section tabulation
    //! points within the EEDF energy range.
    vector<double> m_rateIntegrationEnergy;

    //! Cross sections evaluated on the forward integration grid.
    vector<double> m_rateIntegrationCrossSection;

    //! EEDF grid revision used to build the forward integration cache.
    int m_rateIntegrationLevelNumber = -4;

    //! Energy grid used for super-elastic rate integration.
    //! Cross-section breakpoints are shifted by the threshold.
    vector<double> m_superElasticRateIntegrationEnergy;

    //! Excitation cross sections evaluated on the reverse grid.
    vector<double> m_superElasticRateIntegrationCrossSection;

    //! EEDF grid revision used to build the reverse integration cache.
    int m_superElasticRateIntegrationLevelNumber = -5;

    //! Whether this rate contains validated tabulated cross-section data.
    bool m_hasCrossSectionData = false;

    //! True when threshold was explicitly provided in the input.
    bool m_hasExplicitThreshold = false;

    //! True when the target was determined from a reaction equation.
    bool m_isReactionRate = false;

    bool m_extrapolateCrossSections = false;
    size_t m_tabulatedPointCount = 0;
    bool m_crossSectionsExtrapolated = false;

    //! Validate the contents of a named electron-collision definition.
    void validateCollisionData(const AnyMap& node) const;

    //! Infer the collision threshold when it is not specified explicitly.
    void setDefaultThreshold();

    //! Normalize excitation and ionization cross sections at the threshold.
    void enforceThresholdCrossSections();

    //! Extend tabulated cross sections using a high-energy logarithmic tail.
    void extrapolateCrossSections();

    //! Invalidate rate-integration grids after cross-section changes.
    void invalidateRateIntegrationCaches();
};

}

#endif
