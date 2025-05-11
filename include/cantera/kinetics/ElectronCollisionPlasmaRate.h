//! @file ElectronCollisionPlasmaRate.h   Header for plasma reaction rates parameterized
//!     by electron collision cross section and electron energy distribution.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_ELECTRONCOLLISIONPLASMARATE_H
#define CT_ELECTRONCOLLISIONPLASMARATE_H

#include "cantera/thermo/PlasmaPhase.h"
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
    bool levelChanged;

protected:
    //! integer that is incremented when electron energy distribution changes
    int m_dist_number = -1;

    //! integer that is incremented when electron energy level changes
    int m_level_number = -1;
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
 * where @f$ U @f$ is the threshold energy [eV], equal to the first energy level of the cross
 * section (#m_energyLevels).
 *
 * @ingroup otherRateGroup
 * @since New in %Cantera 3.1.
 */
class ElectronCollisionPlasmaRate : public ReactionRate
{
public:
    ElectronCollisionPlasmaRate() = default;

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

    void modifyRateConstants(const ElectronCollisionPlasmaData& shared_data,
                             double& kf, double& kr);

    //! Evaluate derivative of reaction rate with respect to temperature
    //! divided by reaction rate
    //! @param shared_data  data shared by all reactions of a given type
    double ddTScaledFromStruct(const ElectronCollisionPlasmaData& shared_data) const {
        throw NotImplementedError("ElectronCollisionPlasmaRate::ddTScaledFromStruct");
    }

    //! The kind of the process
    const string& kind() const {
        return m_kind;
    }

    //! The target of the process
    const string& target() const {
        return m_target;
    }

    //! The product of the process
    const string& product() const {
        return m_product;
    }


    //! Set the value of #m_energyLevels [eV]
    void set_energyLevels(vector<double> epsilon) {
        m_energyLevels = epsilon;
    }

    //! Set the value of #m_crossSections [eV]
    void set_crossSections(vector<double> sigma) {
        m_crossSections = sigma;
    }

    //! Set the value of m_threshold [eV]
    void set_threshold(double threshold)
    {
        m_threshold = threshold;
    }

    void set_cs_ok() {
        cs_ok = true;
    }

    const bool get_cs_ok() const {
        return cs_ok;
    }


    //! The value of #m_energyLevels [eV]
    const vector<double>& energyLevels() const {
        return m_energyLevels;
    }

    //! The value of #m_crossSections [m2]
    const vector<double>& crossSections() const {
        return m_crossSections;
    }

    //! The value of #m_crossSectionsInterpolated [m2]
    const vector<double>& crossSectionInterpolated() const {
        return m_crossSectionsInterpolated;
    }

    //! Set the value of #m_crossSectionsInterpolated
    void setCrossSectionInterpolated(vector<double>& cs) {
        m_crossSectionsInterpolated = cs;
    }

    //! Update the value of #m_crossSectionsInterpolated [m2]
    void updateInterpolatedCrossSection(const vector<double>&);

private:
    //! The name of the kind of electron collision
    string m_kind;

    //! The name of the target of electron collision
    string m_target;

    //! The product of electron collision
    string m_product;

    //! The energy threshold of electron collision
    double m_threshold;

    //! Check if a cross-section is define for this rate
    bool cs_ok = false;

    //! electron energy levels [eV]
    vector<double> m_energyLevels;

    //! collision cross sections [m2] at #m_energyLevels
    vector<double> m_crossSections;

    //! collision cross sections [m2] after interpolation
    vector<double> m_crossSectionsInterpolated;

    //! collision cross section [m2] interpolated on #m_energyLevels offset by the
    //! threshold energy (the first energy level).
    //! This is used for the calculation of the super-elastic collision reaction
    //! rate coefficient.
    Eigen::ArrayXd m_crossSectionsOffset;
};

}

#endif
