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
 * energy, @f$ \sigma @f$ [m2] is the reaction collision cross section,
 * @f$ F_0 @f$ [J^(-3/2)] is the normalized electron energy distribution function,
 * and @f$ k @f$ has the unit of [m3/s]. The collision process is treated as a
 * bimolecular reaction and should have units of [m3/kmol/s]. Therefore the
 * forward reaction coefficient becomes,
 *
 *   @f[
 *        k_f = \gamma N_A \int_0^{\infty} \epsilon \sigma F_0 d\epsilon,
 *   @f]
 *
 * where @f$ N_A @f$ [1/kmol] is the Avogadro's number. Since the unit of the
 * electron energy downloaded from LXCat is in [V], the forward reaction
 * coefficient can be written as,
 *
 *   @f[
 *        k_f = \sqrt{\frac{2e}{m_e}} N_A \int_0^{\infty} \epsilon_V \sigma F_0 d\epsilon_V.
 *   @f]
 *
 * For the convenience of calculation, the final form becomes,
 *
 *   @f[
 *        k_f = 0.5 \sqrt{\frac{2e}{m_e}} N_A \int_0^{\infty} \sigma F_0 d{{\epsilon_V}^2}.
 *   @f]
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

    //! Evaluate derivative of reaction rate with respect to temperature
    //! divided by reaction rate
    //! @param shared_data  data shared by all reactions of a given type
    double ddTScaledFromStruct(const ElectronCollisionPlasmaData& shared_data) const {
        throw NotImplementedError("ElectronCollisionPlasmaRate::ddTScaledFromStruct");
    }

    //! The value of #m_energyLevels [eV]
    const vector<double>& energyLevels() const {
        return m_energyLevels;
    }

    //! The value of #m_crossSections [m2]
    const vector<double>& crossSections() const {
        return m_crossSections;
    }

private:
    //! electron energy levels [eV]
    vector<double> m_energyLevels;

    //! collision cross sections [m2] at #m_energyLevels
    vector<double> m_crossSections;

    //! collision cross sections [m2] after interpolation
    vector<double> m_crossSectionsInterpolated;
};

}

#endif
