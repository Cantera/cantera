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
};


//! Electron collision plasma reaction rate type
/*!
 * The electron collision plasma reaction rate uses the electron collision
 * data and the electron energy distribution to calculate the reaction rate.
 *
 *   @f[
 *        k_f =  \gamma \int_0^{\infty} \epsilon \sigma F_0 d\epsilon
 *   @f]
 *
 * where @f$ \gamma = (2e/m)^{1/2} @f$ is a constant, @f$ \epsilon @f$ is
 * the electron energy, @f$ \sigma @f$ is the reaction collision cross
 * section, and @f$ F_0 @f$ is the normalized electron energy distribution
 * function.
 *
 * References:
 *
 * [1] Hagelaar and Pitchford @cite hagelaar2005.
 *
 * @ingroup plasma
 */
class ElectronCollisionPlasmaRate : public ReactionRate
{
public:
    ElectronCollisionPlasmaRate() = default;

    ElectronCollisionPlasmaRate(const AnyMap& node,
                                const UnitStack& rate_units={})
        : ElectronCollisionPlasmaRate()
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
    double evalFromStruct(const ElectronCollisionPlasmaData& shared_data) const;

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
};

}

#endif
