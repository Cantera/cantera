//! @file TwoTempPlasmaRate.h   Header for plasma reaction rates parameterized by two
//!     temperatures (gas and electron).

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_TWOTEMPPLASMARATE_H
#define CT_TWOTEMPPLASMARATE_H

#include "Arrhenius.h"

namespace Cantera
{

//! Data container holding shared data specific to TwoTempPlasmaRate
/**
 * The data container `TwoTempPlasmaData` holds precalculated data common to
 * all `TwoTempPlasmaRate` objects.
 */
struct TwoTempPlasmaData : public ReactionData
{
    TwoTempPlasmaData() = default;

    bool update(const ThermoPhase& phase, const Kinetics& kin) override;
    void update(double T) override;
    void update(double T, double Te) override;
    using ReactionData::update;

    virtual void updateTe(double Te);

    void invalidateCache() override {
        ReactionData::invalidateCache();
        electronTemp = NAN;
    }

    double electronTemp = 1.0; //!< electron temperature
    double logTe = 0.0; //!< logarithm of electron temperature
    double recipTe = 1.0; //!< inverse of electron temperature
};


//! Two temperature plasma reaction rate type depends on both
//! gas temperature and electron temperature.
/*!
 * The form of the two temperature plasma reaction rate coefficient is similar to an
 * Arrhenius reaction rate coefficient. The temperature exponent (b) is applied to
 * the electron temperature instead. In addition, the exponential term with
 * activation energy for electron is included.
 *
 *   @f[
 *        k_f =  A T_e^b \exp (-E_{a,g}/RT) \exp (E_{a,e} (T_e - T)/(R T T_e))
 *   @f]
 *
 * where @f$ T_e @f$ is the electron temperature, @f$ E_{a,g} @f$ is the activation
 * energy for gas, and @f$ E_{a,e} @f$ is the activation energy for electron, see
 * Kossyi, et al. @cite kossyi1992.
 *
 * @ingroup arrheniusGroup
 */
class TwoTempPlasmaRate : public ArrheniusBase
{
public:
    TwoTempPlasmaRate();

    //! Constructor.
    /*!
     *  @param A  Pre-exponential factor. The unit system is (kmol, m, s); actual units
     *      depend on the reaction order and the dimensionality (surface or bulk).
     *  @param b  Temperature exponent (non-dimensional)
     *  @param Ea  Activation energy in energy units [J/kmol]
     *  @param EE  Activation electron energy in energy units [J/kmol]
     */
    TwoTempPlasmaRate(double A, double b, double Ea=0.0, double EE=0.0);

    TwoTempPlasmaRate(const AnyMap& node, const UnitStack& rate_units={});

    unique_ptr<MultiRateBase> newMultiRate() const override {
        return make_unique<MultiRate<TwoTempPlasmaRate, TwoTempPlasmaData>>();
    }

    const string type() const override {
        return "two-temperature-plasma";
    }

    void setContext(const Reaction& rxn, const Kinetics& kin) override;

    //! Evaluate reaction rate
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    double evalFromStruct(const TwoTempPlasmaData& shared_data) const {
        // m_E4_R is the electron activation (in temperature units)
        return m_A * std::exp(m_b * shared_data.logTe -
                              m_Ea_R * shared_data.recipT +
                              m_E4_R * (shared_data.electronTemp - shared_data.temperature)
                              * shared_data.recipTe * shared_data.recipT);
    }

    //! Evaluate derivative of reaction rate with respect to temperature
    //! divided by reaction rate
    /*!
     *  This method does not consider changes of electron temperature.
     *  A corresponding warning is raised.
     *  @param shared_data  data shared by all reactions of a given type
     */
    double ddTScaledFromStruct(const TwoTempPlasmaData& shared_data) const;

    //! Return the electron activation energy *Ea* [J/kmol]
    double activationElectronEnergy() const {
        return m_E4_R * GasConstant;
    }
};

}

#endif
