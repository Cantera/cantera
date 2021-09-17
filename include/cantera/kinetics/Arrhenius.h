// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_ARRHENIUS_H
#define CT_ARRHENIUS_H

#include "cantera/base/ct_defs.h"
#include "cantera/kinetics/ReactionData.h"

namespace Cantera
{

class AnyValue;
class AnyMap;
class Units;
class UnitSystem;

/**
 *  @defgroup arrheniusGroup  Arrhenius-type Parameterizations
 *
 *  This section describes the parameterizations used to describe the standard
 *  Arrhenius rate parameterization and derived models.
 *
 *  @ingroup chemkinetics
 */

//! Base class for Arrhenius-type Parameterizations
class ArrheniusBase
{
public:
    //! Default constructor.
    ArrheniusBase();

    //! Constructor.
    /*!
     *  @param A  Pre-exponential factor. The unit system is (kmol, m, s); actual units
     *      depend on the reaction order and the dimensionality (surface or bulk).
     *  @param b  Temperature exponent (non-dimensional)
     *  @param Ea  Activation energy in energy units [J/kmol]
     */
    ArrheniusBase(double A, double b, double Ea);

    //! Constructor based on AnyValue content
    ArrheniusBase(const AnyValue& rate,
                  const UnitSystem& units, const Units& rate_units)
    {
        setParameters(rate, units, rate_units);
    }

    //! Perform object setup based on AnyValue node information
    /*!
     *  @param node  AnyValue containing rate information
     *  @param units  Unit system
     *  @param rate_units  Unit definitions specific to rate information
     */
    virtual void setParameters(const AnyValue& rate,
                               const UnitSystem& units, const Units& rate_units);

    //! Get parameters
    //! Store the parameters of a ReactionRate needed to reconstruct an identical
    virtual void getParameters(AnyMap& node, const Units& units) const;

    //! Validate the reaction rate expression
    virtual void validate(const std::string& equation);

    //! Return the pre-exponential factor *A* (in m, kmol, s to powers depending
    //! on the reaction order)
    double preExponentialFactor() const {
        return m_A;
    }

    //! Return the temperature exponent *b*
    double temperatureExponent() const {
        return m_b;
    }

    bool allow_negative_pre_exponential_factor; // Flag is directly accessible

    size_t rate_index; //!< Reaction rate index within kinetics evaluator

protected:
    double m_A; //!< Pre-exponential factor
    double m_b; //!< Temperature exponent
    double m_Ea_R; //!< Activation energy (in temperature units)
};


//! Arrhenius reaction rate type depends only on temperature
/*!
 * A reaction rate coefficient of the following form.
 *
 *   \f[
 *        k_f =  A T^b \exp (-Ea/RT)
 *   \f]
 *
 * @todo supersedes Arrhenius: replace existing instances with this class, rename,
 *      and deprecate Arrhenius3.
 *
 * @ingroup arrheniusGroup
 */
class Arrhenius3 : public ArrheniusBase
{
public:
    using ArrheniusBase::ArrheniusBase; // inherit constructors

    //! Constructor based on AnyMap content
    Arrhenius3(const AnyMap& node, const Units& rate_units) {
        setParameters(node, rate_units);
    }

    //! Identifier of reaction rate type
    // const static std::string type()
    virtual std::string type() const {
        return "Arrhenius";
    }

    //! Perform object setup based on AnyMap node information
    /*!
     *  @param node  AnyMap containing rate information
     *  @param rate_units  Unit definitions specific to rate information
     */
    virtual void setParameters(const AnyMap& node, const Units& rate_units);

    virtual void getParameters(AnyMap& node, const Units& units) const;

    //! Update information specific to reaction
    const static bool usesUpdate() {
        return false;
    }

    //! Update information specific to reaction
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    void update(const ArrheniusData& shared_data) {}

    //! Evaluate reaction rate
    //! @param shared_data  data shared by all reactions of a given type
    double eval(const ArrheniusData& shared_data) const {
        return m_A * std::exp(m_b * shared_data.logT - m_Ea_R * shared_data.recipT);
    }

    //! Evaluate derivative of reaction rate with respect to temperature
    //! divided by reaction rate value
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    virtual double ddTscaled(const ArrheniusData& shared_data) const {
        return (m_b + m_Ea_R * shared_data.recipT) * shared_data.recipT;
    }

    //! Return the activation energy *Ea* [J/kmol]
    double activationEnergy() const {
        return m_Ea_R * GasConstant;
    }

    //! Return the activation energy divided by the gas constant (i.e. the
    //! activation temperature) [K]
    double activationEnergy_R() const {
        return m_Ea_R;
    }
};

}

#endif
