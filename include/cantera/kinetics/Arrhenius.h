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

}

#endif
