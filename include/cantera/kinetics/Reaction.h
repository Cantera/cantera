/**
 *  @file Reaction.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTION_H
#define CT_REACTION_H

#include "cantera/base/AnyMap.h"
#include "cantera/base/Units.h"
#include "Arrhenius.h"
#include "ChebyshevRate.h"
#include "Custom.h"
#include "Falloff.h"
#include "InterfaceRate.h"
#include "PlogRate.h"

namespace Cantera
{

class Kinetics;
class ThirdBody;

//! @defgroup reactionGroup Reactions and reaction rates

//! Abstract base class which stores data about a reaction and its rate
//! parameterization so that it can be added to a Kinetics object.
//! @ingroup reactionGroup
class Reaction
{
public:
    Reaction() {}
    Reaction(const Composition& reactants, const Composition& products,
             shared_ptr<ReactionRate> rate={});

    //! Construct a Reaction and it's corresponding ReactionRate based on AnyMap (YAML)
    //! input.
    Reaction(const AnyMap& node, const Kinetics& kin);

    virtual ~Reaction() {}

    //! The reactant side of the chemical equation for this reaction
    std::string reactantString() const;

    //! The product side of the chemical equation for this reaction
    std::string productString() const;

    //! The chemical equation for this reaction
    std::string equation() const;

    //! Set the reactants and products based on the reaction equation. If a Kinetics
    //! object is provided, it is used to check that all reactants and products exist.
    virtual void setEquation(const std::string& equation, const Kinetics* kin=0);

    //! The type of reaction
    virtual std::string type() const;

    //! Calculate the units of the rate constant. These are determined by the units
    //! of the standard concentration of the reactant species' phases and the phase
    //! where the reaction occurs. Sets the value of #rate_units.
    UnitStack calculateRateCoeffUnits(const Kinetics& kin);

    //! Calculate the units of the rate constant.
    //! @deprecated  To be removed after Cantera 3.0. Replaceable by
    //!              calculateRateCoeffUnits.
    UnitStack calculateRateCoeffUnits3(const Kinetics& kin) {
        warn_deprecated("Reaction::calculateRateCoeffUnits3",
            "Deprecated in Cantera 3.0 and to be removed thereafter; replaceable "
            "by calculateRateCoeffUnits.");
        return calculateRateCoeffUnits(kin);
    }

    //! Ensure that the rate constant and other parameters for this reaction are valid.
    //! @since  New in Cantera 3.0.
    virtual void check();

    //! Ensure that the rate constant and other parameters for this reaction are valid.
    //! @deprecated  To be removed after Cantera 3.0. Replaceable by check.
    virtual void validate() {
        warn_deprecated("Reaction::validate",
            "Deprecated in Cantera 3.0 and to be removed thereafter; replaceable "
            "by Reaction::check.");
        check();
    }

    //! Perform validation checks that need access to a complete Kinetics objects, for
    // example to retrieve information about reactant / product species.
    virtual void validate(Kinetics& kin) {
        if (m_rate) {
            m_rate->validate(equation(), kin);
        }
    }

    //! Return the parameters such that an identical Reaction could be reconstructed
    //! using the newReaction() function. Behavior specific to derived classes is
    //! handled by the getParameters() method.
    //! @param withInput  If true, include additional input data fields associated
    //!   with the object, such as user-defined fields from a YAML input file, as
    //!   contained in the #input attribute.
    AnyMap parameters(bool withInput=true) const;

    //! Set up reaction based on AnyMap *node*
    virtual void setParameters(const AnyMap& node, const Kinetics& kin);

    //! Get validity flag of reaction
    bool valid() const {
        return m_valid;
    }

    //! Set validity flag of reaction
    void setValid(bool valid) {
        m_valid = valid;
    }

    //! Check that the specified reaction is balanced (same number of atoms for
    //! each element in the reactants and products). Raises an exception if the
    //! reaction is not balanced. Used by checkSpecies.
    //! @param kin  Kinetics object
    void checkBalance(const Kinetics& kin) const;

    //! Verify that all species involved in the reaction are defined in the Kinetics
    //! object. The function returns true if all species are found, and raises an
    //! exception unless the kinetics object is configured to skip undeclared species,
    //! in which case false is returned.
    //! @param kin  Kinetics object
    bool checkSpecies(const Kinetics& kin) const;

    //! Check whether reaction uses electrochemistry
    //! @param kin  Kinetics object
    bool usesElectrochemistry(const Kinetics& kin) const;

    //! Reactant species and stoichiometric coefficients
    Composition reactants;

    //! Product species and stoichiometric coefficients
    Composition products;

    //! Forward reaction order with respect to specific species. By default,
    //! mass-action kinetics is assumed, with the reaction order equal to each
    //! reactant's stoichiometric coefficient.
    Composition orders;

    //! An identification string for the reaction, used in some filtering
    //! operations
    std::string id;

    //! True if the current reaction is reversible. False otherwise
    bool reversible = true;

    //! True if the current reaction is marked as duplicate
    bool duplicate = false;

    //! True if reaction orders can be specified for non-reactant species.
    //! Default is `false`.
    bool allow_nonreactant_orders = false;

    //! True if negative reaction orders are allowed. Default is `false`.
    bool allow_negative_orders = false;

    //! Input data used for specific models
    AnyMap input;

    //! The units of the rate constant. These are determined by the units of the
    //! standard concentration of the reactant species' phases of the phase
    //! where the reaction occurs.
    Units rate_units{0.0};

    //! Get reaction rate pointer
    shared_ptr<ReactionRate> rate() {
        return m_rate;
    }

    //! Set reaction rate pointer
    void setRate(shared_ptr<ReactionRate> rate);

    //! Get pointer to third-body
    shared_ptr<ThirdBody> thirdBody() {
        return m_third_body;
    }

protected:
    //! Store the parameters of a Reaction needed to reconstruct an identical
    //! object using the newReaction(AnyMap&, Kinetics&) function. Does not
    //! include user-defined fields available in the #input map.
    virtual void getParameters(AnyMap& reactionNode) const;

    //! Flag indicating whether reaction is set up correctly
    bool m_valid = true;

    //! Helper function returning vector of undeclared third body species
    //! and a boolean expression indicating whether the third body is specified.
    //! @note The function is used by the checkSpecies method and only needed as long as
    //! there is no unified approach to handle third body collision partners.
    //! @param kin  Kinetics object
    virtual std::pair<std::vector<std::string>, bool>
        undeclaredThirdBodies(const Kinetics& kin) const;

    //! Reaction rate used by generic reactions
    shared_ptr<ReactionRate> m_rate;

    //! Relative efficiencies of third-body species in enhancing the reaction rate
    //! (if applicable)
    shared_ptr<ThirdBody> m_third_body;
};


//! A class for managing third-body efficiencies, including default values
class ThirdBody
{
public:
    explicit ThirdBody(double default_efficiency=1.0);

    ThirdBody(const AnyMap& node);

    //! Set third-body efficiencies from AnyMap *node*
    void setEfficiencies(const AnyMap& node);

    //! Get the third-body efficiency for species *k*
    double efficiency(const std::string& k) const;

    //! Name or placeholder of third body collider, for example `+ M`
    std::string collider() const;

    //! Map of species to third body efficiency
    Composition efficiencies;

    //! The default third body efficiency for species not listed in #efficiencies.
    double default_efficiency;

    //! Input explicitly specifies collision partner
    bool specified_collision_partner = false;

    //! Third body is used by law of mass action
    //! (`true` for three-body reactions, `false` for falloff reactions)
    bool mass_action = true;
};


//! A reaction with a non-reacting third body "M" that acts to add or remove
//! energy from the reacting species
class ThreeBodyReaction : public Reaction
{
public:
    ThreeBodyReaction();
    ThreeBodyReaction(const Composition& reactants, const Composition& products,
                      const ArrheniusRate& rate, const ThirdBody& tbody);

    ThreeBodyReaction(const AnyMap& node, const Kinetics& kin);

    virtual std::string type() const {
        return "three-body";
    }

    virtual void setEquation(const std::string& equation, const Kinetics* kin=0);
    bool detectEfficiencies();
    virtual void setParameters(const AnyMap& node, const Kinetics& kin);
    virtual void getParameters(AnyMap& reactionNode) const;
};


//! A falloff reaction that is first-order in [M] at low pressure, like a third-body
//! reaction, but zeroth-order in [M] as pressure increases.
//! In addition, the class supports chemically-activated reactions where the rate
//! decreases as pressure increases due to collisional stabilization of a reaction
//! intermediate; in this case, the forward rate constant is written as being
//! proportional to the low-pressure rate constant.
class FalloffReaction : public Reaction
{
public:
    FalloffReaction();
    FalloffReaction(const Composition& reactants, const Composition& products,
                    const ReactionRate& rate, const ThirdBody& tbody);

    FalloffReaction(const AnyMap& node, const Kinetics& kin);

    virtual std::string type() const;

    virtual void setEquation(const std::string& equation, const Kinetics* kin);
    virtual void setParameters(const AnyMap& node, const Kinetics& kin);
    virtual void getParameters(AnyMap& reactionNode) const;
};


//! Create a new empty Reaction object
/*!
 * @param type string identifying type of reaction.
 */
unique_ptr<Reaction> newReaction(const std::string& type);

//! Create a new Reaction object using the specified parameters
/*!
 * @param rxn_node AnyMap node describing reaction.
 * @param kin kinetics manager
 */
unique_ptr<Reaction> newReaction(const AnyMap& rxn_node,
                                 const Kinetics& kin);

//! Create Reaction objects for each item (an AnyMap) in `items`. The species
//! involved in these reactions must exist in the phases associated with the
//! Kinetics object `kinetics`.
std::vector<shared_ptr<Reaction>> getReactions(const AnyValue& items,
                                               Kinetics& kinetics);

//! Parse reaction equation
void parseReactionEquation(Reaction& R, const std::string& equation,
                           const AnyBase& reactionNode, const Kinetics* kin);

using ThreeBodyReaction3 = ThreeBodyReaction; // @todo: remove after Cantera 3.0
using FalloffReaction3 = FalloffReaction; // @todo: remove after Cantera 3.0

}
#endif
