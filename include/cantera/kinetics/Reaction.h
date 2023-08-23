/**
 *  @file Reaction.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTION_H
#define CT_REACTION_H

#include "cantera/base/AnyMap.h"
#include "cantera/base/Units.h"
#include "ReactionRate.h"

namespace Cantera
{

class Kinetics;
class ThirdBody;
class ArrheniusRate; // @todo  Remove after Cantera 3.0
class FalloffRate; // @todo  Remove after Cantera 3.0

//! Abstract base class which stores data about a reaction and its rate
//! parameterization so that it can be added to a Kinetics object.
//! @ingroup reactionGroup
class Reaction
{
public:
    Reaction() {}
    Reaction(const Composition& reactants, const Composition& products,
             shared_ptr<ReactionRate> rate, shared_ptr<ThirdBody> tbody=nullptr);
    Reaction(const string& equation,
             shared_ptr<ReactionRate> rate, shared_ptr<ThirdBody> tbody=nullptr);

    //! Construct a Reaction and corresponding ReactionRate based on AnyMap (YAML)
    //! input.
    Reaction(const AnyMap& node, const Kinetics& kin);

    virtual ~Reaction() {}

    //! The reactant side of the chemical equation for this reaction
    string reactantString() const;

    //! The product side of the chemical equation for this reaction
    string productString() const;

    //! The chemical equation for this reaction
    string equation() const;

    //! Set the reactants and products based on the reaction equation. If a Kinetics
    //! object is provided, it is used to check that all reactants and products exist.
    void setEquation(const string& equation, const Kinetics* kin=0);

    //! The type of reaction, including reaction rate information
    string type() const;

    //! Calculate the units of the rate constant. These are determined by the units
    //! of the standard concentration of the reactant species' phases and the phase
    //! where the reaction occurs. Sets the value of #rate_units.
    UnitStack calculateRateCoeffUnits(const Kinetics& kin);

    //! Calculate the units of the rate constant.
    //! @deprecated To be removed after %Cantera 3.0. Replaceable by
    //!              calculateRateCoeffUnits.
    UnitStack calculateRateCoeffUnits3(const Kinetics& kin) {
        warn_deprecated("Reaction::calculateRateCoeffUnits3",
            "Deprecated in Cantera 3.0 and to be removed thereafter; replaceable "
            "by calculateRateCoeffUnits.");
        return calculateRateCoeffUnits(kin);
    }

    //! Ensure that the rate constant and other parameters for this reaction are valid.
    //! @since New in %Cantera 3.0.
    void check();

    //! Ensure that the rate constant and other parameters for this reaction are valid.
    //! @deprecated To be removed after %Cantera 3.0. Replaceable by check.
    void validate() {
        warn_deprecated("Reaction::validate",
            "Deprecated in Cantera 3.0 and to be removed thereafter; replaceable "
            "by Reaction::check.");
        check();
    }

    //! Perform validation checks that need access to a complete Kinetics objects, for
    // example to retrieve information about reactant / product species.
    void validate(Kinetics& kin) {
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
    void setParameters(const AnyMap& node, const Kinetics& kin);

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
    string id;

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

    //! Get pointer to third-body handler
    shared_ptr<ThirdBody> thirdBody() {
        return m_third_body;
    }

    //! Check whether reaction involves third body collider
    //! @since New in %Cantera 3.0.
    bool usesThirdBody() const {
        return bool(m_third_body);
    }

protected:
    //! Store the parameters of a Reaction needed to reconstruct an identical
    //! object using the newReaction(AnyMap&, Kinetics&) function. Does not
    //! include user-defined fields available in the #input map.
    void getParameters(AnyMap& reactionNode) const;

    //! Flag indicating whether reaction is set up correctly
    bool m_valid = true;

    //! Flag indicating that serialization uses explicit type
    bool m_explicit_type = false;

    //! Flag indicating that object was instantiated from reactant/product compositions
    bool m_from_composition = false;

    //! Reaction rate used by generic reactions
    shared_ptr<ReactionRate> m_rate;

    //! Relative efficiencies of third-body species in enhancing the reaction rate
    //! (if applicable)
    shared_ptr<ThirdBody> m_third_body;
};


//! A class for managing third-body efficiencies, including default values
//! @ingroup reactionGroup
class ThirdBody
{
public:
    explicit ThirdBody() {};
    ThirdBody(const string& third_body);
    ThirdBody(const AnyMap& node);

    //! @deprecated To be removed after %Cantera 3.0; instantiate using string instead
    ThirdBody(double default_efficiency);

    //! Name of the third body collider
    //! @since New in %Cantera 3.0
    string name() const {
        return m_name;
    }

    //! Set name of the third body collider
    //! @since New in %Cantera 3.0
    void setName(const string& third_body);

    //! Set third-body efficiencies from AnyMap *node*
    //! @deprecated To be removed after %Cantera 3.0; renamed to setParameters
    void setEfficiencies(const AnyMap& node);

    //! Set third-body efficiencies from AnyMap *node*
    //! @since New in %Cantera 3.0
    void setParameters(const AnyMap& node);

    //! Get third-body efficiencies from AnyMap *node*
    //! @param node  AnyMap receiving serialized parameters
    //! @since New in %Cantera 3.0
    void getParameters(AnyMap& node) const;

    //! Get the third-body efficiency for species *k*
    double efficiency(const string& k) const;

    //! Suffix representing the third body collider in reaction equation, for example
    //! `+ M` or `(+M)`
    //! @since New in %Cantera 3.0
    string collider() const;

    //! Verify that all species involved in collision efficiencies are defined in the
    //! Kinetics object. The function returns true if all species are found, and raises
    //! an exception unless the Kinetics object is configured to skip undeclared
    //! species, in which case false is returned.
    //! @param rxn  Reaction object
    //! @param kin  Kinetics object
    //! @since New in %Cantera 3.0
    bool checkSpecies(const Reaction& rxn, const Kinetics& kin) const;

    //! Map of species to third body efficiency
    Composition efficiencies;

    //! The default third body efficiency for species not listed in #efficiencies.
    double default_efficiency = 1.;

    //! Third body is used by law of mass action
    //! (`true` for three-body reactions, `false` for falloff reactions)
    bool mass_action = true;

    //! Flag indicating whether third body requires explicit serialization
    bool explicit_3rd = false;

protected:
    //! Name of the third body collider
    string m_name = "M";
};


//! A reaction with a non-reacting third body "M" that acts to add or remove
//! energy from the reacting species
//! @deprecated To be removed after %Cantera 3.0. Merged with Reaction
class ThreeBodyReaction : public Reaction
{
public:
    ThreeBodyReaction();
    ThreeBodyReaction(const Composition& reactants, const Composition& products,
                      const ArrheniusRate& rate, const ThirdBody& tbody);

    ThreeBodyReaction(const AnyMap& node, const Kinetics& kin);
};


//! A falloff reaction that is first-order in [M] at low pressure, like a third-body
//! reaction, but zeroth-order in [M] as pressure increases.
//! In addition, the class supports chemically-activated reactions where the rate
//! decreases as pressure increases due to collisional stabilization of a reaction
//! intermediate; in this case, the forward rate constant is written as being
//! proportional to the low-pressure rate constant.
//! @deprecated To be removed after %Cantera 3.0. Merged with Reaction
class FalloffReaction : public Reaction
{
public:
    FalloffReaction();
    FalloffReaction(const Composition& reactants, const Composition& products,
                    const FalloffRate& rate, const ThirdBody& tbody);

    FalloffReaction(const AnyMap& node, const Kinetics& kin);
};


//! Create a new empty Reaction object
/*!
 * @param type string identifying type of reaction.
 */
unique_ptr<Reaction> newReaction(const string& type);

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
vector<shared_ptr<Reaction>> getReactions(const AnyValue& items, Kinetics& kinetics);

//! Parse reaction equation
void parseReactionEquation(Reaction& R, const string& equation,
                           const AnyBase& reactionNode, const Kinetics* kin);

using ThreeBodyReaction3 = ThreeBodyReaction; // @todo: remove after Cantera 3.0
using FalloffReaction3 = FalloffReaction; // @todo: remove after Cantera 3.0

}
#endif
