/**
 *  @file Reaction.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTION_H
#define CT_REACTION_H

#include "cantera/base/AnyMap.h"
#include "cantera/kinetics/ReactionRate.h"
#include "cantera/kinetics/RxnRates.h"
#include "cantera/base/Units.h"

namespace Cantera
{

class Kinetics;
class Falloff;
class XML_Node;

//! Abstract base class which stores data about a reaction and its rate
//! parameterization so that it can be added to a Kinetics object.
class Reaction
{
public:
    Reaction();
    Reaction(const Composition& reactants,
             const Composition& products);

    //! @deprecated To be removed after Cantera 2.6.
    explicit Reaction(int type);
    //! @deprecated To be removed after Cantera 2.6.
    Reaction(int type, const Composition& reactants,
             const Composition& products);
    virtual ~Reaction() {}

    //! The reactant side of the chemical equation for this reaction
    virtual std::string reactantString() const;

    //! The product side of the chemical equation for this reaction
    virtual std::string productString() const;

    //! The chemical equation for this reaction
    std::string equation() const;

    //! The type of reaction
    virtual std::string type() const = 0; // pure virtual function

    //! Calculate the units of the rate constant. These are determined by the units
    //! of the standard concentration of the reactant species' phases and the phase
    //! where the reaction occurs. Sets the value of #rate_units.
    virtual void calculateRateCoeffUnits(const Kinetics& kin);

    //! Ensure that the rate constant and other parameters for this reaction are
    //! valid.
    virtual void validate();

    //! Return the parameters such that an identical Reaction could be reconstructed
    //! using the newReaction() function. Behavior specific to derived classes is
    //! handled by the getParameters() method.
    //! @param withInput  If true, include additional input data fields associated
    //!   with the object, such as user-defined fields from a YAML input file, as
    //!   contained in the #input attribute.
    AnyMap parameters(bool withInput=true) const;

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

    //! Type of the reaction. The valid types are listed in the file,
    //! reaction_defs.h, with constants ending in `RXN`.
    /*!
     * @deprecated To be removed in Cantera 2.6.
     *             Superseded by Reaction::type().
     */
    int reaction_type;

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
    bool reversible;

    //! True if the current reaction is marked as duplicate
    bool duplicate;

    //! True if reaction orders can be specified for non-reactant species.
    //! Default is `false`.
    bool allow_nonreactant_orders;

    //! True if negative reaction orders are allowed. Default is `false`.
    bool allow_negative_orders;

    //! Input data used for specific models
    AnyMap input;

    //! The units of the rate constant. These are determined by the units of the
    //! standard concentration of the reactant species' phases of the phase
    //! where the reaction occurs.
    Units rate_units;

protected:
    //! Store the parameters of a Reaction needed to reconstruct an identical
    //! object using the newReaction(AnyMap&, Kinetics&) function. Does not
    //! include user-defined fields available in the #input map.
    virtual void getParameters(AnyMap& reactionNode) const;

    //! Flag indicating whether reaction is set up correctly
    bool m_valid;

    //! @internal  Helper function returning vector of undeclared third body species
    //! and a boolean expression indicating whether the third body is specified.
    //! The function is used by the checkSpecies method and only needed as long as
    //! there is no unified approach to handle third body collision partners.
    //! @param kin  Kinetics object
    virtual std::pair<std::vector<std::string>, bool>
        undeclaredThirdBodies(const Kinetics& kin) const;
};


//! An intermediate class used to avoid naming conflicts of 'rate' member
//! variables and getters (see `ElementaryReaction`, `PlogReaction` and
//! `ChebyshevReaction`).
class Reaction3 : public Reaction
{
public:
    Reaction3() : Reaction() {}
    Reaction3(const Composition& reactants, const Composition& products);

    //! Get reaction rate pointer
    shared_ptr<ReactionRateBase> rate() {
        return m_rate;
    }

    //! Set reaction rate pointer
    void setRate(shared_ptr<ReactionRateBase> rate) {
        m_rate = rate;
    }

    //! Set up reaction based on AnyMap node
    virtual bool setParameters(const AnyMap& node, const Kinetics& kin);

protected:
    //! Reaction rate used by generic reactions
    shared_ptr<ReactionRateBase> m_rate;
};


//! A reaction which follows mass-action kinetics with a modified Arrhenius
//! reaction rate.
class ElementaryReaction : public Reaction
{
public:
    ElementaryReaction();
    ElementaryReaction(const Composition& reactants, const Composition products,
                       const Arrhenius& rate);
    virtual void validate();
    virtual void getParameters(AnyMap& reactionNode) const;

    virtual std::string type() const {
        return "elementary-old";
    }

    Arrhenius rate;
    bool allow_negative_pre_exponential_factor;
};

//! A class for managing third-body efficiencies, including default values
class ThirdBody
{
public:
    explicit ThirdBody(double default_efficiency=1.0);

    ThirdBody(const AnyMap& node);

    //! Set third-body efficiencies from AnyMap *node*
    bool setEfficiencies(const AnyMap& node);

    //! Get the third-body efficiency for species *k*
    double efficiency(const std::string& k) const;

    //! Map of species to third body efficiency
    Composition efficiencies;

    //! The default third body efficiency for species not listed in
    //! #efficiencies.
    double default_efficiency;
};

//! A reaction with a non-reacting third body "M" that acts to add or remove
//! energy from the reacting species
class ThreeBodyReaction : public ElementaryReaction
{
public:
    ThreeBodyReaction();
    ThreeBodyReaction(const Composition& reactants, const Composition& products,
                      const Arrhenius& rate, const ThirdBody& tbody);

    virtual std::string type() const {
        return "three-body";
    }

    virtual std::string reactantString() const;
    virtual std::string productString() const;
    virtual void calculateRateCoeffUnits(const Kinetics& kin);
    virtual void getParameters(AnyMap& reactionNode) const;

    //! Relative efficiencies of third-body species in enhancing the reaction
    //! rate.
    ThirdBody third_body;

    bool specified_collision_partner = false; //!< Input specifies collision partner

protected:
    virtual std::pair<std::vector<std::string>, bool>
        undeclaredThirdBodies(const Kinetics& kin) const;
};

//! A reaction that is first-order in [M] at low pressure, like a third-body
//! reaction, but zeroth-order in [M] as pressure increases.
class FalloffReaction : public Reaction
{
public:
    FalloffReaction();
    FalloffReaction(const Composition& reactants, const Composition& products,
                    const Arrhenius& low_rate, const Arrhenius& high_rate,
                    const ThirdBody& tbody);

    virtual std::string type() const {
        return "falloff";
    }

    virtual std::string reactantString() const;
    virtual std::string productString() const;
    virtual void validate();
    virtual void calculateRateCoeffUnits(const Kinetics& kin);
    virtual void getParameters(AnyMap& reactionNode) const;

    //! The rate constant in the low-pressure limit
    Arrhenius low_rate;

    //! The rate constant in the high-pressure limit
    Arrhenius high_rate;

    //! Relative efficiencies of third-body species in enhancing the reaction rate
    ThirdBody third_body;

    //! Falloff function which determines how low_rate and high_rate are
    //! combined to determine the rate constant for the reaction.
    shared_ptr<Falloff> falloff;

    bool allow_negative_pre_exponential_factor;

    //! The units of the low-pressure rate constant. The units of the
    //! high-pressure rate constant are stored in #rate_units.
    Units low_rate_units;

protected:
    virtual std::pair<std::vector<std::string>, bool> undeclaredThirdBodies(
        const Kinetics& kin) const;
};

//! A reaction where the rate decreases as pressure increases due to collisional
//! stabilization of a reaction intermediate. Like a FalloffReaction, except
//! that the forward rate constant is written as being proportional to the low-
//! pressure rate constant.
class ChemicallyActivatedReaction : public FalloffReaction
{
public:
    ChemicallyActivatedReaction();
    ChemicallyActivatedReaction(const Composition& reactants,
        const Composition& products, const Arrhenius& low_rate,
        const Arrhenius& high_rate, const ThirdBody& tbody);

    virtual std::string type() const {
        return "chemically-activated";
    }

    virtual void calculateRateCoeffUnits(const Kinetics& kin);
    virtual void getParameters(AnyMap& reactionNode) const;
};

//! A pressure-dependent reaction parameterized by logarithmically interpolating
//! between Arrhenius rate expressions at various pressures.
class PlogReaction : public Reaction
{
public:
    PlogReaction();
    PlogReaction(const Composition& reactants, const Composition& products,
                 const Plog& rate);

    virtual std::string type() const {
        return "pressure-dependent-Arrhenius";
    }

    virtual void validate();
    virtual void getParameters(AnyMap& reactionNode) const;

    Plog rate;
};

//! A pressure-dependent reaction parameterized by a bi-variate Chebyshev
//! polynomial in temperature and pressure
class ChebyshevReaction : public Reaction
{
public:
    ChebyshevReaction();
    ChebyshevReaction(const Composition& reactants, const Composition& products,
                      const ChebyshevRate& rate);
    virtual void getParameters(AnyMap& reactionNode) const;

    virtual std::string type() const {
        return "Chebyshev";
    }

    ChebyshevRate rate;
};

//! A reaction which follows mass-action kinetics with a modified Arrhenius
//! reaction rate.
class ElementaryReaction3 : public Reaction3
{
public:
    ElementaryReaction3();
    ElementaryReaction3(const Composition& reactants, const Composition& products,
                        const ArrheniusRate& rate);

    ElementaryReaction3(const AnyMap& node, const Kinetics& kin);

    virtual std::string type() const {
        return "elementary";
    }

    virtual bool setParameters(const AnyMap& node, const Kinetics& kin);
    virtual void getParameters(AnyMap& reactionNode) const;

    virtual void validate();

    bool allow_negative_pre_exponential_factor;
};


//! A reaction with a non-reacting third body "M" that acts to add or remove
//! energy from the reacting species
class ThreeBodyReaction3 : public ElementaryReaction3
{
public:
    ThreeBodyReaction3();
    ThreeBodyReaction3(const Composition& reactants, const Composition& products,
                       const ArrheniusRate& rate, const ThirdBody& tbody);

    ThreeBodyReaction3(const AnyMap& node, const Kinetics& kin);

    virtual std::string type() const {
        return "three-body-new";
    }

    virtual void calculateRateCoeffUnits(const Kinetics& kin);
    virtual bool setParameters(const AnyMap& node, const Kinetics& kin);
    virtual void getParameters(AnyMap& reactionNode) const;

    virtual std::string reactantString() const;
    virtual std::string productString() const;

    //! Relative efficiencies of third-body species in enhancing the reaction
    //! rate.
    ThirdBody third_body;

    bool specified_collision_partner = false; //!< Input specifies collision partner
};


//! A reaction which follows mass-action kinetics with a custom reaction rate
//! defined in Python.
class CustomFunc1Reaction : public Reaction3
{
public:
    CustomFunc1Reaction();

    CustomFunc1Reaction(const AnyMap& node, const Kinetics& kin);

    virtual bool setParameters(const AnyMap& node, const Kinetics& kin);

    virtual std::string type() const {
        return "custom-rate-function";
    }
};


//! A reaction which follows mass-action kinetics with a modified Arrhenius
//! reaction rate.
/**
 * Alternative elementary reaction based on ReactionRate.
 *
 * @warning This class is an experimental part of the %Cantera API and
 *    may be changed or removed without notice.
 */
class TestReaction : public Reaction3
{
public:
    TestReaction();

    TestReaction(const AnyMap& node, const Kinetics& kin);

    virtual std::string type() const {
        return "elementary-new";
    }

    virtual bool setParameters(const AnyMap& node, const Kinetics& kin);
};


//! Modifications to an InterfaceReaction rate based on a surface species
//! coverage.
struct CoverageDependency
{
    //! Constructor
    //! @param a_  coefficient for exponential dependence on coverage [dimensionless]
    //! @param E_  modification to the activation energy [K]
    //! @param m_  exponent for power law dependence on coverage [dimensionless]
    CoverageDependency(double a_, double E_, double m_) : a(a_), E(E_), m(m_) {}
    CoverageDependency() {}
    double a; //!< coefficient for exponential dependence on coverage [dimensionless]
    double E; //!< modification to the activation energy [K]
    double m; //!< exponent for power law dependence on coverage [dimensionless]
};

//! A reaction occurring on an interface (i.e. a SurfPhase or an EdgePhase)
class InterfaceReaction : public ElementaryReaction
{
public:
    InterfaceReaction();
    InterfaceReaction(const Composition& reactants, const Composition& products,
                      const Arrhenius& rate, bool isStick=false);
    virtual void calculateRateCoeffUnits(const Kinetics& kin);
    virtual void getParameters(AnyMap& reactionNode) const;

    virtual std::string type() const {
        return "interface";
    }

    //! Adjustments to the Arrhenius rate expression dependent on surface
    //! species coverages. Three coverage parameters (a, E, m) are used for each
    //! species on which the rate depends. See SurfaceArrhenius for details on
    //! the parameterization.
    std::map<std::string, CoverageDependency> coverage_deps;

    //! Set to true if `rate` is a parameterization of the sticking coefficient
    //! rather than the forward rate constant
    bool is_sticking_coefficient;

    //! Set to true if `rate` is a sticking coefficient which should be
    //! translated into a rate coefficient using the correction factor developed
    //! by Motz & Wise for reactions with high (near-unity) sticking
    //! coefficients. Defaults to 'false'.
    bool use_motz_wise_correction;

    //! For reactions with multiple non-surface species, the sticking species
    //! needs to be explicitly identified.
    std::string sticking_species;
};

//! An interface reaction which involves charged species
class ElectrochemicalReaction : public InterfaceReaction
{
public:
    ElectrochemicalReaction();
    ElectrochemicalReaction(const Composition& reactants,
                            const Composition& products, const Arrhenius& rate);
    virtual void getParameters(AnyMap& reactionNode) const;

    //! Forward value of the apparent Electrochemical transfer coefficient
    doublereal beta;

    bool exchange_current_density_formulation;
};

//! A reaction with rate parameters for Blowers-Masel approximation
class BlowersMaselReaction: public Reaction
{
public:
    BlowersMaselReaction();
    BlowersMaselReaction(const Composition& reactants,
                         const Composition& products, const BlowersMasel& rate);
    virtual void getParameters(AnyMap& reactionNode) const;
    virtual void validate();

    virtual std::string type() const {
        return "Blowers-Masel";
    }

    BlowersMasel rate;

    bool allow_negative_pre_exponential_factor;
};

//! A reaction occurring on an interface (i.e. a SurfPhase or an EdgePhase)
//! with the rate calculated with Blowers-Masel approximation.
class BlowersMaselInterfaceReaction : public BlowersMaselReaction
{
public:
    BlowersMaselInterfaceReaction();
    BlowersMaselInterfaceReaction(const Composition& reactants, const Composition& products,
                      const BlowersMasel& rate, bool isStick=false);
    virtual void getParameters(AnyMap& reactionNode) const;
    virtual void calculateRateCoeffUnits(const Kinetics& kin);

    virtual std::string type() const {
        return "surface-Blowers-Masel";
    }
    //! Adjustments to the Arrhenius rate expression dependent on surface
    //! species coverages. Three coverage parameters (a, E, m) are used for each
    //! species on which the rate depends. See SurfaceArrhenius for details on
    //! the parameterization.
    std::map<std::string, CoverageDependency> coverage_deps;

    //! Set to true if `rate` is a parameterization of the sticking coefficient
    //! rather than the forward rate constant
    bool is_sticking_coefficient;

    //! Set to true if `rate` is a sticking coefficient which should be
    //! translated into a rate coefficient using the correction factor developed
    //! by Motz & Wise for reactions with high (near-unity) sticking
    //! coefficients. Defaults to 'false'.
    bool use_motz_wise_correction;

    //! For reactions with multiple non-surface species, the sticking species
    //! needs to be explicitly identified.
    std::string sticking_species;
};

//! Create Reaction objects for all `<reaction>` nodes in an XML document.
//!
//! The `<reaction>` nodes are assumed to be children of the `<reactionData>`
//! node in an XML document with a `<ctml>` root node, as in the case of XML
//! files produced by conversion from CTI files.
//!
//! This function can be used in combination with get_XML_File() and
//! get_XML_from_string() to get Reaction objects from either a file or a
//! string, respectively, where the string or file is formatted as either CTI
//! or XML.
//!
//! If Reaction objects are being created from a CTI definition that does not
//! contain corresponding phase definitions, then one of the following must be
//! true, or the resulting rate constants will be incorrect:
//!
//!   - The rate constants are expressed in (kmol, meter, second) units
//!   - A `units` directive is included **and** all reactions take place in
//!     bulk (e.g. gas) phases
//!
//! @deprecated The XML input format is deprecated and will be removed in
//!     Cantera 3.0.
std::vector<shared_ptr<Reaction> > getReactions(const XML_Node& node);

//! Create Reaction objects for each item (an AnyMap) in `items`. The species
//! involved in these reactions must exist in the phases associated with the
//! Kinetics object `kinetics`.
std::vector<shared_ptr<Reaction>> getReactions(const AnyValue& items,
                                               Kinetics& kinetics);

//! Parse reaction equation
void parseReactionEquation(Reaction& R, const AnyValue& equation,
                           const Kinetics& kin);

// declarations of setup functions
void setupReaction(Reaction& R, const XML_Node& rxn_node);

void setupElementaryReaction(ElementaryReaction&, const XML_Node&);
//! @internal May be changed without notice in future versions
void setupElementaryReaction(ElementaryReaction&, const AnyMap&,
                             const Kinetics&);

void setupThreeBodyReaction(ThreeBodyReaction&, const XML_Node&);
//! @internal May be changed without notice in future versions
void setupThreeBodyReaction(ThreeBodyReaction&, const AnyMap&,
                            const Kinetics&);

void setupFalloffReaction(FalloffReaction&, const XML_Node&);
//! @internal May be changed without notice in future versions
void setupFalloffReaction(FalloffReaction&, const AnyMap&,
                          const Kinetics&);

void setupChemicallyActivatedReaction(ChemicallyActivatedReaction&,
                                      const XML_Node&);

void setupPlogReaction(PlogReaction&, const XML_Node&);
//! @internal May be changed without notice in future versions
void setupPlogReaction(PlogReaction&, const AnyMap&, const Kinetics&);

void setupChebyshevReaction(ChebyshevReaction&, const XML_Node&);
//! @internal May be changed without notice in future versions
void setupChebyshevReaction(ChebyshevReaction&, const AnyMap&,
                            const Kinetics&);

void setupInterfaceReaction(InterfaceReaction&, const XML_Node&);
//! @internal May be changed without notice in future versions
void setupInterfaceReaction(InterfaceReaction&, const AnyMap&,
                            const Kinetics&);

void setupElectrochemicalReaction(ElectrochemicalReaction&,
                                  const XML_Node&);
//! @internal May be changed without notice in future versions
void setupElectrochemicalReaction(ElectrochemicalReaction&,
                                  const AnyMap&, const Kinetics&);

//! @internal May be changed without notice in future versions
void setupBlowersMaselReaction(BlowersMaselReaction&,
                               const AnyMap&, const Kinetics&);
//! @internal May be changed without notice in future versions
void setupBlowersMaselInterfaceReaction(BlowersMaselInterfaceReaction&,
                                        const AnyMap&, const Kinetics&);
}

#endif
