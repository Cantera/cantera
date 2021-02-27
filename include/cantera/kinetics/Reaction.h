/**
 *  @file Reaction.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTION_H
#define CT_REACTION_H

#include "cantera/base/AnyMap.h"
#include "cantera/kinetics/RxnRates.h"

namespace Cantera
{

class Kinetics;
class Falloff;
class XML_Node;

//! Intermediate class which stores data about a reaction and its rate
//! parameterization so that it can be added to a Kinetics object.
class Reaction
{
public:
    explicit Reaction();
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
    virtual std::string type() const {
        return "<reaction>";
    }

    //! Ensure that the rate constant and other parameters for this reaction are
    //! valid.
    virtual void validate();

    //! Get validity flag of reaction
    bool valid() const {
        return m_valid;
    }

    //! Set validity flag of reaction
    void set_valid(bool valid) {
        m_valid = valid;
    }

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

protected:
    //! Flag indicating whether reaction is set up correctly
    bool m_valid;
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

    virtual std::string type() const {
        return "elementary";
    }

    Arrhenius rate;
    bool allow_negative_pre_exponential_factor;
};

//! A class for managing third-body efficiencies, including default values
class ThirdBody
{
public:
    explicit ThirdBody(double default_efficiency=1.0);

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

    //! Relative efficiencies of third-body species in enhancing the reaction
    //! rate.
    ThirdBody third_body;
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

    virtual std::string type() const {
        return "Chebyshev";
    }

    ChebyshevRate rate;
};

//! A reaction which follows mass-action kinetics with a custom reaction rate
//! defined in Python.
class CustomPyReaction : public Reaction
{
public:
    CustomPyReaction();
    CustomPyReaction(const Composition& reactants, const Composition& products,
                     const CustomPyRate& rate);

    virtual std::string type() const {
        return "custom-Python";
    }

    CustomPyRate rate;
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

    //! Forward value of the apparent Electrochemical transfer coefficient
    doublereal beta;

    bool exchange_current_density_formulation;
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
void setupElementaryReaction(ElementaryReaction&, const XML_Node&);
void setupElementaryReaction(ElementaryReaction&, const AnyMap&,
                             const Kinetics&);

void setupThreeBodyReaction(ThreeBodyReaction&, const XML_Node&);
void setupThreeBodyReaction(ThreeBodyReaction&, const AnyMap&,
                            const Kinetics&);

void setupFalloffReaction(FalloffReaction&, const XML_Node&);
void setupFalloffReaction(FalloffReaction&, const AnyMap&,
                          const Kinetics&);

void setupChemicallyActivatedReaction(ChemicallyActivatedReaction&,
                                      const XML_Node&);

void setupPlogReaction(PlogReaction&, const XML_Node&);
void setupPlogReaction(PlogReaction&, const AnyMap&, const Kinetics&);

void setupChebyshevReaction(ChebyshevReaction&, const XML_Node&);
void setupChebyshevReaction(ChebyshevReaction&, const AnyMap&,
                            const Kinetics&);

void setupCustomPyReaction(CustomPyReaction&, const AnyMap&,
                           const Kinetics&);

void setupInterfaceReaction(InterfaceReaction&, const XML_Node&);
void setupInterfaceReaction(InterfaceReaction&, const AnyMap&,
                            const Kinetics&);

void setupElectrochemicalReaction(ElectrochemicalReaction&,
                                  const XML_Node&);
void setupElectrochemicalReaction(ElectrochemicalReaction&,
                                  const AnyMap&, const Kinetics&);

}

#endif
