/**
 *  @file Reaction.h
 */

#ifndef CT_REACTION_H
#define CT_REACTION_H

#include "cantera/base/utilities.h"
#include "cantera/base/smart_ptr.h"
#include "cantera/kinetics/RxnRates.h"
#include "cantera/kinetics/Falloff.h"

namespace Cantera
{

class Kinetics;

//! Intermediate class which stores data about a reaction and its rate
//! parameterization so that it can be added to a Kinetics object.
class Reaction
{
public:
    explicit Reaction(int type);
    Reaction(int type, const Composition& reactants,
             const Composition& products);
    virtual ~Reaction() {}

    //! The reactant side of the chemical equation for this reaction
    virtual std::string reactantString() const;

    //! The product side of the chemical equation for this reaction
    virtual std::string productString() const;

    //! The chemical equation for this reaction
    std::string equation() const;

    //! Ensure that the rate constant and other parameters for this reaction are
    //! valid.
    virtual void validate();

    //! Type of the reaction. The valid types are listed in the file,
    //! reaction_defs.h, with constants ending in `RXN`.
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

    Arrhenius rate;
    bool allow_negative_pre_exponential_factor;
};

//! A class for managing third-body efficiencies, including default values
class ThirdBody
{
public:
    explicit ThirdBody(double default_efficiency=1.0);

    //! Get the third-body efficiency for species *k*
    double efficiency(const std::string& k) const {
        return getValue(efficiencies, k, default_efficiency);
    }

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
};

//! A pressure-dependent reaction parameterized by logarithmically interpolating
//! between Arrhenius rate expressions at various pressures.
class PlogReaction : public Reaction
{
public:
    PlogReaction();
    PlogReaction(const Composition& reactants, const Composition& products,
                 const Plog& rate);
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

    ChebyshevRate rate;
};

//! Modifications to an InterfaceReaction rate based on a surface species
//! coverage.
struct CoverageDependency
{
    //! Constructor
    //! @param a_  modification to the pre-exponential factor [m, kmol, s units]
    //! @param E_  modification to the activation energy [K]
    //! @param m_  modification to the temperature exponent
    CoverageDependency(double a_, double E_, double m_) : a(a_), E(E_), m(m_) {}
    CoverageDependency() {}
    double a; //!< modification to the pre-exponential factor [m, kmol, s units]
    double E; //!< modification to the activation energy [K]
    double m; //!< modification to the temperature exponent
};

//! A reaction occurring on an interface (i.e. a SurfPhase or an EdgePhase)
class InterfaceReaction : public ElementaryReaction
{
public:
    InterfaceReaction();
    InterfaceReaction(const Composition& reactants, const Composition& products,
                      const Arrhenius& rate, bool isStick=false);

    //! Adjustments to the Arrhenius rate expression dependent on surface
    //! species coverages. Three coverage parameters (a, E, m) are used for each
    //! species on which the rate depends. See SurfaceArrhenius for details on
    //! the parameterization.
    std::map<std::string, CoverageDependency> coverage_deps;

    //! Set to true if `rate` is a parameterization of the sticking coefficient
    //! rather than the forward rate constant
    bool is_sticking_coefficient;

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

    //! Film Resistivity value
    /*!
     *  For Butler Volmer reactions, a common addition to the formulation is to
     *  add an electrical resistance to the formulation. The resistance modifies
     *  the electrical current flow in both directions. Only valid for Butler-
     *  Volmer formulations. Units are in ohms m2. Default = 0.0 ohms m2.
     */
    doublereal film_resistivity;

    //! Forward value of the apparent Electrochemical transfer coefficient
    doublereal beta;

    bool exchange_current_density_formulation;
};

//! Create a new Reaction object for the reaction defined in `rxn_node`
shared_ptr<Reaction> newReaction(const XML_Node& rxn_node);

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
std::vector<shared_ptr<Reaction> > getReactions(const XML_Node& node);

}

#endif
