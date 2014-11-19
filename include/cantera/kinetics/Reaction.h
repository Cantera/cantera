/**
 *  @file Reaction.h
 */

#ifndef CT_REACTION_H
#define CT_REACTION_H

#include "cantera/base/utilities.h"
#include "cantera/base/smart_ptr.h"
#include "cantera/kinetics/RxnRates.h"

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

    virtual std::string reactantString() const;
    virtual std::string productString() const;
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


class ThirdBodyReaction : public ElementaryReaction
{
public:
    ThirdBodyReaction();
    ThirdBodyReaction(const Composition& reactants, const Composition& products,
                      const Arrhenius& rate, const ThirdBody& tbody);
    virtual std::string reactantString() const;
    virtual std::string productString() const;

    ThirdBody third_body;
};


class FalloffReaction : public Reaction
{
public:
    FalloffReaction();
    FalloffReaction(const Composition& reactants, const Composition& products,
                    const Arrhenius& low_rate, const Arrhenius& high_rate,
                    const ThirdBody& tbody, int falloff_type,
                    const vector_fp& falloff_params);
    virtual std::string reactantString() const;
    virtual std::string productString() const;
    virtual void validate();

    Arrhenius low_rate;
    Arrhenius high_rate;
    ThirdBody third_body;

    //! Type of falloff parameterization to use. Values are defined in
    //! reaction_defs.h, with names ending in `FALLOFF`.
    int falloff_type;

    //! Values used in the falloff parameterization. Meaning of each parameter
    //! depends on #falloff_type.
    vector_fp falloff_parameters;
};


class ChemicallyActivatedReaction : public FalloffReaction
{
public:
    ChemicallyActivatedReaction();
    ChemicallyActivatedReaction(const Composition& reactants,
        const Composition& products, const Arrhenius& low_rate,
        const Arrhenius& high_rate, const ThirdBody& tbody, int falloff_type,
        const vector_fp& falloff_params);
};


class PlogReaction : public Reaction
{
public:
    PlogReaction();
    PlogReaction(const Composition& reactants, const Composition& products,
                 const Plog& rate);
    virtual void validate();
    Plog rate;
};


class ChebyshevReaction : public Reaction
{
public:
    ChebyshevReaction();
    ChebyshevReaction(const Composition& reactants, const Composition& products,
                      const ChebyshevRate& rate);

    ChebyshevRate rate;
};


struct CoverageDependency
{
    CoverageDependency(double a_, double E_, double m_) : a(a_), E(E_), m(m_) {}
    CoverageDependency() {}
    double a;
    double E;
    double m;
};


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

    //! Power of the equilibrium constant within the Affinity representation
    /*!
     *  Only valid for Affinity representation. Default = 1.0.
     */
    doublereal equilibrium_constant_power;

    //! Power of the "One minus Affinity" term within the Affinity representation
    /*!
     *   Only valid for Affinity representation. Default = 1.0.
     */
    doublereal affinity_power;

    //! Forward value of the apparent Electrochemical transfer coefficient
    doublereal beta;

    bool exchange_current_density_formulation;
};

//! Create a new Reaction object for the reaction defined in `rxn_node`
shared_ptr<Reaction> newReaction(const XML_Node& rxn_node);
}

#endif
