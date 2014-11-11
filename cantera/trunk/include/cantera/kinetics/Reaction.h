/**
 *  @file Reaction.h
 */

#ifndef CT_REACTION_H
#define CT_REACTION_H

#include "cantera/base/utilities.h"
#include "cantera/kinetics/RxnRates.h"

namespace Cantera
{

class Kinetics;

//! Intermediate class which stores data about a reaction and its rate
//! parameterization so that it can be added to a Kinetics object.
class Reaction
{
public:
    Reaction(int type, const Composition& reactants,
             const Composition& products);
    virtual ~Reaction() {}

    friend class Kinetics;

    virtual std::string reactantString() { return ""; } //!< @todo: implement
    virtual std::string productString() { return ""; } //!< @todo: implement
    std::string equation() { return ""; } //!< @todo: implement

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

    bool validate; //!< Perform validation of the rate coefficient data

    //! True if the current reaction is marked as duplicate
    bool duplicate;
};


//! A reaction which follows mass-action kinetics with a modified Arrhenius
//! reaction rate.
class ElementaryReaction : public Reaction
{
public:
    ElementaryReaction(const Composition& reactants, const Composition products,
                       const Arrhenius& rate);

    Arrhenius rate;
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
    //! #third_body_efficiencies.
    double default_efficiency;
};


class ThirdBodyReaction : public ElementaryReaction
{
public:
    ThirdBodyReaction(const Composition& reactants, const Composition& products,
                      const Arrhenius& rate, const ThirdBody& tbody);
    virtual std::string reactantString();
    virtual std::string productString();

    ThirdBody third_body;
};


class FalloffReaction : public Reaction
{
public:
    FalloffReaction(const Composition& reactants, const Composition& products,
                    const Arrhenius& low_rate, const Arrhenius& high_rate,
                    const ThirdBody& tbody, int falloff_type,
                    const vector_fp& falloff_params);
    virtual std::string reactantString();
    virtual std::string productString();

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
    ChemicallyActivatedReaction(const Composition& reactants,
        const Composition& products, const Arrhenius& low_rate,
        const Arrhenius& high_rate, const ThirdBody& tbody, int falloff_type,
        const vector_fp& falloff_params);
};


class PlogReaction : public Reaction
{
public:
    PlogReaction(const Composition& reactants, const Composition& products,
                 const Plog& rate);

    Plog rate;
};


class ChebyshevReaction : public Reaction
{
public:
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


class InterfaceReaction : public Reaction
{
public:
    InterfaceReaction(const Composition& reactants, const Composition& products,
                      const Arrhenius& rate);

    //! Adjustments to the Arrhenius rate expression dependent on surface
    //! species coverages. Three coverage parameters (a, E, m) are used for each
    //! species on which the rate depends. See SurfaceArrhenius for details on
    //! the parameterization.
    std::map<std::string, CoverageDependency> coverage_deps;

    //! The rate coefficient, without taking into account the coverage
    //! dependencies.
    Arrhenius rate;
};


class ElectrochemicalReaction : public InterfaceReaction
{
public:
    ElectrochemicalReaction(const Composition& reactants,
                            const Composition& products, const Arrhenius& rate);

    //! Film Resistivity value
    /*!
     *  Only valid for Butler-Volmer formulations. Units are in ohms m2.
     *  Default = 0.0 ohms m2.
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
};

}

#endif
