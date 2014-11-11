#include "cantera/kinetics/Reaction.h"

namespace Cantera
{

Reaction::Reaction(int type, const Composition& reactants_,
                   const Composition& products_)
    : reaction_type(type)
    , reactants(reactants_)
    , products(products_)
    , reversible(true)
    , validate(true)
    , duplicate(false)
{
}

ElementaryReaction::ElementaryReaction(const Composition& reactants_,
                                       const Composition products_,
                                       const Arrhenius& rate_)
    : Reaction(ELEMENTARY_RXN, reactants_, products_)
    , rate(rate_)
{
}

ThirdBody::ThirdBody(double default_eff)
    : default_efficiency(default_eff)
{
}

ThirdBodyReaction::ThirdBodyReaction(const Composition& reactants_,
                                     const Composition& products_,
                                     const Arrhenius& rate_,
                                     const ThirdBody& tbody)
    : ElementaryReaction(reactants_, products_, rate_)
    , third_body(tbody)
{
    reaction_type = THREE_BODY_RXN;
}

std::string ThirdBodyReaction::reactantString() {
    return ElementaryReaction::reactantString() + " + M";
}

std::string ThirdBodyReaction::productString() {
    return ElementaryReaction::productString() + " + M";
}

FalloffReaction::FalloffReaction(
        const Composition& reactants_, const Composition& products_,
        const Arrhenius& low_rate_, const Arrhenius& high_rate_,
        const ThirdBody& tbody, int type,
        const vector_fp& params)
    : Reaction(FALLOFF_RXN, reactants_, products_)
    , low_rate(low_rate_)
    , high_rate(high_rate_)
    , third_body(tbody)
    , falloff_type(type)
    , falloff_parameters(params)
{
}

std::string FalloffReaction::reactantString() {
    return Reaction::reactantString() + " (+M)";
}

std::string FalloffReaction::productString() {
    return Reaction::productString() + " (+M)";
}

ChemicallyActivatedReaction::ChemicallyActivatedReaction(
        const Composition& reactants_, const Composition& products_,
        const Arrhenius& low_rate_, const Arrhenius& high_rate_,
        const ThirdBody& tbody, int falloff_type,
        const vector_fp& falloff_params)
    : FalloffReaction(reactants_, products_, low_rate, high_rate, tbody,
                      falloff_type, falloff_params)
{
    reaction_type = CHEMACT_RXN;
}

PlogReaction::PlogReaction(const Composition& reactants_,
                           const Composition& products_, const Plog& rate_)
    : Reaction(PLOG_RXN, reactants_, products_)
    , rate(rate_)
{
}

ChebyshevReaction::ChebyshevReaction(const Composition& reactants_,
                                     const Composition& products_,
                                     const ChebyshevRate& rate_)
    : Reaction(CHEBYSHEV_RXN, reactants_, products_)
    , rate(rate_)
{
}

InterfaceReaction::InterfaceReaction(const Composition& reactants_,
                                     const Composition& products_,
                                     const Arrhenius& rate_)
    : Reaction(INTERFACE_RXN, reactants_, products_)
    , rate(rate_)
{
}

ElectrochemicalReaction::ElectrochemicalReaction(const Composition& reactants_,
                                                 const Composition& products_,
                                                 const Arrhenius& rate_)
    : InterfaceReaction(reactants_, products_, rate_)
    , film_resistivity(0.0)
    , equilibrium_constant_power(1.0)
    , affinity_power(1.0)
    , beta(0.0)
{
}

}
