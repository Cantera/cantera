/**
 * @file MultiRate.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MULTIRATE_H
#define CT_MULTIRATE_H

#include "cantera/kinetics/ReactionRate.h"

namespace Cantera
{

//! An abstract base class for evaluating all reactions of a particular type.
/**
 * Because this class has no template parameters, the `Kinetics` object
 * can store all of these rate coefficient evaluators as a
 * `vector<shared_ptr<MultiRateBase>>`.
 *
 * @todo At the moment, implemented methods are specific to `BulkKinetics`,
 *     which can be updated using information of a single `ThermoPhase`.
 *     `InterfaceKinetics` will require access to an entire `Kinetics` object
 *     or the underlying `vector<ThermoPhase*>` vector (e.g. `m_thermo`).
 */
class MultiRateBase
{
public:
    virtual ~MultiRateBase() {}

    //! Add reaction rate object to the evaluator
    //! @param rxn_index  index of reaction
    //! @param rate  reaction rate object
    virtual void add(const size_t rxn_index,
                     std::shared_ptr<ReactionRateBase> rate) = 0;

    //! Replace reaction rate object handled by the evaluator
    //! @param rxn_index  index of reaction
    //! @param rate  reaction rate object
    virtual bool replace(const size_t rxn_index,
                         std::shared_ptr<ReactionRateBase> rate) = 0;

    //! Access reaction rate object handled by evaluator
    //! @param rxn_index  index of reaction
    virtual ReactionRateBase& rate(size_t rxn_index) = 0;

    //! Update number of species
    //! @param n_species  number of species
    virtual void resizeSpecies(size_t n_species) = 0;

    //! Evaluate all rate constants handled by the evaluator
    //! @param bulk  object representing bulk phase
    //! @param kf  array of rate constants
    virtual void getRateConstants(const ThermoPhase& bulk,
                                  double* kf, double* concm) const = 0;

    //! Update data common to reaction rates of a specific type
    //! @param bulk  object representing bulk phase
    //! @param concm  effective third-body concentrations
    //! @TODO  enable more generic handling of non-trivial concentration dependencies
    virtual void update(const ThermoPhase& bulk, double* concm) = 0;
};


//! A class template handling all reaction rates specific to `BulkKinetics`.
template <class RateType, class DataType>
class MultiBulkRates final : public MultiRateBase
{
public:
    virtual ~MultiBulkRates() {
        for (size_t i = 0; i < m_bases.size(); i++) {
            m_bases[i]->releaseEvaluator();
        }
    }

    virtual void add(const size_t rxn_index,
                     std::shared_ptr<ReactionRateBase> rate) override
    {
        if (typeid(*rate) != typeid(RateType)) {
            throw CanteraError("MultiBulkRate::add",
                 "Wrong type: cannot add rate object of type '{}'",
                 rate->type());
        } else if (rate->linked()) {
            // ensure there are no dangling objects (unlikely)
            throw CanteraError("MultiBulkRate::add",
                 "Reaction rate is already linked to a reaction rate evaluator;\n"
                 "object needs to release link before it can be added again.");
        }
        size_t j = m_rxn_rates.size();
        m_indices[rxn_index] = j;
        m_rxn_rates.emplace_back(rxn_index, dynamic_cast<RateType&>(*rate));
        // keep link to original object
        m_bases.push_back(rate);
    }

    virtual bool replace(const size_t rxn_index,
                         std::shared_ptr<ReactionRateBase> rate) override
    {
        if (!m_rxn_rates.size()) {
            throw CanteraError("MultiBulkRate::replace",
                 "Invalid operation: cannot replace rate object "
                 "in empty rate handler.");
        } else if (typeid(*rate) != typeid(RateType)) {
            throw CanteraError("MultiBulkRate::replace",
                 "Invalid operation: cannot replace rate object of type '{}' "
                 "with a new rate of type '{}'.",
                 m_rxn_rates.at(0).second.type(), rate->type());
        } else if (rate->linked()) {
            // ensure there are no dangling objects (unlikely)
            throw CanteraError("MultiBulkRate::replace",
                 "Reaction rate is already linked to a reaction rate "
                 "evaluator. Re-linked reaction rates are not allowed");
        }
        if (m_indices.find(rxn_index) != m_indices.end()) {
            size_t j = m_indices[rxn_index];
            m_rxn_rates.at(j).second = dynamic_cast<RateType&>(*rate);
            // release evaluator from previously used rate object and update
            m_bases[j]->releaseEvaluator();
            m_bases[j] = rate;
            return true;
        }
        return false;
    }

    virtual ReactionRateBase& rate(size_t rxn_index) override
    {
        size_t j = m_indices[rxn_index];
        return dynamic_cast<ReactionRateBase&>(m_rxn_rates.at(j).second);
    }

    virtual void resizeSpecies(size_t n_species)
    {
        m_shared.resizeSpecies(n_species);
    }

    virtual void getRateConstants(const ThermoPhase& bulk,
                                  double* kf, double* concm) const override
    {
        for (const auto& rxn : m_rxn_rates) {
            kf[rxn.first] = rxn.second.eval(m_shared, concm[rxn.first]);
        }
    }

    virtual void update(const ThermoPhase& bulk, double* concm) override
    {
        // update common data once for each reaction type
        m_shared.update(bulk);
        if (RateType::usesUpdate()) {
            // update reaction-specific data for each reaction. This loop
            // is efficient as all function calls are de-virtualized, and
            // all of the rate objects are contiguous in memory
            for (auto& rxn : m_rxn_rates) {
                rxn.second.update(m_shared, concm[rxn.first]);
            }
        }
    }

protected:
    //! Pointers to reaction rate objects managed by Reaction object
    std::vector<shared_ptr<ReactionRateBase>> m_bases;

    //! Vector of pairs of reaction rates indices and reaction rates
    std::vector<std::pair<size_t, RateType>> m_rxn_rates;
    std::map<size_t, size_t> m_indices; //! Mapping of indices
    DataType m_shared;
};

}

#endif
