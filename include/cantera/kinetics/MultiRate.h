/**
 * @file MultiRate.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MULTIRATE_H
#define CT_MULTIRATE_H

#include "ReactionRate.h"
#include "MultiRateBase.h"

namespace Cantera
{

//! A class template handling all reaction rates specific to `BulkKinetics`.
template <class RateType, class DataType>
class MultiBulkRate final : public MultiRateBase
{
public:
    virtual std::string type() override {
        if (!m_rxn_rates.size()) {
            throw CanteraError("MultiBulkRate::type",
                 "Cannot determine type of empty rate handler.");
        }
        return m_rxn_rates.at(0).second.type();
    }

    virtual void add(const size_t rxn_index, ReactionRate& rate) override {
        m_indices[rxn_index] = m_rxn_rates.size();
        m_rxn_rates.emplace_back(rxn_index, dynamic_cast<RateType&>(rate));
    }

    virtual bool replace(const size_t rxn_index, ReactionRate& rate) override {
        if (!m_rxn_rates.size()) {
            throw CanteraError("MultiBulkRate::replace",
                 "Invalid operation: cannot replace rate object "
                 "in empty rate handler.");
        }
        if (rate.type() != type()) {
            throw CanteraError("MultiBulkRate::replace",
                 "Invalid operation: cannot replace rate object of type '{}' "
                 "with a new rate of type '{}'.", type(), rate.type());
        }
        if (m_indices.find(rxn_index) != m_indices.end()) {
            size_t j = m_indices[rxn_index];
            m_rxn_rates.at(j).second = dynamic_cast<RateType&>(rate);
            return true;
        }
        return false;
    }

    virtual void resize(size_t n_species, size_t n_reactions) override {
        m_shared.resize(n_species, n_reactions);
    }

    virtual void getRateConstants(double* kf) override {
        for (auto& rxn : m_rxn_rates) {
            kf[rxn.first] = rxn.second.eval(m_shared);
        }
    }

    virtual void update(double T) override {
        // update common data once for each reaction type
        m_shared.update(T);
        if (!m_rxn_rates.empty() && m_rxn_rates[0].second.usesUpdate()) {
            // update reaction-specific data for each reaction. This loop
            // is efficient as all function calls are de-virtualized, and
            // all of the rate objects are contiguous in memory
            for (auto& rxn : m_rxn_rates) {
                rxn.second.update(m_shared);
            }
        }
    }

    virtual void update(double T, double P) override {
        // update common data once for each reaction type
        m_shared.update(T, P);
        if (!m_rxn_rates.empty() && m_rxn_rates[0].second.usesUpdate()) {
            // update reaction-specific data for each reaction. This loop
            // is efficient as all function calls are de-virtualized, and
            // all of the rate objects are contiguous in memory
            for (auto& rxn : m_rxn_rates) {
                rxn.second.update(m_shared);
            }
        }
    }

    virtual void update(const ThermoPhase& bulk, const Kinetics& kin) override {
        // update common data once for each reaction type
        m_shared.update(bulk, kin);
        if (!m_rxn_rates.empty() && m_rxn_rates[0].second.usesUpdate()) {
            // update reaction-specific data for each reaction. This loop
            // is efficient as all function calls are de-virtualized, and
            // all of the rate objects are contiguous in memory
            for (auto& rxn : m_rxn_rates) {
                rxn.second.update(m_shared);
            }
        }
    }

    virtual double evalSingle(ReactionRate& rate) override
    {
        RateType& R = static_cast<RateType&>(rate);
        if (R.usesUpdate()) {
            R.update(m_shared);
        }
        return R.eval(m_shared);
    }

    virtual double ddTSingle(ReactionRate& rate) override
    {
        RateType& R = static_cast<RateType&>(rate);
        if (R.usesUpdate()) {
            R.update(m_shared);
        }
        return R.ddT(m_shared);
    }

protected:
    //! Vector of pairs of reaction rates indices and reaction rates
    std::vector<std::pair<size_t, RateType>> m_rxn_rates;
    std::map<size_t, size_t> m_indices; //! Mapping of indices
    DataType m_shared;
};

}

#endif
