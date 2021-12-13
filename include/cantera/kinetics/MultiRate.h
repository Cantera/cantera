/**
 * @file MultiRate.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MULTIRATE_H
#define CT_MULTIRATE_H

#include "ReactionRate.h"
#include "MultiRateBase.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

//! A class template handling all reaction rates specific to `BulkKinetics`.
template <class RateType, class DataType>
class MultiBulkRate final : public MultiRateBase
{
    CT_DEFINE_HAS_MEMBER(has_ddT, ddTScaledFromStruct)

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
        m_shared.invalidateCache();
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
        m_shared.invalidateCache();
        if (m_indices.find(rxn_index) != m_indices.end()) {
            size_t j = m_indices[rxn_index];
            m_rxn_rates.at(j).second = dynamic_cast<RateType&>(rate);
            return true;
        }
        return false;
    }

    virtual void resize(size_t n_species, size_t n_reactions) override {
        m_shared.resize(n_species, n_reactions);
        m_shared.invalidateCache();
    }

    virtual void getRateConstants(double* kf) override {
        for (auto& rxn : m_rxn_rates) {
            kf[rxn.first] = rxn.second.evalFromStruct(m_shared);
        }
    }

    virtual void update(double T) override {
        m_shared.update(T);
    }

    virtual void update(double T, double extra) override {
        m_shared.update(T, extra);
    }

    virtual bool update(const ThermoPhase& bulk, const Kinetics& kin) override {
        return m_shared.update(bulk, kin);
    }

    virtual double evalSingle(ReactionRate& rate) override
    {
        RateType& R = static_cast<RateType&>(rate);
        return R.evalFromStruct(m_shared);
    }

protected:
    //! Vector of pairs of reaction rates indices and reaction rates
    std::vector<std::pair<size_t, RateType>> m_rxn_rates;
    std::map<size_t, size_t> m_indices; //! Mapping of indices
    DataType m_shared;
};

}

#endif
