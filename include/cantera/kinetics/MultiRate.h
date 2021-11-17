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
    CT_DEFINE_HAS_MEMBER(has_update, updateFromStruct)
    CT_DEFINE_HAS_MEMBER(has_ddT, ddTFromStruct)

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
            kf[rxn.first] = rxn.second.evalFromStruct(m_shared);
        }
    }

    virtual void update(double T) override {
        // update common data once for each reaction type
        m_shared.update(T);
        _updateRates();
    }

    virtual void update(double T, double P) override {
        // update common data once for each reaction type
        m_shared.update(T, P);
        _updateRates();
    }

    virtual void update(const ThermoPhase& bulk, const Kinetics& kin) override {
        // update common data once for each reaction type
        m_shared.update(bulk, kin);
        _updateRates();
    }

    virtual double evalSingle(ReactionRate& rate) override
    {
        RateType& R = static_cast<RateType&>(rate);
        _updateRate(R);
        return R.evalFromStruct(m_shared);
    }

    virtual double ddTSingle(ReactionRate& rate) override
    {
        RateType& R = static_cast<RateType&>(rate);
        _updateRate(R);
        return _get_ddT(R);
    }

protected:
    template <typename T=RateType, typename std::enable_if<has_update<T>::value, bool>::type = true>
    void _updateRates() {
        // Update reaction-specific data for each reaction. This loop is efficient as
        // all calls are de-virtualized and all rate objects are contiguous in memory
        for (auto& rxn : m_rxn_rates) {
            rxn.second.updateFromStruct(m_shared);
        }
    }

    template <typename T=RateType, typename std::enable_if<!has_update<T>::value, bool>::type = true>
    void _updateRates() {
        // Do nothing if the rate has no update method
    }

    template <typename T=RateType, typename std::enable_if<has_update<T>::value, bool>::type = true>
    void _updateRate(RateType& rate) {
        rate.updateFromStruct(m_shared);
    }

    template <typename T=RateType, typename std::enable_if<!has_update<T>::value, bool>::type = true>
    void _updateRate(RateType& rate) {
        // Do nothing if the rate has no update method
    }

    template <typename T=RateType, typename std::enable_if<has_ddT<T>::value, bool>::type = true>
    double _get_ddT(RateType& rate) {
        return rate.ddTFromStruct(m_shared);
    }

    template <typename T=RateType, typename std::enable_if<!has_ddT<T>::value, bool>::type = true>
    double _get_ddT(RateType& rate) {
        throw NotImplementedError("ReactionRate::ddTFromStruct", "For rate of type {}", rate.type());
    }

    //! Vector of pairs of reaction rates indices and reaction rates
    std::vector<std::pair<size_t, RateType>> m_rxn_rates;
    std::map<size_t, size_t> m_indices; //! Mapping of indices
    DataType m_shared;
};

}

#endif
