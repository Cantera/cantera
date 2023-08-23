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

//! A class template handling ReactionRate specializations.
//! @ingroup rateEvaluators
template <class RateType, class DataType>
class MultiRate final : public MultiRateBase
{
    CT_DEFINE_HAS_MEMBER(has_update, updateFromStruct)
    CT_DEFINE_HAS_MEMBER(has_ddT, ddTScaledFromStruct)
    CT_DEFINE_HAS_MEMBER(has_ddP, perturbPressure)
    CT_DEFINE_HAS_MEMBER(has_ddM, perturbThirdBodies)

public:
    string type() override {
        if (!m_rxn_rates.size()) {
            throw CanteraError("MultiRate::type",
                 "Cannot determine type of empty rate handler.");
        }
        return m_rxn_rates.at(0).second.type();
    }

    void add(size_t rxn_index, ReactionRate& rate) override {
        m_indices[rxn_index] = m_rxn_rates.size();
        m_rxn_rates.emplace_back(rxn_index, dynamic_cast<RateType&>(rate));
        m_shared.invalidateCache();
    }

    bool replace(size_t rxn_index, ReactionRate& rate) override {
        if (!m_rxn_rates.size()) {
            throw CanteraError("MultiRate::replace",
                 "Invalid operation: cannot replace rate object "
                 "in empty rate handler.");
        }
        if (rate.type() != type()) {
            throw CanteraError("MultiRate::replace",
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

    void resize(size_t nSpecies, size_t nReactions, size_t nPhases) override {
        m_shared.resize(nSpecies, nReactions, nPhases);
        m_shared.invalidateCache();
    }

    void getRateConstants(double* kf) override {
        for (auto& [iRxn, rate] : m_rxn_rates) {
            kf[iRxn] = rate.evalFromStruct(m_shared);
        }
    }

    void processRateConstants_ddT(double* rop, const double* kf, double deltaT) override
    {
        if constexpr (has_ddT<RateType>::value) {
            for (const auto& [iRxn, rate] : m_rxn_rates) {
                rop[iRxn] *= rate.ddTScaledFromStruct(m_shared);
            }
        } else {
            // perturb conditions
            double dTinv = 1. / (m_shared.temperature * deltaT);
            m_shared.perturbTemperature(deltaT);
            _update();

            // apply numerical derivative
            for (auto& [iRxn, rate] : m_rxn_rates) {
                if (kf[iRxn] != 0.) {
                    double k1 = rate.evalFromStruct(m_shared);
                    rop[iRxn] *= dTinv * (k1 / kf[iRxn] - 1.);
                } // else not needed: derivative is already zero
            }

            // revert changes
            m_shared.restore();
            _update();
        }
    }

    void processRateConstants_ddP(double* rop, const double* kf, double deltaP) override
    {
        if constexpr (has_ddP<DataType>::value) {
            double dPinv = 1. / (m_shared.pressure * deltaP);
            m_shared.perturbPressure(deltaP);
            _update();

            for (auto& [iRxn, rate] : m_rxn_rates) {
                if (kf[iRxn] != 0.) {
                    double k1 = rate.evalFromStruct(m_shared);
                    rop[iRxn] *= dPinv * (k1 / kf[iRxn] - 1.);
                } // else not needed: derivative is already zero
            }

            // revert changes
            m_shared.restore();
            _update();
        } else {
            for (const auto& [iRxn, rate] : m_rxn_rates) {
                rop[iRxn] = 0.;
            }
        }
    }

    void processRateConstants_ddM(double* rop, const double* kf, double deltaM,
                                  bool overwrite=true) override
    {
        if constexpr (has_ddM<DataType>::value) {
            double dMinv = 1. / deltaM;
            m_shared.perturbThirdBodies(deltaM);
            _update();

            for (auto& [iRxn, rate] : m_rxn_rates) {
                if (kf[iRxn] != 0. && m_shared.conc_3b[iRxn] > 0.) {
                    double k1 = rate.evalFromStruct(m_shared);
                    rop[iRxn] *= dMinv * (k1 / kf[iRxn] - 1.);
                    rop[iRxn] /= m_shared.conc_3b[iRxn];
                } else {
                    rop[iRxn] = 0.;
                }
            }

            // revert changes
            m_shared.restore();
            _update();
        } else {
            if (!overwrite) {
                // do not overwrite existing entries
                return;
            }
            for (const auto& [iRxn, rate] : m_rxn_rates) {
                rop[iRxn] = 0.;
            }
        }
    }

    void update(double T) override {
        m_shared.update(T);
        _update();
    }

    void update(double T, double extra) override {
        m_shared.update(T, extra);
        _update();
    }

    void update(double T, const vector<double>& extra) override {
        m_shared.update(T, extra);
        _update();
    }

    bool update(const ThermoPhase& phase, const Kinetics& kin) override {
        bool changed = m_shared.update(phase, kin);
        if (changed) {
            // call helper function only if needed: implementation depends on whether
            // ReactionRate::updateFromStruct is defined
            _update();
        }
        return changed;
    }

    double evalSingle(ReactionRate& rate) override {
        RateType& R = static_cast<RateType&>(rate);
        if constexpr (has_update<RateType>::value) {
            R.updateFromStruct(m_shared);
        }
        return R.evalFromStruct(m_shared);
    }

    //! Access the underlying shared data object. Used for setting up
    //! ReactionDataDelegator instances.
    DataType& sharedData() {
        return m_shared;
    }

protected:
    //! Helper function to process updates
    void _update() {
        if constexpr (has_update<RateType>::value) {
            for (auto& [i, rxn] : m_rxn_rates) {
                rxn.updateFromStruct(m_shared);
            }
        }
    }

    //! Vector of pairs of reaction rates indices and reaction rates
    vector<pair<size_t, RateType>> m_rxn_rates;
    map<size_t, size_t> m_indices; //! Mapping of indices
    DataType m_shared;
};

}

#endif
