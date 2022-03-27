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
template <class RateType, class DataType>
class MultiRate final : public MultiRateBase
{
    CT_DEFINE_HAS_MEMBER(has_update, updateFromStruct)
    CT_DEFINE_HAS_MEMBER(has_ddT, ddTScaledFromStruct)
    CT_DEFINE_HAS_MEMBER(has_ddP, perturbPressure)
    CT_DEFINE_HAS_MEMBER(has_ddM, perturbThirdBodies)

public:
    virtual std::string type() override {
        if (!m_rxn_rates.size()) {
            throw CanteraError("MultiRate::type",
                 "Cannot determine type of empty rate handler.");
        }
        return m_rxn_rates.at(0).second.type();
    }

    virtual void add(size_t rxn_index, ReactionRate& rate) override {
        m_indices[rxn_index] = m_rxn_rates.size();
        m_rxn_rates.emplace_back(rxn_index, dynamic_cast<RateType&>(rate));
        m_shared.invalidateCache();
    }

    virtual bool replace(size_t rxn_index, ReactionRate& rate) override {
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

    virtual void resize(size_t nSpecies, size_t nReactions, size_t nPhases) override {
        m_shared.resize(nSpecies, nReactions, nPhases);
        m_shared.invalidateCache();
    }

    virtual void getRateConstants(double* kf) override {
        for (auto& rxn : m_rxn_rates) {
            kf[rxn.first] = rxn.second.evalFromStruct(m_shared);
        }
    }

    virtual void processRateConstants_ddT(double* rop,
                                          const double* kf,
                                          double deltaT) override
    {
        // call helper function: implementation of derivative depends on whether
        // ReactionRate::ddTFromStruct is defined
        _process_ddT(rop, kf, deltaT);
    }

    virtual void processRateConstants_ddP(double* rop,
                                          const double* kf,
                                          double deltaP) override
    {
        // call helper function: implementation of derivative depends on whether
        // ReactionData::perturbPressure is defined
        _process_ddP(rop, kf, deltaP);
    }

    virtual void processRateConstants_ddM(double* rop,
                                          const double* kf,
                                          double deltaM,
                                          bool overwrite=true) override
    {
        // call helper function: implementation of derivative depends on whether
        // ReactionRate::thirdBodyConcentration is defined
        _process_ddM(rop, kf, deltaM, overwrite);
    }

    virtual void update(double T) override {
        m_shared.update(T);
        _update();
    }

    virtual void update(double T, double extra) override {
        m_shared.update(T, extra);
        _update();
    }

    virtual void update(double T, const vector_fp& extra) override {
        m_shared.update(T, extra);
        _update();
    }

    virtual bool update(const ThermoPhase& phase, const Kinetics& kin) override {
        bool changed = m_shared.update(phase, kin);
        if (changed) {
            // call helper function only if needed: implementation depends on whether
            // ReactionRate::updateFromStruct is defined
            _update();
        }
        return changed;
    }

    virtual double evalSingle(ReactionRate& rate) override {
        RateType& R = static_cast<RateType&>(rate);
        _updateRate(R);
        return R.evalFromStruct(m_shared);
    }

protected:
    //! Helper function to process updates for rate types that implement the
    //! `updateFromStruct` method.
    template <typename T=RateType,
        typename std::enable_if<has_update<T>::value, bool>::type = true>
    void _update() {
        for (auto& rxn : m_rxn_rates) {
            rxn.second.updateFromStruct(m_shared);
        }
    }

    //! Helper function for rate types that do not implement `updateFromStruct`.
    //! Does nothing, but exists to allow generic implementations of update().
    template <typename T=RateType,
        typename std::enable_if<!has_update<T>::value, bool>::type = true>
    void _update() {
    }

    //! Helper function to update a single rate that has an `updateFromStruct` method`.
    template <typename T=RateType,
        typename std::enable_if<has_update<T>::value, bool>::type = true>
    void _updateRate(RateType& rate) {
        rate.updateFromStruct(m_shared);
    }

    //! Helper function for single rate that does not implement `updateFromStruct`.
    //! Exists to allow generic implementations of `evalSingle` and `ddTSingle`.
    template <typename T=RateType,
        typename std::enable_if<!has_update<T>::value, bool>::type = true>
    void _updateRate(RateType& rate) {
    }

    //! Helper function to process temperature derivatives for rate types that
    //! implement the `ddTScaledFromStruct` method.
    template <typename T=RateType,
        typename std::enable_if<has_ddT<T>::value, bool>::type = true>
    void _process_ddT(double* rop, const double* kf, double deltaT) {
        for (const auto& rxn : m_rxn_rates) {
            rop[rxn.first] *= rxn.second.ddTScaledFromStruct(m_shared);
        }
    }

    //! Helper function for rate types that do not implement `ddTScaledFromStruct`
    template <typename T=RateType,
        typename std::enable_if<!has_ddT<T>::value, bool>::type = true>
    void _process_ddT(double* rop, const double* kf, double deltaT) {

        // perturb conditions
        double dTinv = 1. / (m_shared.temperature * deltaT);
        m_shared.perturbTemperature(deltaT);
        _update();

        // apply numerical derivative
        for (auto& rxn : m_rxn_rates) {
            if (kf[rxn.first] != 0.) {
                double k1 = rxn.second.evalFromStruct(m_shared);
                rop[rxn.first] *= dTinv * (k1 / kf[rxn.first] - 1.);
            } // else not needed: derivative is already zero
        }

        // revert changes
        m_shared.restore();
        _update();
    }

    //! Helper function to process third-body derivatives for rate data that
    //! implement the `perturbThirdBodies` method.
    template <typename T=RateType, typename D=DataType,
        typename std::enable_if<has_ddM<D>::value, bool>::type = true>
    void _process_ddM(double* rop, const double* kf, double deltaM, bool overwrite) {
        double dMinv = 1. / deltaM;
        m_shared.perturbThirdBodies(deltaM);
        _update();

        for (auto& rxn : m_rxn_rates) {
            if (kf[rxn.first] != 0. && m_shared.conc_3b[rxn.first] > 0.) {
                double k1 = rxn.second.evalFromStruct(m_shared);
                rop[rxn.first] *= dMinv * (k1 / kf[rxn.first] - 1.);
                rop[rxn.first] /= m_shared.conc_3b[rxn.first];
            } else {
                rop[rxn.first] = 0.;
            }
        }

        // revert changes
        m_shared.restore();
        _update();
    }

    //! Helper function for rate data that do not implement `perturbThirdBodies`
    template <typename T=RateType, typename D=DataType,
        typename std::enable_if<!has_ddM<D>::value, bool>::type = true>
    void _process_ddM(double* rop, const double* kf, double deltaM, bool overwrite) {
        if (!overwrite) {
            // do not overwrite existing entries
            return;
        }
        for (const auto& rxn : m_rxn_rates) {
            rop[rxn.first] = 0.;
        }
    }

    //! Helper function to process pressure derivatives for rate data that
    //! implement the `perturbPressure` method.
    template <typename T=RateType, typename D=DataType,
        typename std::enable_if<has_ddP<D>::value, bool>::type = true>
    void _process_ddP(double* rop, const double* kf, double deltaP) {
        double dPinv = 1. / (m_shared.pressure * deltaP);
        m_shared.perturbPressure(deltaP);
        _update();

        for (auto& rxn : m_rxn_rates) {
            if (kf[rxn.first] != 0.) {
                double k1 = rxn.second.evalFromStruct(m_shared);
                rop[rxn.first] *= dPinv * (k1 / kf[rxn.first] - 1.);
            } // else not needed: derivative is already zero
        }

        // revert changes
        m_shared.restore();
        _update();
    }

    //! Helper function for rate data that do not implement `perturbPressure`
    template <typename T=RateType, typename D=DataType,
        typename std::enable_if<!has_ddP<D>::value, bool>::type = true>
    void _process_ddP(double* rop, const double* kf, double deltaP) {
        for (const auto& rxn : m_rxn_rates) {
            rop[rxn.first] = 0.;
        }
    }

    //! Vector of pairs of reaction rates indices and reaction rates
    std::vector<std::pair<size_t, RateType>> m_rxn_rates;
    std::map<size_t, size_t> m_indices; //! Mapping of indices
    DataType m_shared;
};

}

#endif
