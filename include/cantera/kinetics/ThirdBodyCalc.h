/**
 *  @file ThirdBodyCalc.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_THIRDBODYCALC_H
#define CT_THIRDBODYCALC_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

//! Calculate and apply third-body effects on reaction rates, including non-
//! unity third-body efficiencies.
//! Used by legacy reactions
class ThirdBodyCalc
{
public:
    void install(size_t rxnNumber, const std::map<size_t, double>& enhanced,
                 double dflt=1.0, size_t rxnIndex=npos) {
        m_reaction_index.push_back(rxnNumber);
        m_default.push_back(dflt);

        m_species.emplace_back();
        m_eff.emplace_back();
        for (const auto& eff : enhanced) {
            AssertTrace(eff.first != npos);
            m_species.back().push_back(eff.first);
            m_eff.back().push_back(eff.second - dflt);
        }
        if (rxnIndex == npos) {
            m_true_index.push_back(rxnNumber);
        } else {
            m_true_index.push_back(rxnIndex);
        }
    }

    void update(const vector_fp& conc, double ctot, double* work) {
        for (size_t i = 0; i < m_species.size(); i++) {
            double sum = 0.0;
            for (size_t j = 0; j < m_species[i].size(); j++) {
                sum += m_eff[i][j] * conc[m_species[i][j]];
            }
            work[i] = m_default[i] * ctot + sum;
        }
    }

    //! Update third-body concentrations in full vector
    void copy(const vector_fp& work, double* concm) {
        for (size_t i = 0; i < m_true_index.size(); i++) {
            concm[m_true_index[i]] = work[i];
        }
    }

    void multiply(double* output, const double* work) {
        for (size_t i = 0; i < m_reaction_index.size(); i++) {
            output[m_reaction_index[i]] *= work[i];
        }
    }

    size_t workSize() {
        return m_reaction_index.size();
    }

protected:
    //! Indices of reactions that use third-bodies within vector of concentrations
    std::vector<size_t> m_reaction_index;

    //! Actual index of reaction within the full reaction array
    std::vector<size_t> m_true_index;

    //! m_species[i][j] is the index of the j-th species in reaction i.
    std::vector<std::vector<size_t> > m_species;

    //! m_eff[i][j] is the efficiency of the j-th species in reaction i.
    std::vector<vector_fp> m_eff;

    //! The default efficiency for each reaction
    vector_fp m_default;
};


//! Calculate and apply third-body effects on reaction rates, including non-
//! unity third-body efficiencies.
class ThirdBodyCalc3
{
public:
    //! Install reaction that uses third-body effects in ThirdBodyCalc3 manager
    void install(size_t rxnNumber, const std::map<size_t, double>& efficiencies,
                  double default_efficiency, bool mass_action) {
        m_reaction_index.push_back(rxnNumber);
        m_default.push_back(default_efficiency);

        if (mass_action) {
            m_mass_action.push_back(true);
            m_mass_action_index.push_back(m_reaction_index.size() - 1);
        } else {
            m_mass_action.push_back(false);
        }

        m_species.emplace_back();
        m_eff.emplace_back();
        for (const auto& eff : efficiencies) {
            AssertTrace(eff.first != npos);
            m_species.back().push_back(eff.first);
            m_eff.back().push_back(eff.second - default_efficiency);
            m_efficiencyList.emplace_back(
                rxnNumber, eff.first, eff.second - default_efficiency);
        }
    }

    void resizeCoeffs(size_t nSpc, size_t nRxn)
    {
        // Sparse Efficiency coefficient matrix
        Eigen::SparseMatrix<double> efficiencies;
        efficiencies.setZero();
        efficiencies.resize(nRxn, nSpc);
        efficiencies.reserve(m_efficiencyList.size());
        efficiencies.setFromTriplets(
            m_efficiencyList.begin(), m_efficiencyList.end());

        // Jacobian matrix multipliers
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(nRxn * nSpc);
        for (size_t i = 0; i < m_default.size(); i++) {
            if (m_default[i] != 0) {
                for (size_t j = 0; j < nSpc; j++) {
                    triplets.emplace_back(m_reaction_index[i], j, m_default[i]);
                }
            }
        }
        Eigen::SparseMatrix<double> defaults(nRxn, nSpc);
        defaults.reserve(triplets.size());
        defaults.setFromTriplets(triplets.begin(), triplets.end());
        m_multipliers = efficiencies + defaults;
    }

    //! Update third-body concentrations in full vector
    void update(const vector_fp& conc, double ctot, double* concm) const {
        for (size_t i = 0; i < m_reaction_index.size(); i++) {
            double sum = 0.0;
            for (size_t j = 0; j < m_species[i].size(); j++) {
                sum += m_eff[i][j] * conc[m_species[i][j]];
            }
            concm[m_reaction_index[i]] = m_default[i] * ctot + sum;
        }
    }

    //! Multiply output with effective third-body concentration
    void multiply(double* output, const double* concm) {
        for (size_t i = 0; i < m_mass_action_index.size(); i++) {
            size_t ix = m_reaction_index[m_mass_action_index[i]];
            output[ix] *= concm[ix];
        }
    }

    //! Calculate derivatives with respect to species concentrations.
    /*!
     *  @param product   Product of law of mass action and rate terms.
     */
    Eigen::SparseMatrix<double> jacobian(const double* product)
    {
        Eigen::Map<const Eigen::VectorXd> mapped(product, m_multipliers.rows());
        return mapped.asDiagonal() * m_multipliers;
    }

    //! Scale input by reaction order with respect to third-body concentrations
    void scaleOrder(const double* in, double* out) const
    {
        for (const auto& ix : m_reaction_index) {
            if (m_mass_action[ix]) {
                out[ix] = in[ix];
            }
        }
    }

protected:
    //! Indices of reactions that use third-bodies within vector of concentrations
    std::vector<size_t> m_reaction_index;

    //! Vector of flags indicating whether reactions consider third-body effects
    std::vector<bool> m_mass_action;

    //! Indices within m_reaction_index of reactions that consider third-body effects
    //! in the law of mass action
    std::vector<size_t> m_mass_action_index;

    //! m_species[i][j] is the index of the j-th species in reaction i.
    std::vector<std::vector<size_t> > m_species;

    //! m_eff[i][j] is the efficiency of the j-th species in reaction i.
    std::vector<vector_fp> m_eff;

    //! The default efficiency for each reaction
    vector_fp m_default;

    //! Sparse efficiency matrix (compensated for defaults)
    std::vector<Eigen::Triplet<double>> m_efficiencyList;

    //! Sparse Jacobian multiplier matrix
    Eigen::SparseMatrix<double> m_multipliers;
};


}

#endif
