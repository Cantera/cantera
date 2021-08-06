/**
 *  @file ThirdBodyCalc.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_THIRDBODYCALC_H
#define CT_THIRDBODYCALC_H

#include "cantera/base/ct_defs.h"
#include <cassert>
#include <numeric>

namespace Cantera
{

//! Calculate and apply third-body effects on reaction rates, including non-
//! unity third-body efficiencies.
class ThirdBodyCalc
{
public:
    void install(size_t rxnNumber, const std::map<size_t, double>& enhanced,
                 double dflt=1.0) {
        m_reaction_index.push_back(rxnNumber);
        m_default.push_back(dflt);

        m_species.emplace_back();
        m_eff.emplace_back();
        for (const auto& eff : enhanced) {
            assert(eff.first != npos);
            m_species.back().push_back(eff.first);
            m_eff.back().push_back(eff.second - dflt);
            m_efficiencyList.emplace_back(rxnNumber, eff.first, eff.second - dflt);
        }
    }

    void finalizeSetup(size_t nSpc, size_t nRxn)
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

    void update(const vector_fp& conc, double ctot, double* work) {
        for (size_t i = 0; i < m_species.size(); i++) {
            double sum = 0.0;
            for (size_t j = 0; j < m_species[i].size(); j++) {
                sum += m_eff[i][j] * conc[m_species[i][j]];
            }
            work[i] = m_default[i] * ctot + sum;
        }
    }

    void multiply(double* output, const double* work) {
        // @TODO: Tests reveal that multiply is usually called exactly once after
        //      each update.
        for (size_t i = 0; i < m_reaction_index.size(); i++) {
            output[m_reaction_index[i]] *= work[i];
        }
    }

    //! Calculate derivatives with respect to species concentrations.
    /*!
     *  @param product   Product of law of mass action and rate terms.
     */
    Eigen::SparseMatrix<double> speciesDerivatives(const double* product)
    {
        ConstMappedVector mapped(product, m_multipliers.rows());
        return mapped.asDiagonal() * m_multipliers;
    }

    size_t workSize() {
        return m_reaction_index.size();
    }

protected:
    //! Indices of third-body reactions within the full reaction array
    std::vector<size_t> m_reaction_index;

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
