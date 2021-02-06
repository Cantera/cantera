/**
 *  @file ThirdBodyCalc.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_THIRDBODYCALC_H
#define CT_THIRDBODYCALC_H

#include "cantera/base/ct_defs.h"
#include <cassert>

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

    void multiply(double* output, const double* work) {
        for (size_t i = 0; i < m_reaction_index.size(); i++) {
            output[m_reaction_index[i]] *= work[i];
        }
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
};

}

#endif
