/**
 *  @file ThirdBodyCalc.h
 */

#ifndef CT_THIRDBODYCALC_H
#define CT_THIRDBODYCALC_H

#include "cantera/base/utilities.h"
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

        m_species.push_back(std::vector<size_t>());
        m_eff.push_back(vector_fp());
        for (std::map<size_t, double>::const_iterator iter = enhanced.begin();
             iter != enhanced.end();
             ++iter)
        {
            assert(iter->first != npos);
            m_species.back().push_back(iter->first);
            m_eff.back().push_back(iter->second - dflt);
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
        scatter_mult(work, work + m_reaction_index.size(),
                     output, m_reaction_index.begin());
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
