/**
 *  @file FalloffMgr.h
 *
 *  @deprecated  Deprecated in Cantera 2.6 and removed thereafter. Replaced by
 *      FalloffRate objects managed by MultiRate evaluators.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_FALLOFFMGR_H
#define CT_FALLOFFMGR_H

#include "Falloff.h"
#include "cantera/base/global.h"

namespace Cantera
{


/**
 *  A falloff manager that implements any set of falloff functions.
 *  @ingroup falloffGroup
 *
 *  @deprecated  Deprecated in Cantera 2.6 and removed thereafter. Replaced by
 *      FalloffRate objects managed by MultiRate evaluators.
 */
class FalloffMgr
{
public:
    //! Constructor.
    FalloffMgr() :
        m_worksize(0) {
    }

    //! Install a new falloff function calculator.
    /*
     * @param rxn Index of the falloff reaction. This will be used to
     *     determine which array entry is modified in method pr_to_falloff.
     * @param reactionType Either `FALLOFF_RXN` or `CHEMACT_RXN`
     * @param f The falloff function.
     *
     * @deprecated To be removed after Cantera 2.6.
     */
    void install(size_t rxn, int reactionType, shared_ptr<Falloff> f) {
        warn_deprecated("FalloffMgr::install()",
            "To be removed after Cantera 2.6. Specify reaction type using "
            "string instead.");

        m_rxn.push_back(rxn);
        m_offset.push_back(m_worksize);
        m_worksize += f->workSize();
        m_falloff.push_back(f);
        m_isfalloff.push_back(reactionType == FALLOFF_RXN);
        m_indices[rxn] = m_falloff.size()-1;
    }

    //! Install a new falloff function calculator.
    /*
     * @param rxn Index of the falloff reaction. This will be used to
     *     determine which array entry is modified in method pr_to_falloff.
     * @param type Reaction type identifier.
     * @param f The falloff function.
     */
     void install(size_t rxn, std::string type, shared_ptr<Falloff> f) {
        m_rxn.push_back(rxn);
        m_offset.push_back(m_worksize);
        m_worksize += f->workSize();
        m_falloff.push_back(f);
        m_isfalloff.push_back(type == "falloff-legacy");
        m_indices[rxn] = m_falloff.size()-1;
    }

    /*!
     * Replace an existing falloff function calculator
     *
     * @param rxn   External reaction index
     * @param f     New falloff function, of the same kind as the existing one
     */
    void replace(size_t rxn, shared_ptr<Falloff> f) {
        m_falloff[m_indices[rxn]] = f;
    }

    //! Size of the work array required to store intermediate results.
    size_t workSize() {
        return m_worksize;
    }

    /**
     * Update the cached temperature-dependent intermediate
     * results for all installed falloff functions.
     * @param t Temperature [K].
     * @param work Work array. Must be dimensioned at least workSize().
     */
    void updateTemp(doublereal t, doublereal* work) {
        for (size_t i = 0; i < m_rxn.size(); i++) {
            m_falloff[i]->updateTemp(t, work + m_offset[i]);
        }
    }

    /**
     * Given a vector of reduced pressures for each falloff reaction,
     * replace each entry by the value of the falloff function.
     */
    void pr_to_falloff(doublereal* values, const doublereal* work) {
        for (size_t i = 0; i < m_rxn.size(); i++) {
            double pr = values[m_rxn[i]];
            if (m_isfalloff[i]) {
                // Pr / (1 + Pr) * F
                values[m_rxn[i]] *=
                    m_falloff[i]->F(pr, work + m_offset[i]) / (1.0 + pr);
            } else {
                // 1 / (1 + Pr) * F
                values[m_rxn[i]] =
                    m_falloff[i]->F(pr, work + m_offset[i]) / (1.0 + pr);
            }
        }
    }

protected:
    std::vector<size_t> m_rxn;
    std::vector<shared_ptr<Falloff> > m_falloff;
    vector_int m_loc;
    std::vector<vector_fp::difference_type> m_offset;
    size_t m_worksize;

    //! Distinguish between falloff and chemically activated reactions
    std::vector<bool> m_isfalloff;

    //! map of external reaction index to local index
    std::map<size_t, size_t> m_indices;
};
}

#endif
