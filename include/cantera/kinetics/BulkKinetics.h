/**
 * @file BulkKinetics.h
 * @ingroup chemkinetics
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_BULKKINETICS_H
#define CT_BULKKINETICS_H

#include "Kinetics.h"
#include "RateCoeffMgr.h"
#include "ThirdBodyCalc.h"
#include "cantera/kinetics/MultiRate.h"

namespace Cantera
{

class ElementaryReaction2;

//! Partial specialization of Kinetics for chemistry in a single bulk phase
class BulkKinetics : public Kinetics
{
public:
    BulkKinetics(ThermoPhase* thermo = 0);

    virtual void resizeReactions();

    virtual bool isReversible(size_t i);

    virtual void getDeltaGibbs(doublereal* deltaG);
    virtual void getDeltaEnthalpy(doublereal* deltaH);
    virtual void getDeltaEntropy(doublereal* deltaS);

    virtual void getDeltaSSGibbs(doublereal* deltaG);
    virtual void getDeltaSSEnthalpy(doublereal* deltaH);
    virtual void getDeltaSSEntropy(doublereal* deltaS);

    virtual void getRevRateConstants(double* krev,
                                     bool doIrreversible = false);

    virtual bool addReaction(shared_ptr<Reaction> r, bool resize=true);
    virtual void modifyReaction(size_t i, shared_ptr<Reaction> rNew);

    virtual void resizeSpecies();

    virtual void setMultiplier(size_t i, double f);
    virtual void invalidateCache();

    void addThirdBody(shared_ptr<Reaction> r);

protected:
    virtual void addElementaryReaction(ElementaryReaction2& r);
    virtual void modifyElementaryReaction(size_t i, ElementaryReaction2& rNew);

    //! Vector of rate handlers
    std::vector<unique_ptr<MultiRateBase>> m_bulk_rates;
    std::map<std::string, size_t> m_bulk_types; //!< Mapping of rate handlers

    Rate1<Arrhenius2> m_rates; //!< @deprecated (legacy only)
    std::vector<size_t> m_revindex; //!< Indices of reversible reactions
    std::vector<size_t> m_irrev; //!< Indices of irreversible reactions

    //! Difference between the global reactants order and the global products
    //! order. Of type "double" to account for the fact that we can have real-
    //! valued stoichiometries.
    vector_fp m_dn;

    ThirdBodyCalc3 m_multi_concm; //!< used with MultiRate evaluator

    //! Third body concentrations
    vector_fp m_concm;

    //! Activity concentrations, as calculated by ThermoPhase::getActivityConcentrations
    vector_fp m_act_conc;

    //! Physical concentrations, as calculated by ThermoPhase::getConcentrations
    vector_fp m_phys_conc;

    vector_fp m_grt;

    bool m_ROP_ok;
    doublereal m_temp;
};

}

#endif
