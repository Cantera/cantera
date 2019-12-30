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

namespace Cantera
{

class ElementaryReaction;

//! Partial specialization of Kinetics for chemistry in a single bulk phase
class BulkKinetics : public Kinetics
{
public:
    BulkKinetics(thermo_t* thermo = 0);

    virtual bool isReversible(size_t i);

    virtual void getDeltaGibbs(doublereal* deltaG);
    virtual void getDeltaEnthalpy(doublereal* deltaH);
    virtual void getDeltaEntropy(doublereal* deltaS);

    virtual void getDeltaSSGibbs(doublereal* deltaG);
    virtual void getDeltaSSEnthalpy(doublereal* deltaH);
    virtual void getDeltaSSEntropy(doublereal* deltaS);

    virtual void getRevRateConstants(doublereal* krev,
                                     bool doIrreversible = false);

    virtual bool addReaction(shared_ptr<Reaction> r);
    virtual void resizeSpecies();

    virtual void setMultiplier(size_t i, double f);
    virtual void invalidateCache();

protected:
    virtual void addElementaryReaction(ElementaryReaction& r);
    virtual void modifyElementaryReaction(size_t i, ElementaryReaction& rNew);

    Rate1<Arrhenius> m_rates;
    std::vector<size_t> m_revindex; //!< Indices of reversible reactions
    std::vector<size_t> m_irrev; //!< Indices of irreversible reactions

    //! Difference between the global reactants order and the global products
    //! order. Of type "double" to account for the fact that we can have real-
    //! valued stoichiometries.
    vector_fp m_dn;

    vector_fp m_conc;
    vector_fp m_grt;

    bool m_ROP_ok;
    doublereal m_temp;
};

}

#endif
