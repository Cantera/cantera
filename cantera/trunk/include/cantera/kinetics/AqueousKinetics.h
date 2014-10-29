/**
 * @file AqueousKinetics.h
 *
 * @ingroup chemkinetics
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_AQUEOUSKINETICS_H
#define CT_AQUEOUSKINETICS_H

#include "Kinetics.h"
#include "ReactionStoichMgr.h"
#include "RateCoeffMgr.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

// forward references

class ReactionData;

/**
 * Kinetics manager for elementary aqueous-phase chemistry. This
 * kinetics manager implements standard mass-action reaction rate
 * expressions for liquids
 *
 *
 *   Concentration
 *
 * @ingroup kinetics
 * @deprecated Not actually implemented
 */
class AqueousKinetics : public Kinetics
{

public:

    //! @name Constructors
    //! @{

    /// Constructor. Creates an empty reaction mechanism.
    AqueousKinetics(thermo_t* thermo = 0);

    AqueousKinetics(const AqueousKinetics& right);

    AqueousKinetics& operator=(const AqueousKinetics& right);

    //! Duplication routine for objects which inherit from Kinetics
    /*!
     *  This virtual routine can be used to duplicate %Kinetics objects
     *  inherited from %Kinetics even if the application only has
     *  a pointer to %Kinetics to work with.
     *
     *  These routines are basically wrappers around the derived copy  constructor.
     *
     * @param  tpVector Vector of shallow pointers to ThermoPhase objects. this is the
     *                  m_thermo vector within this object
     */
    virtual Kinetics* duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const;
    //@}

    virtual int type() const {
        return cAqueousKinetics;
    }

    virtual doublereal reactantStoichCoeff(size_t k, size_t i) const {
        return getValue(m_rrxn[k], i, 0.0);
    }

    virtual doublereal productStoichCoeff(size_t k, size_t i) const {
        return getValue(m_prxn[k], i, 0.0);
    }

    //! @name Reaction Rates Of Progress
    //@{
    virtual void getEquilibriumConstants(doublereal* kc);

    virtual void getDeltaGibbs(doublereal* deltaG);
    virtual void getDeltaEnthalpy(doublereal* deltaH);
    virtual void getDeltaEntropy(doublereal* deltaS);

    virtual void getDeltaSSGibbs(doublereal* deltaG);
    virtual void getDeltaSSEnthalpy(doublereal* deltaH);
    virtual void getDeltaSSEntropy(doublereal* deltaS);

    //! @}
    //! @name Reaction Mechanism Informational Query Routines
    //! @{

    virtual bool isReversible(size_t i) {
        if (std::find(m_revindex.begin(), m_revindex.end(), i)
                < m_revindex.end()) {
            return true;
        } else {
            return false;
        }
    }

    virtual void getFwdRateConstants(doublereal* kfwd);
    virtual void getRevRateConstants(doublereal* krev,
                                     bool doIrreversible = false);

    //! @}
    //! @name Reaction Mechanism Setup Routines
    //! @{

    virtual void init();
    virtual void addReaction(ReactionData& r);
    virtual void finalize();
    virtual bool ready() const;

    virtual void update_T();

    virtual void update_C();

    void updateROP();

    /*!
     * Update temperature-dependent portions of reaction rates and
     * falloff functions.
     */
    void _update_rates_T();

    /*!
     * Update properties that depend on concentrations. Currently only
     * the enhanced collision partner concentrations are updated here.
     */
    void _update_rates_C();

    //@}

protected:

    size_t m_nfall;

    Rate1<Arrhenius>                    m_rates;

    std::vector<size_t> m_irrev;

    std::vector<size_t> m_fwdOrder;

    size_t m_nirrev;
    size_t m_nrev;

    /**
     * Difference between the input global reactants order
     * and the input global products order. Changed to a double
     * to account for the fact that we can have real-valued
     * stoichiometries.
     */
    vector_fp  m_dn;
    std::vector<size_t> m_revindex;

    vector_fp m_conc;
    vector_fp m_grt;

    //! @name Aqueous kinetics data
    //!@{
    bool m_ROP_ok;

    doublereal m_temp;
    vector_fp  m_rfn;

    vector_fp m_rkcn;
    //!@}

private:
    void addElementaryReaction(ReactionData& r);

    void installReagents(const ReactionData& r);

    /**
     * Update the equilibrium constants in molar units.
     */
    void updateKc();

    bool m_finalized;
};
}

#endif
