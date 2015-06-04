/**
 * @file AqueousKinetics.h
 * @ingroup chemkinetics
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_AQUEOUSKINETICS_H
#define CT_AQUEOUSKINETICS_H

#include "BulkKinetics.h"

namespace Cantera
{

/**
 * Kinetics manager for elementary aqueous-phase chemistry. This kinetics
 * manager implements standard mass-action reaction rate expressions for liquids
 *
 * @ingroup kinetics
 */
class AqueousKinetics : public BulkKinetics
{
public:
    /// Constructor. Creates an empty reaction mechanism.
    AqueousKinetics(thermo_t* thermo = 0);

    //! Duplication routine for objects which inherit from Kinetics
    /*!
     *  This virtual routine can be used to duplicate Kinetics objects
     *  inherited from Kinetics even if the application only has
     *  a pointer to Kinetics to work with.
     *
     *  These routines are basically wrappers around the derived copy  constructor.
     *
     * @param  tpVector Vector of shallow pointers to ThermoPhase objects. this is the
     *                  m_thermo vector within this object
     */
    virtual Kinetics* duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const;

    virtual int type() const {
        return cAqueousKinetics;
    }

    virtual void getEquilibriumConstants(doublereal* kc);
    virtual void getFwdRateConstants(doublereal* kfwd);
    void updateROP();

    //! Update temperature-dependent portions of reaction rates
    void _update_rates_T();

    //! Update properties that depend on concentrations.
    void _update_rates_C();

    //! Update the equilibrium constants in molar units.
    void updateKc();

    virtual void addReaction(ReactionData& r);
    virtual bool addReaction(shared_ptr<Reaction> r);
    virtual void modifyReaction(size_t i, shared_ptr<Reaction> rNew);
};
}

#endif
