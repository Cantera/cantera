/**
 * @file MultiRateBase.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MULTIRATEBASE_H
#define CT_MULTIRATEBASE_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

class ReactionRate;
class ThermoPhase;
class Kinetics;

//! An abstract base class for evaluating all reactions of a particular type.
/**
 * Because this class has no template parameters, the Kinetics object can store all of
 * these rate coefficient evaluators as a `vector<shared_ptr<MultiRateBase>>`. All of
 * the actual implementation for this capability is done in the MultiBulkRate class.
 *
 * @todo At the moment, implemented methods are specific to BulkKinetics,
 *     which can be updated using information of a single ThermoPhase.
 *     InterfaceKinetics will require access to an entire Kinetics object
 *     or the underlying `vector<ThermoPhase*>` vector (e.g. `m_thermo`).
 */
class MultiRateBase
{
public:
    virtual ~MultiRateBase() {}

    //! Identifier of reaction rate type
    virtual std::string type() = 0;

    //! Add reaction rate object to the evaluator
    //! @param rxn_index  index of reaction
    //! @param rate  reaction rate object
    virtual void add(size_t rxn_index, ReactionRate& rate) = 0;

    //! Replace reaction rate object handled by the evaluator
    //! @param rxn_index  index of reaction
    //! @param rate  reaction rate object
    virtual bool replace(size_t rxn_index, ReactionRate& rate) = 0;

    //! Update number of species and reactions
    //! @param n_species  number of species
    //! @param n_reactions  number of reactions
    virtual void resize(size_t n_species, size_t n_reactions) = 0;

    //! Evaluate all rate constants handled by the evaluator
    //! @param kf  array of rate constants
    virtual void getRateConstants(double* kf) = 0;

    //! Update common reaction rate data based on temperature.
    //! Only used in conjunction with evalSingle and ReactionRate::eval
    //! @param T  temperature [K]
    virtual void update(double T) = 0;

    //! Update common reaction rate data based on temperature and extra parameter.
    //! Only used in conjunction with evalSingle and ReactionRate::eval
    //! @param T  temperature [K]
    //! @param extra  extra parameter (depends on parameterization)
    virtual void update(double T, double extra) = 0;

    //! Update data common to reaction rates of a specific type.
    //! This update mechanism is used by Kinetics reaction rate evaluators.
    //! @param bulk  object representing bulk phase
    //! @param kin  object representing kinetics
    //! @returns flag indicating reaction rates need to be re-evaluated
    virtual bool update(const ThermoPhase& bulk, const Kinetics& kin) = 0;

    //! Get the rate for a single reaction. Used to implement ReactionRate::eval,
    //! which allows for the evaluation of a reaction rate expression outside of
    //! Kinetics reaction rate evaluators. Mainly used for testing purposes.
    virtual double evalSingle(ReactionRate& rate) = 0;
};

} // end namespace Cantera

#endif
