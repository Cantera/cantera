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
 * the actual implementation for this capability is done in the MultiRate class.
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
    //! @param nSpecies  number of species
    //! @param nReactions  number of reactions
    //! @param nPhases  number of phases
    virtual void resize(size_t nSpecies, size_t nReactions, size_t nPhases) = 0;

    //! Evaluate all rate constants handled by the evaluator
    //! @param kf  array of rate constants
    virtual void getRateConstants(double* kf) = 0;

    //! Evaluate all rate constant temperature derivatives handled by the evaluator;
    //! which are multiplied with the array of rate-of-progress variables.
    //! Depending on the implementation of a rate object, either an exact derivative or
    //! a numerical approximation may be used.
    //! @param[in,out] rop  array of rop, which is modified by the method;
    //!     contains rop on input, and d(rop)/dT on output
    //! @param kf  array of forward rate constants (numerical derivative only)
    //! @param deltaT  relative temperature perturbation (numerical derivative only)
    virtual void processRateConstants_ddT(double* rop,
                                          const double* kf,
                                          double deltaT) = 0;

    //! Evaluate all rate constant pressure derivatives handled by the evaluator;
    //! which are multiplied with the array of rate-of-progress variables.
    //! @param[in,out] rop  array of rop, which is modified by the method;
    //!     contains rop on input, and d(rop)/dP on output
    //! @param kf  array of forward rate constants
    //! @param deltaP  relative pressure perturbation
    virtual void processRateConstants_ddP(double* rop,
                                          const double* kf,
                                          double deltaP) = 0;

    //! Evaluate all rate constant third-body derivatives handled by the evaluator;
    //! which are multiplied with the array of rate-of-progress variables.
    //! @param[in,out] rop  array of rop, which is modified by the method;
    //!     contains rop on input, and d(rop)/dM on output
    //! @param kf  array of forward rate constants
    //! @param deltaM  relative perturbation of third-body concentrations
    //! @param overwrite  if `true`, rop entries not affected by M are set to zero
    virtual void processRateConstants_ddM(double* rop,
                                          const double* kf,
                                          double deltaM,
                                          bool overwrite=true) = 0;

    //! Update common reaction rate data based on temperature.
    //! Only used in conjunction with evalSingle and ReactionRate::eval
    //! @param T  temperature [K]
    virtual void update(double T) = 0;

    //! Update common reaction rate data based on temperature and extra parameter.
    //! Only used in conjunction with evalSingle and ReactionRate::eval
    //! @param T  temperature [K]
    //! @param extra  extra parameter (depends on parameterization)
    virtual void update(double T, double extra) = 0;

    //! Update common reaction rate data based on temperature and extra parameter.
    //! Only used in conjunction with evalSingle and ReactionRate::eval
    //! @param T  temperature [K]
    //! @param extra  extra vector parameter (depends on parameterization)
    //! @warning  This method is an experimental part of the %Cantera API and
    //!     may be changed or removed without notice.
    virtual void update(double T, const vector_fp& extra) = 0;

    //! Update data common to reaction rates of a specific type.
    //! This update mechanism is used by Kinetics reaction rate evaluators.
    //! @param phase  object representing reacting phase
    //! @param kin  object representing kinetics
    //! @returns  flag indicating whether reaction rates need to be re-evaluated
    virtual bool update(const ThermoPhase& phase, const Kinetics& kin) = 0;

    //! Get the rate for a single reaction. Used to implement ReactionRate::eval,
    //! which allows for the evaluation of a reaction rate expression outside of
    //! Kinetics reaction rate evaluators. Mainly used for testing purposes.
    virtual double evalSingle(ReactionRate& rate) = 0;
};

} // end namespace Cantera

#endif
