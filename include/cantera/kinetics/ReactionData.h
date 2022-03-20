/**
 * @file ReactionData.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTIONDATA_H
#define CT_REACTIONDATA_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

class ThermoPhase;
class Kinetics;


//! Data container holding shared data used for ReactionRate calculation
/**
 * The base class defines variables and methods used by all specializations.
 */
struct ReactionData
{
    ReactionData() : temperature(1.), logT(0.), recipT(1.), m_temperature_buf(-1.) {}

    //! Update data container based on temperature *T*
    /**
     * Only used in conjunction with MultiRateBase::evalSingle / ReactionRate::eval.
     * This method allows for testing of a reaction rate expression outside of
     * Kinetics reaction rate evaluators.
     */
    virtual void update(double T) {
        temperature = T;
        logT = std::log(T);
        recipT = 1./T;
    }

    //! Update data container based on temperature *T* and an *extra* parameter
    /**
     * Only used in conjunction with MultiRateBase::evalSingle / ReactionRate::eval.
     * This method allows for testing of a reaction rate expression outside of
     * Kinetics reaction rate evaluators.
     */
    virtual void update(double T, double extra) {
        throw NotImplementedError("ReactionData::update",
            "ReactionData type does not use extra scalar argument.");
    }

    //! Update data container based on temperature *T* and a vector parameter *extra*
    /**
     * Only used in conjunction with MultiRateBase::evalSingle / ReactionRate::eval.
     * This method allows for testing of a reaction rate expression outside of
     * Kinetics reaction rate evaluators.
     *
     * @warning  This method is an experimental part of the %Cantera API and
     *      may be changed or removed without notice.
     */
    virtual void update(double T, const vector_fp& extra) {
        throw NotImplementedError("ReactionData::update",
            "ReactionData type does not use extra vector argument.");
    }

    //! Update data container based on thermodynamic phase state
    /**
     * This update mechanism is used by Kinetics reaction rate evaluators.
     * @returns  A boolean element indicating whether the `evalFromStruct` method
     *      needs to be called (assuming previously-calculated values were cached)
     */
    virtual bool update(const ThermoPhase& phase, const Kinetics& kin) = 0;

    //! Perturb temperature of data container
    /**
     * The method is used for the evaluation of numerical derivatives.
     * @param  deltaT  relative temperature perturbation
     */
    void perturbTemperature(double deltaT) {
        if (m_temperature_buf > 0.) {
            throw CanteraError("ReactionData::perturbTemperature",
                "Cannot apply another perturbation as state is already perturbed.");
        }
        m_temperature_buf = temperature;
        ReactionData::update(temperature * (1. + deltaT));
    }

    //! Restore data container after a perturbation
    virtual void restore() {
        // only restore if there is a valid buffered value
        if (m_temperature_buf < 0.) {
            return;
        }
        ReactionData::update(m_temperature_buf);
        m_temperature_buf = -1.;
    }

    //! Update number of species, reactions and phases
    virtual void resize(size_t nSpecies, size_t nReactions, size_t nPhases) {}

    //! Force shared data and reaction rates to be updated next time. This is called by
    //! functions that change quantities affecting rate calculations that are normally
    //! assumed to be constant, like the reaction rate parameters or the number of
    //! reactions.
    virtual void invalidateCache() {
        temperature = NAN;
    }

    double temperature; //!< temperature
    double logT; //!< logarithm of temperature
    double recipT; //!< inverse of temperature

protected:
    double m_temperature_buf; //!< buffered temperature
};

}

#endif
