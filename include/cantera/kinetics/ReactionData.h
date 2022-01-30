/**
 * @file ReactionData.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTIONDATA_H
#define CT_REACTIONDATA_H

#include "cantera/base/ct_defs.h"

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
    virtual void update(double T, double extra);

    //! Update data container based on thermodynamic phase state
    /**
     * This update mechanism is used by Kinetics reaction rate evaluators.
     * @returns  A boolean element indicating whether the `evalFromStruct` method
     *      needs to be called (assuming previously-calculated values were cached)
     *
     * @todo  Remove Kinetics argument
     */
    virtual bool update(const ThermoPhase& phase, const Kinetics& kin) = 0;

    //! Perturb temperature of data container
    /**
     * The method is used for the evaluation of numerical derivatives.
     * @param  deltaT  relative temperature perturbation
     */
    void perturbTemperature(double deltaT);

    //! Restore data container after a perturbation
    virtual void restore();

    //! Update number of species and reactions
    virtual void resize(size_t n_species, size_t n_reactions) {}

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


//! Data container holding shared data specific to ArrheniusRate
/**
 * The data container `ArrheniusData` holds precalculated data common to
 * all `ArrheniusRate` objects.
 */
struct ArrheniusData : public ReactionData
{
    virtual bool update(const ThermoPhase& phase, const Kinetics& kin);
    using ReactionData::update;
};


//! Data container holding shared data specific to TwoTempPlasmaRate
/**
 * The data container `TwoTempPlasmaData` holds precalculated data common to
 * all `TwoTempPlasmaRate` objects.
 */
struct TwoTempPlasmaData : public ReactionData
{
    TwoTempPlasmaData() : electronTemp(1.), logTe(0.), recipTe(1.) {}

    virtual bool update(const ThermoPhase& phase, const Kinetics& kin) override;

    virtual void update(double T) override;

    virtual void update(double T, double Te) override;

    virtual void updateTe(double Te);

    virtual void invalidateCache() override {
        ReactionData::invalidateCache();
        electronTemp = NAN;
    }

    double electronTemp; //!< electron temperature
    double logTe; //!< logarithm of electron temperature
    double recipTe; //!< inverse of electron temperature
};


//! Data container holding shared data specific to BlowersMaselRate
/**
 * The data container `BlowersMaselData` holds precalculated data common to
 * all `BlowersMaselRate` objects.
 */
struct BlowersMaselData : public ReactionData
{
    BlowersMaselData();

    virtual bool update(const ThermoPhase& phase, const Kinetics& kin) override;

    virtual void update(double T) override;

    virtual void update(double T, double deltaH) override;

    virtual void resize(size_t n_species, size_t n_reactions) override {
        grt.resize(n_species, 0.);
        ready = true;
    }

    bool ready; //!< boolean indicating whether vectors are accessible
    double density; //!< used to determine if updates are needed
    double dH_direct; //!< enthalpy change for each reaction (for testing)
    vector_fp grt; //!< partial molar enthalpies

protected:
    int m_state_mf_number; //!< integer that is incremented when composition changes
};


//! Data container holding shared data specific to Falloff rates
/**
 * The data container `FalloffData` holds precalculated data common to
 * all Falloff related reaction rate classes.
 */
struct FalloffData : public ReactionData
{
    FalloffData();

    virtual bool update(const ThermoPhase& phase, const Kinetics& kin) override;

    virtual void update(double T) override;

    virtual void update(double T, double M) override;

    //! Perturb third-body concentration vector of data container
    /**
     * The method is used for the evaluation of numerical derivatives.
     * @param deltaM  relative third-body perturbation
     */
    void perturbThirdBodies(double deltaM);

    virtual void restore() override;

    virtual void resize(size_t n_species, size_t n_reactions) override {
        conc_3b.resize(n_reactions, NAN);
        m_conc_3b_buf.resize(n_reactions, NAN);
        ready = true;
    }

    virtual void invalidateCache() override {
        ReactionData::invalidateCache();
        molar_density = NAN;
    }

    bool ready; //!< boolean indicating whether vectors are accessible
    double molar_density; //!< used to determine if updates are needed
    vector_fp conc_3b; //!< vector of effective third-body concentrations

protected:
    int m_state_mf_number; //!< integer that is incremented when composition changes
    bool m_perturbed; //!< boolean indicating whether 3-rd body values are perturbed
    vector_fp m_conc_3b_buf; //!< buffered third-body concentrations
};


//! Data container holding shared data specific to PlogRate
/**
 * The data container `PlogData` holds precalculated data common to
 * all `PlogRate` objects.
 */
struct PlogData : public ReactionData
{
    PlogData() : pressure(NAN), logP(0.), m_pressure_buf(-1.) {}

    virtual void update(double T) override;

    virtual void update(double T, double P) override {
        ReactionData::update(T);
        pressure = P;
        logP = std::log(P);
    }

    virtual bool update(const ThermoPhase& phase, const Kinetics& kin) override;

    //! Perturb pressure of data container
    /**
     * The method is used for the evaluation of numerical derivatives.
     * @param  deltaP  relative pressure perturbation
     */
    void perturbPressure(double deltaP);

    virtual void restore() override;

    virtual void invalidateCache() override {
        ReactionData::invalidateCache();
        pressure = NAN;
    }

    double pressure; //!< pressure
    double logP; //!< logarithm of pressure

protected:
    double m_pressure_buf; //!< buffered pressure
};


//! Data container holding shared data specific to ChebyshevRate
/**
 * The data container `ChebyshevData` holds precalculated data common to
 * all `ChebyshevRate` objects.
 */
struct ChebyshevData : public ReactionData
{
    ChebyshevData() : pressure(NAN), log10P(0.), m_pressure_buf(-1.) {}

    virtual void update(double T) override;

    virtual void update(double T, double P) override {
        ReactionData::update(T);
        pressure = P;
        log10P = std::log10(P);
    }

    virtual bool update(const ThermoPhase& phase, const Kinetics& kin) override;

    //! Perturb pressure of data container
    /**
     * The method is used for the evaluation of numerical derivatives.
     * @param  deltaP  relative pressure perturbation
     */
    void perturbPressure(double deltaP);

    virtual void restore() override;

    virtual void invalidateCache() override {
        ReactionData::invalidateCache();
        pressure = NAN;
    }

    double pressure; //!< pressure
    double log10P; //!< base 10 logarithm of pressure

protected:
    double m_pressure_buf; //!< buffered pressure
};

}

#endif
