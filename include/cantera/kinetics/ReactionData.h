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


//! Data container holding shared data specific to ArrheniusRate
/**
 * The data container `ArrheniusData` holds precalculated data common to
 * all `ArrheniusRate` objects.
 */
struct ArrheniusData
{
    ArrheniusData() : temperature(1.), logT(0.), recipT(1.) {}

    //! Update data container based on temperature *T*
    void update(double T) {
        temperature = T;
        logT = std::log(T);
        recipT = 1./T;
    }

    //! Update data container based on temperature *T* and pressure *P*
    void update(double T, double P) { update(T); }

    //! Update data container based on *bulk* phase state
    void update(const ThermoPhase& bulk, const Kinetics& kin);

    //! Update number of species and reactions; unused
    void resize(size_t n_species, size_t n_reactions) {}

    double temperature; //!< temperature
    double logT; //!< logarithm of temperature
    double recipT;  //!< inverse of temperature
};


//! Data container holding shared data specific to BlowersMaselRate
/**
 * The data container `BlowersMaselData` holds precalculated data common to
 * all `BlowersMaselRate` objects.
 */
struct BlowersMaselData
{
    BlowersMaselData()
        : temperature(1.), logT(0.), recipT(1.), finalized(false) {}

    //! Update data container based on temperature *T*
    void update(double T);

    //! Update data container based on temperature *T* and pressure *P*
    void update(double T, double P) {
        update(T);
    }

    //! Update data container based on *bulk* phase state and *kin* kinetics
    void update(const ThermoPhase& bulk, const Kinetics& kin);

    //! Finalize setup
    void resize(size_t n_species, size_t n_reactions) {
        m_grt.resize(n_species, 0.);
        dH.resize(n_reactions, 0.);
        finalized = true;
    }

    double temperature; //!< temperature
    double logT; //!< logarithm of temperature
    double recipT; //!< inverse of temperature

    bool finalized; //!< boolean indicating whether vectors are accessible
    vector_fp dH; //!< enthalpy change for each reaction

protected:
    vector_fp m_grt; //!< work vector holding partial molar enthalpies
};


//! Data container holding shared data specific to PlogRate
/**
 * The data container `PlogData` holds precalculated data common to
 * all `PlogRate` objects.
 */
struct PlogData
{
    PlogData() : temperature(1.), logT(0.), recipT(1.), logP(0.) {}

    //! Update data container based on temperature *T* (raises exception)
    void update(double T);

    //! Update data container based on temperature *T* and *P*
    void update(double T, double P) {
        temperature = T;
        logT = std::log(T);
        recipT = 1./T;
        logP = std::log(P);
    }

    //! Update data container based on *bulk* phase state
    void update(const ThermoPhase& bulk, const Kinetics& kin);

    //! Update number of species and reactions; unused
    void resize(size_t n_species, size_t n_reactions) {}

    double temperature; //!< temperature
    double logT; //!< logarithm of temperature
    double recipT; //!< inverse of temperature
    double logP; //!< logarithm of pressure
};


//! Data container holding shared data specific to ChebyshevRate
/**
 * The data container `ChebyshevData` holds precalculated data common to
 * all `ChebyshevRate3` objects.
 */
struct ChebyshevData
{
    ChebyshevData() : temperature(1.), recipT(1.), log10P(0.) {}

    //! Update data container based on temperature *T* (raises exception)
    void update(double T);

    //! Update data container based on temperature *T* and *P*
    void update(double T, double P)
    {
        temperature = T;
        recipT = 1./T;
        log10P = std::log10(P);
    }

    //! Update data container based on *bulk* phase state
    void update(const ThermoPhase& bulk, const Kinetics& kin);

    //! Update number of species and reactions; unused
    void resize(size_t n_species, size_t n_reactions) {}

    double temperature; //!< temperature
    double recipT; //!< inverse of temperature
    double log10P; //!< base 10 logarithm of pressure
};


//! Data container holding shared data specific to CustomFunc1Rate
struct CustomFunc1Data
{
    CustomFunc1Data() : temperature(1.) {}

    //! Update data container based on temperature *T*
    void update(double T) { temperature = T; }

    //! Update data container based on temperature *T* and pressure *P*
    void update(double T, double P) { update(T); }

    //! Update data container based on *bulk* phase state
    void update(const ThermoPhase& bulk, const Kinetics& kin);

    //! Update number of species and reactions; unused
    void resize(size_t n_species, size_t n_reactions) {}

    double temperature; //!< temperature
};

}

#endif
