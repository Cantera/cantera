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
    ArrheniusData() : m_temperature(1.), m_logT(0.), m_recipT(1.) {}

    //! Update data container based on temperature *T*
    void update(double T)
    {
        m_temperature = T;
        m_logT = std::log(T);
        m_recipT = 1./T;
    }

    //! Update data container based on temperature *T* and pressure *P*
    void update(double T, double P) { update(T); }

    //! Update data container based on *bulk* phase state
    void update(const ThermoPhase& bulk, const Kinetics& kin);

    //! Update number of species and reactions; unused
    void resize(size_t n_species, size_t n_reactions) {}

    double m_temperature; //!< temperature
    double m_logT; //!< logarithm of temperature
    double m_recipT;  //!< inverse of temperature
};


//! Data container holding shared data specific to PlogRate
/**
 * The data container `PlogData` holds precalculated data common to
 * all `PlogRate` objects.
 */
struct PlogData
{
    PlogData() : m_temperature(1.), m_logT(0.), m_recipT(1.), m_logP(0.) {}

    //! Update data container based on temperature *T* (raises exception)
    void update(double T);

    //! Update data container based on temperature *T* and *P*
    void update(double T, double P)
    {
        m_temperature = T;
        m_logT = std::log(T);
        m_recipT = 1./T;
        m_logP = std::log(P);
   }

    //! Update data container based on *bulk* phase state
    void update(const ThermoPhase& bulk, const Kinetics& kin);

    //! Update number of species and reactions; unused
    void resize(size_t n_species, size_t n_reactions) {}

    double m_temperature; //!< temperature
    double m_logT; //!< logarithm of temperature
    double m_recipT; //!< inverse of temperature
    double m_logP; //!< logarithm of pressure
};


//! Data container holding shared data specific to ChebyshevRate
/**
 * The data container `ChebyshevData` holds precalculated data common to
 * all `ChebyshevRate3` objects.
 */
struct ChebyshevData
{
    ChebyshevData() : m_temperature(1.), m_recipT(1.), m_log10P(0.) {}

    //! Update data container based on temperature *T* (raises exception)
    void update(double T);

    //! Update data container based on temperature *T* and *P*
    void update(double T, double P)
    {
        m_temperature = T;
        m_recipT = 1./T;
        m_log10P = std::log10(P);
    }

    //! Update data container based on *bulk* phase state
    void update(const ThermoPhase& bulk, const Kinetics& kin);

    //! Update number of species and reactions; unused
    void resize(size_t n_species, size_t n_reactions) {}

    double m_temperature; //!< temperature
    double m_recipT; //!< inverse of temperature
    double m_log10P; //!< base 10 logarithm of pressure
};


//! Data container holding shared data specific to CustomFunc1Rate
struct CustomFunc1Data
{
    CustomFunc1Data() : m_temperature(1.) {}

    //! Update data container based on temperature *T*
    void update(double T) { m_temperature = T; }

    //! Update data container based on temperature *T* and pressure *P*
    void update(double T, double P) { update(T); }

    //! Update data container based on *bulk* phase state
    void update(const ThermoPhase& bulk, const Kinetics& kin);

    //! Update number of species and reactions; unused
    void resize(size_t n_species, size_t n_reactions) {}

    double m_temperature; //!< temperature
};

}

#endif
