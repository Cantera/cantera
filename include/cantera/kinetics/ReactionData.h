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


//! Data container holding shared data specific to ArrheniusRate
/**
 * The data container `ArrheniusData` holds precalculated data common to
 * all `ArrheniusRate` objects.
 */
struct ArrheniusData
{
    ArrheniusData() : m_temperature(1.), m_logT(0.), m_recipT(1.) {}

    //! Constructor based on temperature *T*
    ArrheniusData(double T) { update(T); }

    //! Constructor based on temperature *T* and pressure *P*
    ArrheniusData(double T, double P) { update(T); };

    //! Constructor accessing *bulk* phase definitions
    ArrheniusData(const ThermoPhase& bulk) { update(bulk); }

    void update(double T) {
        m_temperature = T;
        m_logT = std::log(T);
        m_recipT = 1./T;
    }

    void update(const ThermoPhase& bulk);

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

    //! Constructor based on temperature *T* and pressure *P*
    PlogData(double T);

    //! Constructor based on temperature *T* and pressure *P*
    PlogData(double T, double P) { update(T, P); };

    //! Constructor accessing *bulk* phase definitions
    PlogData(const ThermoPhase& bulk) { update(bulk); }

    void update(double T);

    void update(double T, double P) {
        m_temperature = T;
        m_logT = std::log(T);
        m_recipT = 1./T;
        m_logP = std::log(P);
   }

    void update(const ThermoPhase& bulk);

    //! Pointer to logP (required by Plog::update_C)
    const double* logP() const { return &m_logP; }

    double m_temperature; //!< temperature
    double m_logT; //!< logarithm of temperature
    double m_recipT;  //!< inverse of temperature
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

    //! Constructor based on temperature *T* and pressure *P*
    ChebyshevData(double T);

    //! Constructor based on temperature *T* and pressure *P*
    ChebyshevData(double T, double P) { update(T, P); };

    //! Constructor accessing *bulk* phase definitions
    ChebyshevData(const ThermoPhase& bulk) { update(bulk); }

    void update(double T);

    void update(double T, double P) {
        m_temperature = T;
        m_recipT = 1./T;
        m_log10P = std::log10(P);
   }

    void update(const ThermoPhase& bulk);

    //! Pointer to logP (required by Chebyshev::update_C)
    const double* log10P() const { return &m_log10P; }

    double m_temperature; //!< temperature
    double m_recipT;  //!< inverse of temperature
    double m_log10P; //!< base 10 logarithm of pressure
};


//! Data container holding shared data specific to CustomFunc1Rate
struct CustomFunc1Data
{
    CustomFunc1Data() : m_temperature(1.) {}

    //! Constructor based on temperature *T*
    CustomFunc1Data(double T) { update(T); }

    //! Constructor based on temperature *T* and pressure *P*
    CustomFunc1Data(double T, double P) { update(T); };

    //! Constructor accessing *bulk* phase definitions
    CustomFunc1Data(const ThermoPhase& bulk) { update(bulk); }

    void update(double T) { m_temperature = T; }

    void update(const ThermoPhase& bulk);

    double m_temperature; //!< temperature
};

}

#endif
