/**
 * @file ReactionData.h
 *
 * @warning This file is an experimental part of the %Cantera API and
 *    may be changed or removed without notice.
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
 *
 * @warning This class is an experimental part of the %Cantera API and
 *    may be changed or removed without notice.
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


//! Data container holding shared data specific to CustomFunc1Rate
/**
 * @warning This class is an experimental part of the %Cantera API and
 *    may be changed or removed without notice.
 */
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
