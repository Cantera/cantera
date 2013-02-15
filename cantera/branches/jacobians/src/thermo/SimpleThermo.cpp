/**
 *  @file ShomatePoly.h
 *  Header for a single-species standard state object derived
 *  from \link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType\endlink  based
 *  on the Shomate temperature polynomial form applied to one temperature region
 *  (see \ref spthermo and class \link Cantera::ShomatePoly ShomatePoly\endlink and
 *   \link Cantera::ShomatePoly2 ShomatePoly2\endlink).
 *    Shomate polynomial expressions.
 */
// Copyright 2001  California Institute of Technology
#include "cantera/thermo/SimpleThermo.h"
#include "cantera/Cantera.h"

namespace Cantera {

// Assignment operator
/*
 * @param right Object to be copied
 */
template<typename ValAndDerivType>
SimpleThermo<ValAndDerivType>& SimpleThermo<ValAndDerivType>::operator=(const SimpleThermo<ValAndDerivType>& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    m_loc = right.m_loc;
    m_index = right.m_index;
    m_tlow_max = right.m_tlow_max;
    m_thigh_min = right.m_thigh_min;
    m_tlow = right.m_tlow;
    m_thigh = right.m_thigh;
    m_t0 = right.m_t0;
    m_logt0 = right.m_logt0;
    m_h0_R = right.m_h0_R;
    m_s0_R = right.m_s0_R;
    m_cp0_R = right.m_cp0_R;
    m_p0 = right.m_p0;
    m_nspData = right.m_nspData;

    return *this;
}

// Assignment operator
/*
 * @param right Object to be copied
 */
template<typename ValAndDerivType>
template<typename ValAndDerivType2>
SimpleThermo<ValAndDerivType>& SimpleThermo<ValAndDerivType>::operator=(const SimpleThermo<ValAndDerivType2>& right)
{
    /*
     * Check for self assignment.
     */
    if (this == (SimpleThermo<ValAndDerivType> *) &right) {
        return *this;
    }

    m_loc = right.m_loc;
    m_index = right.m_index;
    m_tlow_max = right.m_tlow_max;
    m_thigh_min = right.m_thigh_min;
    m_tlow = right.m_tlow;
    m_thigh = right.m_thigh;
    m_t0 = right.m_t0;
    m_logt0 = right.m_logt0;
    m_h0_R = right.m_h0_R;
    m_s0_R = right.m_s0_R;
    m_cp0_R = right.m_cp0_R;
    m_p0 = right.m_p0;
    m_nspData = right.m_nspData;

    return *this;
}

// Assignment operator
/*
 * @param right Object to be copied
 */
template<>
template<>
SimpleThermo<doublereal>& SimpleThermo<doublereal>::operator=(const SimpleThermo<doubleFAD>& right)
{
    /*
     * Check for self assignment.
     */
    if (this == (SimpleThermo<doublereal> *) &right) {
        return *this;
    }

    m_loc = right.m_loc;
    m_index = right.m_index;
    m_tlow_max = right.m_tlow_max;
    m_thigh_min = right.m_thigh_min;
    m_tlow = right.m_tlow;
    m_thigh = right.m_thigh;
    m_t0 = right.m_t0;
    m_logt0 = right.m_logt0;
    m_h0_R = right.m_h0_R;
    m_s0_R = right.m_s0_R;
    m_cp0_R = right.m_cp0_R;
    m_p0 = right.m_p0;
    m_nspData = right.m_nspData;

    return *this;
}


template class SimpleThermo<doublereal>;
#ifdef INDEPENDENT_VARIABLE_DERIVATIVES
#ifdef HAS_SACADO
template class SimpleThermo<doubleFAD>;
#endif
#endif


}
