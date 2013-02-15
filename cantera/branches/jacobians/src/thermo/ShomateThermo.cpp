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
#include "ShomateThermo.h"

namespace Cantera {

//================================================================================================================================

template<typename ValAndDerivType>
ShomateThermo<ValAndDerivType>& ShomateThermo<ValAndDerivType>::operator=(const ShomateThermo<ValAndDerivType>& right)
{
    if (&right == this) {
        return *this;
    }

    m_high = right.m_high;
    m_low = right.m_low;
    m_index = right.m_index;
    m_tmid = right.m_tmid;
    m_tlow_max = right.m_tlow_max;
    m_thigh_min = right.m_thigh_min;
    m_tlow = right.m_tlow;
    m_thigh = right.m_thigh;
    m_p0 = right.m_p0;
    m_ngroups = right.m_ngroups;
    m_t = right.m_t;
    m_group_map = right.m_group_map;
    m_posInGroup_map = right.m_posInGroup_map;

    return *this;
}

//================================================================================================================================

template<typename ValAndDerivType>
template<typename ValAndDerivType2>
ShomateThermo<ValAndDerivType>& ShomateThermo<ValAndDerivType>::operator=(const ShomateThermo<ValAndDerivType2>& right)
{
    if ((ShomateThermo<ValAndDerivType> *) &right == this) {
        return *this;
    }

    m_high = right.m_high;
    m_low = right.m_low;
    m_index = right.m_index;
    m_tmid = right.m_tmid;
    m_tlow_max = right.m_tlow_max;
    m_thigh_min = right.m_thigh_min;
    m_tlow = right.m_tlow;
    m_thigh = right.m_thigh;
    m_p0 = right.m_p0;
    m_ngroups = right.m_ngroups;
    m_t = right.m_t;
    m_group_map = right.m_group_map;
    m_posInGroup_map = right.m_posInGroup_map;

    return *this;
}

template<>
template<>
ShomateThermo<doublereal>& ShomateThermo<doublereal>::operator=(const ShomateThermo<doubleFAD>& right)
{
    if ((ShomateThermo<doublereal> *) &right == this) {
        return *this;
    }

    assignVectorVectorClass(m_high, right.m_high);
    assignVectorVectorClass(m_low, right.m_low);

    m_index = right.m_index;
    m_tmid = right.m_tmid;
    m_tlow_max = right.m_tlow_max;
    m_thigh_min = right.m_thigh_min;
    m_tlow = right.m_tlow;
    m_thigh = right.m_thigh;
    m_p0 = right.m_p0;
    m_ngroups = right.m_ngroups;

    assignVectorFadToDouble(m_t, right.m_t);

    m_group_map = right.m_group_map;
    m_posInGroup_map = right.m_posInGroup_map;

    return *this;
}


template class ShomateThermo<doublereal>;
#ifdef INDEPENDENT_VARIABLE_DERIVATIVES
#ifdef HAS_SACADO
template class ShomateThermo<doubleFAD>;
#endif
#endif


}
