/**
 *  @file NasaPoly.h
 *  Header for a single-species standard state object derived
 *  from \link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType\endlink  based
 *  on the Nasa temperature polynomial form applied to one temperature region
 *  (see \ref spthermo and class \link Cantera::NasaPoly NasaPoly\endlink and
 *   \link Cantera::NasaPoly2 NasaPoly2\endlink).
 *    Nasa polynomial expressions.
 */
// Copyright 2001  California Institute of Technology
#include "NasaThermo.h"

namespace Cantera {

//================================================================================================================================

template<typename ValAndDerivType>
NasaThermo<ValAndDerivType>& NasaThermo<ValAndDerivType>::operator=(const NasaThermo<ValAndDerivType>& right)
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
    m_name = right.m_name;

    return *this;
}
//================================================================================================================================
template<typename ValAndDerivType>
template<typename ValAndDerivType2>
NasaThermo<ValAndDerivType>& NasaThermo<ValAndDerivType>::operator=(const NasaThermo<ValAndDerivType2>& right)
{
    if ((NasaThermo<ValAndDerivType> *) &right == this) {
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
    m_name = right.m_name;

    return *this;
}

template<>
template<>
NasaThermo<doublereal>&
NasaThermo<doublereal>::operator=(const NasaThermo<doubleFAD>& right)
{
    if ((NasaThermo<doublereal> *) &right == this) {
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
    m_name = right.m_name;

    return *this;
}


template class NasaThermo<doublereal>;
#ifdef INDEPENDENT_VARIABLE_DERIVATIVES
#ifdef HAS_SACADO
template class NasaThermo<doubleFAD>;
template NasaThermo<doublereal>& NasaThermo<doublereal>::operator=(const NasaThermo<doubleFAD>& right);

#endif
#endif


}
