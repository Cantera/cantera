/**
 *  @file NasaPoly2.cpp
 *  \link Cantera::NasaPoly2 NasaPoly2\endlink
 */

// Copyright 2007  Sandia National Laboratories
#include "NasaPoly2.h"

#include <vector>
#include <map>

namespace Cantera {

// Update the properties for this species, given a temperature polynomial
/*
 * This method is called with a pointer to an array containing the functions of
 * temperature needed by this  parameterization, and three pointers to arrays where the
 * computed property values should be written. This method updates only one value in
 * each array.
 *
 * Temperature Polynomial:
 *  tt[0] = t;
 *  tt[1] = t*t;
 *  tt[2] = m_t[1]*t;
 *  tt[3] = m_t[2]*t;
 *  tt[4] = 1.0/t;
 *  tt[5] = std::log(t);
 *
 * @param tt  vector of temperature polynomials
 * @param cp_R    Vector of Dimensionless heat capacities.
 *                (length m_kk).
 * @param h_RT    Vector of Dimensionless enthalpies.
 *                (length m_kk).
 * @param s_R     Vector of Dimensionless entropies.
 *                (length m_kk).
 */
template<typename ValAndDerivType>
void NasaPoly2<ValAndDerivType>::updateProperties(const ValAndDerivType* tt, ValAndDerivType* cp_R, ValAndDerivType* h_RT, ValAndDerivType* s_R) const
{

    double T = tt[0];
    if (T <= m_midT) {
        mnp_low.updateProperties(tt, cp_R, h_RT, s_R);
    } else {
        mnp_high.updateProperties(tt, cp_R, h_RT, s_R);
    }
}

// Overriding the particular template

template<> void NasaPoly2<doubleFAD>::updateProperties(const doubleFAD* tt, doubleFAD* cp_R, doubleFAD* h_RT, doubleFAD* s_R) const
{

    double T = tt[0].val();
    if (T <= m_midT) {
        mnp_low.updateProperties(tt, cp_R, h_RT, s_R);
    } else {
        mnp_high.updateProperties(tt, cp_R, h_RT, s_R);
    }
}

/*
 *  HKM For some reason this tidbit was needed to complete the linking.
 */
template void NasaPoly2<doublereal>::updateProperties(const doublereal* tt, doublereal* cp_R, doublereal* h_RT, doublereal* s_R) const;

}
