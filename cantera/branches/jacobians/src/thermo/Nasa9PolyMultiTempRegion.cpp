/**
 *  @file Nasa9PolyMultiTempRegion.cpp
 *  Definitions for a single-species standard state object derived
 *  from \link Cantera::SpeciesThermoInterpType
 *    SpeciesThermoInterpType\endlink  based
 *  on the NASA 9 coefficient temperature polynomial form
 *  applied to one temperature region
 *  (see \ref spthermo and class
 *   \link Cantera::Nasa9Poly1 Nasa9Poly1\endlink).
 *
 *  This parameterization has one NASA temperature region.
 */
// Copyright 2007  Sandia National Laboratories

#include "cantera/base/global.h"
#include "cantera/base/ctexceptions.h"
#include "Nasa9PolyMultiTempRegion.h"

using namespace std;

namespace Cantera
{

// The NASA 9 polynomial parameterization for a single species
// encompassing multiple temperature regions.
/*
 *  This parameterization expresses the heat capacity via a
 *  7 coefficient polynomial.
 *  Note that this is the form used in the
 *  2002 NASA equilibrium program. A reference to the form is
 *  provided below:
 *
 *  "NASA Glenn Coefficients for Calculating Thermodynamic
 *  Properties of Individual Species,"
 *  B. J. McBride, M. J. Zehe, S. Gordon
 *  NASA/TP-2002-211556, Sept. 2002
 *
 * Nine coefficients \f$(a_0,\dots,a_6)\f$ are used to represent
 * \f$ C_p^0(T)\f$, \f$ H^0(T)\f$, and \f$ S^0(T) \f$ as
 * polynomials in \f$ T \f$ :
 * \f[
 * \frac{c_p(T)}{R} = a_0 T^{-2} + a_1 T^{-1} + a_2 + a_3 T
 *                  + a_4 T^2 + a_5 T^3 + a_6 T^4
 * \f]
 *
 * \f[
 * \frac{H^0(T)}{RT} = - a_0 T^{-2} + a_1 \frac{\ln(T)}{T} + a_2
 * + a_3 T + a_4 T^2  + a_5 T^3 + a_6 T^4 + \frac{a_7}{T}
 * \f]
 *
 * \f[
 * \frac{s^0(T)}{R} = - \frac{a_0}{2} T^{-2} - a_1 T^{-1} + a_2 \ln(T)
 +    + a_3 T  \frac{a_4}{2} T^2 + \frac{a_5}{3} T^3
 *    + \frac{a_6}{4} T^4 + a_8
 * \f]
 *
 *  The standard state is assumed to be an ideal gas at the
 *  standard pressure of 1 bar, for gases.
 *  For condensed species, the standard state is the
 *  pure crystalline or liquid substance at the standard
 *  pressure of 1 atm.
 *
 * These NASA representations may have multiple temperature regions
 * through the use of this %Nasa9PolyMultiTempRegion object, which uses
 * multiple copies of the Nasa9Poly1 object to handle multiple temperature
 * regions.
 *
 * @ingroup spthermo
 */


//! Empty constructor
template<typename ValAndDerivType>
Nasa9PolyMultiTempRegion<ValAndDerivType>::Nasa9PolyMultiTempRegion() :
    m_lowT(0.0),
    m_highT(0.0),
    m_Pref(0.0),
    m_index(0),
    m_numTempRegions(0),
    m_currRegion(0)
{
}


// Constructor used in templated instantiations
/*
 * @param regionPts Vector of pointers to Nasa9Poly1 objects. These
 *                  objects all refer to the temperature regions for the
 *                  same species. The vector must be in increasing
 *                  temperature region format.  Together they
 *                  represent the reference temperature parameterization
 *                  for a single species.
 */
template<typename ValAndDerivType>
Nasa9PolyMultiTempRegion<ValAndDerivType>::
Nasa9PolyMultiTempRegion(std::vector<Cantera::Nasa9Poly1<ValAndDerivType> * > &regionPts) :
    m_lowT(0.0),
    m_highT(0.0),
    m_Pref(0.0),
    m_index(0),
    m_numTempRegions(0),
    m_currRegion(0)
{
    m_numTempRegions = regionPts.size();
    // Do a shallow copy of the pointers. From now on, we will
    // own these pointers and be responsible for deleting them.
    m_regionPts = regionPts;
    m_lowerTempBounds.resize(m_numTempRegions);
    m_lowT = m_regionPts[0]->minTemp();
    m_highT = m_regionPts[m_numTempRegions-1]->maxTemp();
    m_Pref = m_regionPts[0]->refPressure();
    m_index = m_regionPts[0]->speciesIndex();
    for (size_t i = 0; i < m_numTempRegions; i++) {
        m_lowerTempBounds[i] = m_regionPts[i]->minTemp();
        if (m_regionPts[i]->speciesIndex() != m_index) {
            throw CanteraError("Nasa9PolyMultiTempRegion::Nasa9PolyMultiTempRegion",
                               "m_index inconsistency");
        }
        if (fabs(m_regionPts[i]->refPressure() - m_Pref) > 0.0001) {
            throw CanteraError("Nasa9PolyMultiTempRegion::Nasa9PolyMultiTempRegion",
                               "refPressure inconsistency");
        }
        if (i > 0) {
            if (m_lowerTempBounds[i-1] >= m_lowerTempBounds[i]) {
                throw CanteraError("Nasa9PolyMultiTempRegion::Nasa9PolyMultiTempRegion",
                                   "minTemp bounds inconsistency");
            }
            if (fabs(m_regionPts[i-1]->maxTemp() - m_lowerTempBounds[i]) > 0.0001) {
                throw CanteraError("Nasa9PolyMultiTempRegion::Nasa9PolyMultiTempRegion",
                                   "Temp bounds inconsistency");
            }
        }
    }
}

// copy constructor
/*
 * @param b object to be copied
 */
template<typename ValAndDerivType>
Nasa9PolyMultiTempRegion<ValAndDerivType>::
Nasa9PolyMultiTempRegion(const Nasa9PolyMultiTempRegion<ValAndDerivType>& b) :
    m_lowT(b.m_lowT),
    m_highT(b.m_highT),
    m_Pref(b.m_Pref),
    m_index(b.m_index),
    m_numTempRegions(b.m_numTempRegions),
    m_lowerTempBounds(b.m_lowerTempBounds),
    m_currRegion(b.m_currRegion)
{
    m_regionPts.resize(m_numTempRegions);
    for (size_t i = 0; i < m_numTempRegions; i++) {
        Nasa9Poly1* dptr = b.m_regionPts[i];
        m_regionPts[i] = new Nasa9Poly1(*dptr);
    }
}

// assignment operator
/*
 * @param b object to be copied
 */
template<typename ValAndDerivType>
Nasa9PolyMultiTempRegion<ValAndDerivType>&
Nasa9PolyMultiTempRegion<ValAndDerivType>::operator=(const Nasa9PolyMultiTempRegion<ValAndDerivType>& b)
{
    if (&b != this) {
        for (size_t i = 0; i < m_numTempRegions; i++) {
            delete m_regionPts[i];
            m_regionPts[i] = 0;
        }
        m_lowT   = b.m_lowT;
        m_highT  = b.m_highT;
        m_Pref   = b.m_Pref;
        m_index  = b.m_index;
        m_numTempRegions = b.m_numTempRegions;
        m_lowerTempBounds = b.m_lowerTempBounds;
        m_currRegion = b.m_currRegion;
        m_regionPts.resize(m_numTempRegions);
        for (size_t i = 0; i < m_numTempRegions; i++) {
            m_regionPts[i] = new Nasa9Poly1(*(b.m_regionPts[i]));
        }
    }
    return *this;
}

// Destructor
template<typename ValAndDerivType>
Nasa9PolyMultiTempRegion<ValAndDerivType>::~Nasa9PolyMultiTempRegion()
{
    for (size_t i = 0; i < m_numTempRegions; i++) {
        delete m_regionPts[i];
        m_regionPts[i] = 0;
    }
}

// duplicator
template<typename ValAndDerivType>
SpeciesThermoInterpType<ValAndDerivType> *
Nasa9PolyMultiTempRegion<ValAndDerivType>::duplMyselfAsSpeciesThermoInterpType() const
{
    return new Nasa9PolyMultiTempRegion<ValAndDerivType>(*this);
}

// Returns the minimum temperature that the thermo
// parameterization is valid
template<typename ValAndDerivType>
doublereal Nasa9PolyMultiTempRegion<ValAndDerivType>::minTemp() const
{
    return m_lowT;
}

// Returns the maximum temperature that the thermo
// parameterization is valid
template<typename ValAndDerivType>
doublereal Nasa9PolyMultiTempRegion<ValAndDerivType>::maxTemp() const
{
    return m_highT;
}

// Returns the reference pressure (Pa)
template<typename ValAndDerivType>
doublereal Nasa9PolyMultiTempRegion<ValAndDerivType>::refPressure() const
{
    return m_Pref;
}

// Returns an integer representing the type of parameterization
template<typename ValAndDerivType>
int Nasa9PolyMultiTempRegion<ValAndDerivType>::reportType() const
{
    return NASA9MULTITEMP;
}


// Returns an integer representing the species index
template<typename ValAndDerivType>
size_t Nasa9PolyMultiTempRegion<ValAndDerivType>::speciesIndex() const
{
    return m_index;
}


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
 *  tt[2] = t*t*t;
 *  tt[3] = t*t*t*t;
 *  tt[4] = 1.0/t;
 *  tt[5] = 1.0/(t*t);
 *  tt[6] = std::log(t);
 *
 * @param tt      vector of temperature polynomials
 * @param cp_R    Vector of Dimensionless heat capacities.
 *                (length m_kk).
 * @param h_RT    Vector of Dimensionless enthalpies.
 *                (length m_kk).
 * @param s_R     Vector of Dimensionless entropies.
 *                (length m_kk).
 */
template<typename ValAndDerivType>
void Nasa9PolyMultiTempRegion<ValAndDerivType>::updateProperties(const doublereal* tt,
        ValAndDerivType* cp_R,
        ValAndDerivType* h_RT,
        ValAndDerivType* s_R) const
{
    // Let's put some additional debugging here.
    // This is an external routine
#ifdef DEBUG_HKM
    double temp = tt[0];
    if (temp < m_regionPts[m_currRegion]->minTemp()) {
        if (m_currRegion != 0) {
            throw CanteraError("Nasa9PolyMultiTempRegion::updateProperties",
                               "region problem");
        }
    }
    if (temp > m_regionPts[m_currRegion]->maxTemp()) {
        if (m_currRegion != m_numTempRegions - 1) {
            throw CanteraError("Nasa9PolyMultiTempRegion::updateProperties",
                               "region problem");
        }
    }
#endif
    (m_regionPts[m_currRegion])->updateProperties(tt, cp_R, h_RT, s_R);
}


// Compute the reference-state property of one species
/*
 * Given temperature T in K, this method updates the values of
 * the non-dimensional heat capacity at constant pressure,
 * enthalpy, and entropy, at the reference pressure, Pref
 * of one of the species. The species index is used
 * to reference into the cp_R, h_RT, and s_R arrays.
 *
 * Temperature Polynomial:
 *  tt[0] = t;
 *  tt[1] = t*t;
 *  tt[2] = t*t*t;
 *  tt[3] = t*t*t*t;
 *  tt[4] = 1.0/t;
 *  tt[5] = 1.0/(t*t);
 *  tt[6] = std::log(t);
 *
 * @param temp    Temperature (Kelvin)
 * @param cp_R    Vector of Dimensionless heat capacities.
 *                (length m_kk).
 * @param h_RT    Vector of Dimensionless enthalpies.
 *                (length m_kk).
 * @param s_R     Vector of Dimensionless entropies.
 *                (length m_kk).
 */
template<typename ValAndDerivType>
void Nasa9PolyMultiTempRegion<ValAndDerivType>::updatePropertiesTemp(const doublereal temp,
        ValAndDerivType* cp_R, ValAndDerivType* h_RT,
        ValAndDerivType* s_R) const
{
    double tPoly[7];
    tPoly[0]  = temp;
    tPoly[1]  = temp * temp;
    tPoly[2]  = tPoly[1] * temp;
    tPoly[3]  = tPoly[2] * temp;
    tPoly[4]  = 1.0 / temp;
    tPoly[5]  = tPoly[4] / temp;
    tPoly[6]  = std::log(temp);
    // Now find the region
    m_currRegion = 0;
    for (size_t i = 1; i < m_numTempRegions; i++) {
        if (temp < m_lowerTempBounds[i]) {
            break;
        }
        m_currRegion++;
    }

    updateProperties(tPoly, cp_R, h_RT, s_R);
}

//This utility function reports back the type of
// parameterization and all of the parameters for the
// species, index.
/*
 * All parameters are output variables
 *
 * @param n         Species index
 * @param type      Integer type of the standard type
 * @param tlow      output - Minimum temperature
 * @param thigh     output - Maximum temperature
 * @param pref      output - reference pressure (Pa).
 * @param coeffs    Vector of coefficients used to set the
 *                  parameters for the standard state.
 */
template<typename ValAndDerivType>
void Nasa9PolyMultiTempRegion<ValAndDerivType>::reportParameters(size_t& n, int& type,
        doublereal& tlow, doublereal& thigh,
        doublereal& pref,
        doublereal* const coeffs) const
{
    n = m_index;
    type = NASA9MULTITEMP;
    tlow = m_lowT;
    thigh = m_highT;
    pref = m_Pref;
    double ctmp[12];
    coeffs[0] = double(m_numTempRegions);
    int index = 1;
    size_t n_tmp = 0;
    int type_tmp = 0;
    double pref_tmp = 0.0;
    for (size_t iReg = 0; iReg < m_numTempRegions; iReg++) {
        m_regionPts[iReg]->reportParameters(n_tmp, type_tmp,
                                            coeffs[index], coeffs[index+1],
                                            pref_tmp, ctmp);
        for (int i = 0; i < 9; i++) {
            coeffs[index+2+i] = ctmp[3+i];
        }
        index += 11;
    }

}

// Modify parameters for the standard state
/*
 * @param coeffs   Vector of coefficients used to set the
 *                 parameters for the standard state.
 */
template<typename ValAndDerivType>
void Nasa9PolyMultiTempRegion<ValAndDerivType>::modifyParameters(doublereal* coeffs)
{
    int index = 3;
    for (size_t iReg = 0; iReg < m_numTempRegions; iReg++) {
        m_regionPts[iReg]->modifyParameters(coeffs + index);
        index += 11;
    }
}

}
