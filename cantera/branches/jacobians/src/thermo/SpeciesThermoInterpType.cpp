/**
 *  @file SpeciesThermoInterpType.cpp
 *  Definitions for a
 */
// Copyright 2007  Sandia National Laboratories
#include "cantera/thermo/SpeciesThermoInterpType.h"
#include "cantera/thermo/VPSSMgr.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/ct_defs.h"

namespace Cantera {

template<typename ValAndDerivType>
SpeciesThermoInterpType<ValAndDerivType>::SpeciesThermoInterpType()
{
}

template<typename ValAndDerivType>
SpeciesThermoInterpType<ValAndDerivType>::~SpeciesThermoInterpType()
{
}

template<typename ValAndDerivType>
void SpeciesThermoInterpType<ValAndDerivType>::updateProperties(const ValAndDerivType* tempPoly, ValAndDerivType* cp_R,
                                                                ValAndDerivType* h_RT, ValAndDerivType* s_R) const
{
    double T = tempPoly[0];
    updatePropertiesTemp(T, cp_R, h_RT, s_R);
}

// Specialized instantiation of the updateProperties() member function for the DFad template case.
/*
 * This is necessary because we need to get the plain value of temperature for updatePropertiesTemp() function
 */
template<> void SpeciesThermoInterpType<doubleFAD>::updateProperties(const doubleFAD* tempPoly, doubleFAD* cp_R,
                                                                doubleFAD* h_RT, doubleFAD* s_R) const
{
    double T = tempPoly[0].val();
    updatePropertiesTemp(T, cp_R, h_RT, s_R);
}

#ifdef H298MODIFY_CAPABILITY

template<typename ValAndDerivType>
doublereal SpeciesThermoInterpType<ValAndDerivType>::reportHf298(doublereal* const h298) const
{
    throw CanteraError("SpeciesThermoInterpType::reportHf298", "Not implemented");
}

template<typename ValAndDerivType>
void SpeciesThermoInterpType<ValAndDerivType>::modifyOneHf298(const size_t k, const doublereal Hf298New)
{
    throw CanteraError("SpeciesThermoInterpType::modifyOneHf298", "Not implemented");
}

#endif

/***************************************************************************************************/

//! Constructor
template<typename ValAndDerivType>
STITbyPDSS<ValAndDerivType>::STITbyPDSS() :
        m_speciesIndex(npos)
{
}

template<typename ValAndDerivType>
STITbyPDSS<ValAndDerivType>::STITbyPDSS(size_t k, VPSSMgr* vpssmgr_ptr, PDSS* PDSS_ptr) :
        m_vpssmgr_ptr(vpssmgr_ptr),
        m_PDSS_ptr(PDSS_ptr),
        m_speciesIndex(k)
{
}

template<typename ValAndDerivType>
STITbyPDSS<ValAndDerivType>::STITbyPDSS(const STITbyPDSS& b) :
        m_vpssmgr_ptr(b.m_vpssmgr_ptr),
        m_PDSS_ptr(b.m_PDSS_ptr),
        m_speciesIndex(b.m_speciesIndex)
{
}

//! Destructor
template<typename ValAndDerivType>
STITbyPDSS<ValAndDerivType>::~STITbyPDSS()
{
}

//! duplicator
template<typename ValAndDerivType>
SpeciesThermoInterpType<ValAndDerivType> *
STITbyPDSS<ValAndDerivType>::duplMyselfAsSpeciesThermoInterpType() const
{
    return new STITbyPDSS<ValAndDerivType>(*this);
}

template<typename ValAndDerivType>
void STITbyPDSS<ValAndDerivType>::initAllPtrs(size_t speciesIndex, VPSSMgr* vpssmgr_ptr, PDSS* PDSS_ptr)
{
    AssertThrow(speciesIndex == m_speciesIndex, "STITbyPDSS::initAllPtrs internal confusion");
    m_vpssmgr_ptr = vpssmgr_ptr;
    m_PDSS_ptr = PDSS_ptr;
}

//! Returns the minimum temperature that the thermo
//! parameterization is valid
template<typename ValAndDerivType>
doublereal STITbyPDSS<ValAndDerivType>::minTemp() const
{
    return m_PDSS_ptr->minTemp();
}

//! Returns the maximum temperature that the thermo
//! parameterization is valid
template<typename ValAndDerivType>
doublereal STITbyPDSS<ValAndDerivType>::maxTemp() const
{
    return m_PDSS_ptr->maxTemp();
}

//! Returns the reference pressure (Pa)
template<typename ValAndDerivType>
doublereal STITbyPDSS<ValAndDerivType>::refPressure() const
{
    return m_PDSS_ptr->refPressure();
}

//! Returns an integer representing the type of parameterization
template<typename ValAndDerivType>
int STITbyPDSS<ValAndDerivType>::reportType() const
{
    return PDSS_TYPE;
}

//! Returns an integer representing the species index
template<typename ValAndDerivType>
size_t STITbyPDSS<ValAndDerivType>::speciesIndex() const
{
    return m_speciesIndex;
}

//! Update the properties for this species, given a temperature
//! polynomial
/*!
 * This method is called with a pointer to an array containing the functions of
 * temperature needed by this  parameterization, and three pointers to arrays where the
 * computed property values should be written. This method updates only one value in
 * each array.
 *
 * The form and length of the Temperature Polynomial may vary depending on the
 * parameterization.
 *
 * @param tempPoly  vector of temperature polynomials
 * @param cp_R    Vector of Dimensionless heat capacities.
 *                (length m_kk).
 * @param h_RT    Vector of Dimensionless enthalpies.
 *                (length m_kk).
 * @param s_R     Vector of Dimensionless entropies.
 *                (length m_kk).
 */
template<typename ValAndDerivType>
void STITbyPDSS<ValAndDerivType>::updateProperties(const ValAndDerivType* tempPoly, ValAndDerivType* cp_R, ValAndDerivType* h_RT,
                                                   ValAndDerivType* s_R) const
{
    doublereal T = tempPoly[0];
    updatePropertiesTemp(T, cp_R, h_RT, s_R);
}

//! Compute the reference-state property of one species
/*!
 * Given temperature T in K, this method updates the values of
 * the non-dimensional heat capacity at constant pressure,
 * enthalpy, and entropy, at the reference pressure, Pref
 * of one of the species. The species index is used
 * to reference into the cp_R, h_RT, and s_R arrays.
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
void STITbyPDSS<ValAndDerivType>::updatePropertiesTemp(const doublereal temp, ValAndDerivType* cp_R, ValAndDerivType* h_RT,
                                                       ValAndDerivType* s_R) const
{
    //m_vpssmgr_ptr->setState_T(temp);
    m_PDSS_ptr->setTemperature(temp);
    AssertThrowMsg(m_speciesIndex != npos, "STITbyPDSS::updatePropertiesTemp", "object was probably not installed correctly");
    h_RT[m_speciesIndex] = m_PDSS_ptr->enthalpy_RT_ref();
    cp_R[m_speciesIndex] = m_PDSS_ptr->cp_R_ref();
    s_R[m_speciesIndex] = m_PDSS_ptr->entropy_R_ref();
}

//!This utility function reports back the type of
//! parameterization and all of the parameters for the
//! species, index.
/*!
 * All parameters are output variables
 *
 * @param index     Species index
 * @param type      Integer type of the standard type
 * @param minTemp   output - Minimum temperature
 * @param maxTemp   output - Maximum temperature
 * @param refPressure output - reference pressure (Pa).
 * @param coeffs    Vector of coefficients used to set the
 *                  parameters for the standard state.
 */
template<typename ValAndDerivType>
void STITbyPDSS<ValAndDerivType>::reportParameters(size_t& index, int& type, doublereal& minTemp, doublereal& maxTemp,
                                                   doublereal& refPressure, doublereal* const coeffs) const
{
    index = m_speciesIndex;
    type = PDSS_TYPE;
    minTemp = m_vpssmgr_ptr->minTemp(m_speciesIndex);
    maxTemp = m_vpssmgr_ptr->maxTemp(m_speciesIndex);
    refPressure = m_PDSS_ptr->refPressure();
}

//! Modify parameters for the standard state
/*!
 * @param coeffs   Vector of coefficients used to set the
 *                 parameters for the standard state.
 */
template<typename ValAndDerivType>
void STITbyPDSS<ValAndDerivType>::modifyParameters(doublereal* coeffs)
{
}

// ExplicitInstantiation Section
template class SpeciesThermoInterpType<doublereal>;
#ifdef INDEPENDENT_VARIABLE_DERIVATIVES
#ifdef HAS_SACADO
template class SpeciesThermoInterpType<doubleFAD>;
#endif
#endif

}
