/**
 *  @file SpeciesThermoInterpType.cpp
 *  Definitions for a
 */
// Copyright 2007  Sandia National Laboratories

#include "cantera/thermo/SpeciesThermoInterpType.h"
#include "cantera/thermo/VPSSMgr.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

SpeciesThermoInterpType::SpeciesThermoInterpType()
{
}

SpeciesThermoInterpType::~SpeciesThermoInterpType()
{
}

void SpeciesThermoInterpType::updateProperties(const doublereal* tempPoly,
        doublereal* cp_R, doublereal* h_RT,
        doublereal* s_R) const
{
    double T = tempPoly[0];
    updatePropertiesTemp(T, cp_R, h_RT, s_R);
}

#ifdef H298MODIFY_CAPABILITY

doublereal SpeciesThermoInterpType::reportHf298(doublereal* const h298) const
{
    throw CanteraError("SpeciesThermoInterpType::reportHf298",
                       "Not implemented");
}

void SpeciesThermoInterpType::modifyOneHf298(const int k, const doublereal Hf298New)
{
    throw CanteraError("SpeciesThermoInterpType::modifyOneHf298",
                       "Not implemented");
}

#endif

/***************************************************************************************************/

//! Constructor
STITbyPDSS::STITbyPDSS() :
    m_speciesIndex(npos)
{
}

STITbyPDSS::STITbyPDSS(size_t k, VPSSMgr* vpssmgr_ptr, PDSS* PDSS_ptr) :
    m_vpssmgr_ptr(vpssmgr_ptr),
    m_PDSS_ptr(PDSS_ptr),
    m_speciesIndex(k)
{
}

STITbyPDSS::STITbyPDSS(const STITbyPDSS& b) :
    m_vpssmgr_ptr(b.m_vpssmgr_ptr),
    m_PDSS_ptr(b.m_PDSS_ptr),
    m_speciesIndex(b.m_speciesIndex)
{
}

//! Destructor
STITbyPDSS::~STITbyPDSS()
{
}

//! duplicator
SpeciesThermoInterpType*
STITbyPDSS::duplMyselfAsSpeciesThermoInterpType() const
{
    return new STITbyPDSS(*this);
}


void STITbyPDSS::initAllPtrs(size_t speciesIndex, VPSSMgr* vpssmgr_ptr, PDSS* PDSS_ptr)
{
    AssertThrow(speciesIndex == m_speciesIndex, "STITbyPDSS::initAllPtrs internal confusion");
    m_vpssmgr_ptr = vpssmgr_ptr;
    m_PDSS_ptr = PDSS_ptr;
}

//! Returns the minimum temperature that the thermo
//! parameterization is valid
doublereal  STITbyPDSS::minTemp() const
{
    return m_PDSS_ptr->minTemp();
}

//! Returns the maximum temperature that the thermo
//! parameterization is valid
doublereal  STITbyPDSS::maxTemp() const
{
    return m_PDSS_ptr->maxTemp();
}

//! Returns the reference pressure (Pa)
doublereal  STITbyPDSS::refPressure() const
{
    return m_PDSS_ptr->refPressure();
}

//! Returns an integer representing the type of parameterization
int  STITbyPDSS::reportType() const
{
    return PDSS_TYPE;
}

//! Returns an integer representing the species index
size_t STITbyPDSS::speciesIndex() const
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
void  STITbyPDSS::updateProperties(const doublereal* tempPoly,
                                   doublereal* cp_R, doublereal* h_RT,
                                   doublereal* s_R) const
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
void  STITbyPDSS::updatePropertiesTemp(const doublereal temp,
                                       doublereal* cp_R,
                                       doublereal* h_RT,
                                       doublereal* s_R) const
{
    //m_vpssmgr_ptr->setState_T(temp);
    m_PDSS_ptr->setTemperature(temp);
    AssertThrowMsg(m_speciesIndex != npos, "STITbyPDSS::updatePropertiesTemp",
                   "object was probably not installed correctly");
    h_RT[m_speciesIndex] = m_PDSS_ptr->enthalpy_RT_ref();
    cp_R[m_speciesIndex] = m_PDSS_ptr->cp_R_ref();
    s_R[m_speciesIndex]  = m_PDSS_ptr->entropy_R_ref();
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
void  STITbyPDSS::reportParameters(size_t& index, int& type,
                                   doublereal& minTemp, doublereal& maxTemp,
                                   doublereal& refPressure,
                                   doublereal* const coeffs) const
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
void  STITbyPDSS::modifyParameters(doublereal* coeffs)
{
}


}
