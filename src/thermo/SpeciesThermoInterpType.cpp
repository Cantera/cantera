/**
 *  @file SpeciesThermoInterpType.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/SpeciesThermoInterpType.h"
#include "cantera/thermo/PDSS.h"

namespace Cantera
{

SpeciesThermoInterpType::SpeciesThermoInterpType() :
    m_lowT(0.0),
    m_highT(0.0),
    m_Pref(0.0)
{
}

SpeciesThermoInterpType::SpeciesThermoInterpType(double tlow,
                                                 double thigh,
                                                 double pref) :
    m_lowT(tlow),
    m_highT(thigh),
    m_Pref(pref)
{
}

void SpeciesThermoInterpType::updateProperties(const double* tempPoly,
        double* cp_R, double* h_RT, double* s_R) const
{
    double T = tempPoly[0];
    updatePropertiesTemp(T, cp_R, h_RT, s_R);
}

double SpeciesThermoInterpType::reportHf298(double* const h298) const
{
    throw CanteraError("SpeciesThermoInterpType::reportHf298",
                       "Not implemented");
}

void SpeciesThermoInterpType::modifyOneHf298(const size_t k,
                                             const double Hf298New)
{
    throw CanteraError("SpeciesThermoInterpType::modifyOneHf298",
                       "Not implemented");
}

//=============================================================================

STITbyPDSS::STITbyPDSS(PDSS* PDSS_ptr) :
    m_PDSS_ptr(PDSS_ptr)
{
}

double STITbyPDSS::minTemp() const
{
    return m_PDSS_ptr->minTemp();
}

double STITbyPDSS::maxTemp() const
{
    return m_PDSS_ptr->maxTemp();
}

double STITbyPDSS::refPressure() const
{
    return m_PDSS_ptr->refPressure();
}

int STITbyPDSS::reportType() const
{
    return PDSS_TYPE;
}

void STITbyPDSS::updateProperties(const double* tempPoly,
                                  double* cp_R, double* h_RT,
                                  double* s_R) const
{
    double T = tempPoly[0];
    updatePropertiesTemp(T, cp_R, h_RT, s_R);
}

void STITbyPDSS::updatePropertiesTemp(const double temp,
                                      double* cp_R,
                                      double* h_RT,
                                      double* s_R) const
{
    m_PDSS_ptr->setTemperature(temp);
    *h_RT = m_PDSS_ptr->enthalpy_RT_ref();
    *cp_R = m_PDSS_ptr->cp_R_ref();
    *s_R = m_PDSS_ptr->entropy_R_ref();
}

void STITbyPDSS::reportParameters(size_t& index, int& type,
                                  double& minTemp, double& maxTemp,
                                  double& refPressure,
                                  double* const coeffs) const
{
    index = 0;
    type = PDSS_TYPE;
    minTemp = m_PDSS_ptr->minTemp();
    maxTemp = m_PDSS_ptr->maxTemp();
    refPressure = m_PDSS_ptr->refPressure();
}

}
