/**
 *  @file SpeciesThermoInterpType.cpp
 */
// Copyright 2007  Sandia National Laboratories

#include "cantera/thermo/SpeciesThermoInterpType.h"
#include "cantera/thermo/VPSSMgr.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/base/ctexceptions.h"

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

SpeciesThermoInterpType::SpeciesThermoInterpType(const SpeciesThermoInterpType &b) :
    m_lowT(b.m_lowT),
    m_highT(b.m_highT),
    m_Pref(b.m_Pref)
{
}

void SpeciesThermoInterpType::updateProperties(const doublereal* tempPoly,
        doublereal* cp_R, doublereal* h_RT, doublereal* s_R) const
{
    double T = tempPoly[0];
    updatePropertiesTemp(T, cp_R, h_RT, s_R);
}

doublereal SpeciesThermoInterpType::reportHf298(doublereal* const h298) const
{
    throw CanteraError("SpeciesThermoInterpType::reportHf298",
                       "Not implemented");
}

void SpeciesThermoInterpType::modifyOneHf298(const size_t k,
                                             const doublereal Hf298New)
{
    throw CanteraError("SpeciesThermoInterpType::modifyOneHf298",
                       "Not implemented");
}

//=============================================================================

STITbyPDSS::STITbyPDSS()
{
}

STITbyPDSS::STITbyPDSS(VPSSMgr* vpssmgr_ptr, PDSS* PDSS_ptr) :
    SpeciesThermoInterpType(),
    m_vpssmgr_ptr(vpssmgr_ptr),
    m_PDSS_ptr(PDSS_ptr)
{
}

STITbyPDSS::STITbyPDSS(const STITbyPDSS& b) :
    SpeciesThermoInterpType(b),
    m_vpssmgr_ptr(b.m_vpssmgr_ptr),
    m_PDSS_ptr(b.m_PDSS_ptr)
{
}

SpeciesThermoInterpType* STITbyPDSS::duplMyselfAsSpeciesThermoInterpType() const
{
    return new STITbyPDSS(*this);
}

void STITbyPDSS::initAllPtrs(size_t speciesIndex, VPSSMgr* vpssmgr_ptr,
                             PDSS* PDSS_ptr)
{
    m_vpssmgr_ptr = vpssmgr_ptr;
    m_PDSS_ptr = PDSS_ptr;
}

doublereal STITbyPDSS::minTemp() const
{
    return m_PDSS_ptr->minTemp();
}

doublereal STITbyPDSS::maxTemp() const
{
    return m_PDSS_ptr->maxTemp();
}

doublereal STITbyPDSS::refPressure() const
{
    return m_PDSS_ptr->refPressure();
}

int STITbyPDSS::reportType() const
{
    return PDSS_TYPE;
}

void STITbyPDSS::updateProperties(const doublereal* tempPoly,
                                  doublereal* cp_R, doublereal* h_RT,
                                  doublereal* s_R) const
{
    doublereal T = tempPoly[0];
    updatePropertiesTemp(T, cp_R, h_RT, s_R);
}

void STITbyPDSS::updatePropertiesTemp(const doublereal temp,
                                      doublereal* cp_R,
                                      doublereal* h_RT,
                                      doublereal* s_R) const
{
    m_PDSS_ptr->setTemperature(temp);
    *h_RT = m_PDSS_ptr->enthalpy_RT_ref();
    *cp_R = m_PDSS_ptr->cp_R_ref();
    *s_R = m_PDSS_ptr->entropy_R_ref();
}

void STITbyPDSS::reportParameters(size_t& index, int& type,
                                  doublereal& minTemp, doublereal& maxTemp,
                                  doublereal& refPressure,
                                  doublereal* const coeffs) const
{
    index = 0;
    type = PDSS_TYPE;
    minTemp = m_vpssmgr_ptr->minTemp();
    maxTemp = m_vpssmgr_ptr->maxTemp();
    refPressure = m_PDSS_ptr->refPressure();
}

}
