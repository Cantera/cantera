/**
 * @file PDSS_IdealGas.cpp
 * Implementation of a pressure dependent standard state
 * virtual function.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/base/ctml.h"
#include "cantera/thermo/PDSS_IdealGas.h"
#include "cantera/thermo/VPStandardStateTP.h"

using namespace std;

namespace Cantera
{
PDSS_IdealGas::PDSS_IdealGas(VPStandardStateTP* tp, int spindex) :
    PDSS(tp, spindex)
{
    m_pdssType = cPDSS_IDEALGAS;
}

PDSS_IdealGas::PDSS_IdealGas(VPStandardStateTP* tp, int spindex,
                             const std::string& inputFile, const std::string& id) :
    PDSS(tp, spindex)
{
    warn_deprecated("PDSS_IdealGas constructor from XML input file",
                    "To be removed after Cantera 2.3.");
    m_pdssType = cPDSS_IDEALGAS;
    constructPDSSFile(tp, spindex, inputFile, id);
}

PDSS_IdealGas::PDSS_IdealGas(VPStandardStateTP* tp, size_t spindex, const XML_Node& speciesNode,
                             const XML_Node& phaseRoot, bool spInstalled) :
    PDSS(tp, spindex)
{
    if (!spInstalled) {
        throw CanteraError("PDSS_IdealGas", "sp installing not done yet");
    }
    m_pdssType = cPDSS_IDEALGAS;
    constructPDSSXML(tp, spindex, phaseRoot, "");
}

PDSS_IdealGas::PDSS_IdealGas(const PDSS_IdealGas& b) :
    PDSS(b)
{
    // Use the assignment operator to do the brunt of the work for the copy
    // constructor.
    *this = b;
}

PDSS_IdealGas& PDSS_IdealGas::operator=(const PDSS_IdealGas& b)
{
    if (&b == this) {
        return *this;
    }
    PDSS::operator=(b);
    return *this;
}

PDSS* PDSS_IdealGas::duplMyselfAsPDSS() const
{
    return new PDSS_IdealGas(*this);
}

void PDSS_IdealGas::constructPDSSXML(VPStandardStateTP* tp, size_t spindex,
                                     const XML_Node& phaseNode, const std::string& id)
{
}

void PDSS_IdealGas::constructPDSSFile(VPStandardStateTP* tp, size_t spindex,
                                      const std::string& inputFile,
                                      const std::string& id)
{
    warn_deprecated("PDSS_IdealGas::constructPDSSFile",
                    "To be removed after Cantera 2.3.");
    if (inputFile.size() == 0) {
        throw CanteraError("PDSS_IdealGas::constructPDSSFile",
                           "input file is null");
    }

    // The phase object automatically constructs an XML object. Use this object
    // to store information.
    XML_Node fxml;
    fxml.build(findInputFile(inputFile));
    XML_Node* fxml_phase = findXMLPhase(&fxml, id);
    if (!fxml_phase) {
        throw CanteraError("PDSS_IdealGas::constructPDSSFile",
                           "ERROR: Can not find phase named " +
                           id + " in file named " + inputFile);
    }
    constructPDSSXML(tp, spindex, *fxml_phase, id);
}

void PDSS_IdealGas::initThermo()
{
    PDSS::initThermo();
    m_p0 = m_tp->speciesThermo().refPressure(m_spindex);
    m_minTemp = m_spthermo->minTemp(m_spindex);
    m_maxTemp = m_spthermo->maxTemp(m_spindex);
}

doublereal PDSS_IdealGas::enthalpy_RT() const
{
    return m_h0_RT_ptr[m_spindex];
}

doublereal PDSS_IdealGas::intEnergy_mole() const
{
    return (m_h0_RT_ptr[m_spindex] - 1.0) * GasConstant * m_temp;
}

doublereal PDSS_IdealGas::entropy_R() const
{
    return m_s0_R_ptr[m_spindex] - log(m_pres/m_p0);
}

doublereal PDSS_IdealGas::gibbs_RT() const
{
    return m_g0_RT_ptr[m_spindex] + log(m_pres/m_p0);
}

doublereal PDSS_IdealGas::cp_R() const
{
    return m_cp0_R_ptr[m_spindex];
}

doublereal PDSS_IdealGas::molarVolume() const
{
    return GasConstant * m_temp / m_pres;
}

doublereal PDSS_IdealGas::density() const
{
    return m_pres * m_mw / (GasConstant * m_temp);
}

doublereal PDSS_IdealGas::cv_mole() const
{
    return cp_mole() - GasConstant;
}

doublereal PDSS_IdealGas::gibbs_RT_ref() const
{
    return m_g0_RT_ptr[m_spindex];
}

doublereal PDSS_IdealGas::enthalpy_RT_ref() const
{
    return m_h0_RT_ptr[m_spindex];
}

doublereal PDSS_IdealGas::entropy_R_ref() const
{
    return m_s0_R_ptr[m_spindex];
}

doublereal PDSS_IdealGas::cp_R_ref() const
{
    return cp_R();
}

doublereal PDSS_IdealGas::molarVolume_ref() const
{
    return GasConstant * m_temp / m_p0;
}

doublereal PDSS_IdealGas::pressure() const
{
    throw CanteraError("PDSS_IdealGas::pressure()", "unimplemented");
}

void PDSS_IdealGas::setPressure(doublereal p)
{
    m_sss_R_ptr[m_spindex] = m_s0_R_ptr[m_spindex] + log(m_pres/m_p0);
    m_gss_RT_ptr[m_spindex] = m_hss_RT_ptr[m_spindex] - m_sss_R_ptr[m_spindex];
    m_Vss_ptr[m_spindex] = GasConstant * m_temp / m_pres;
}

doublereal PDSS_IdealGas::temperature() const
{
    m_temp = m_vpssmgr_ptr->temperature();
    return m_temp;
}

void PDSS_IdealGas::setTemperature(doublereal temp)
{
    m_temp = temp;
    m_spthermo->update_one(m_spindex, temp,
                           m_cp0_R_ptr, m_h0_RT_ptr, m_s0_R_ptr);
    m_g0_RT_ptr[m_spindex] = m_h0_RT_ptr[m_spindex] - m_s0_R_ptr[m_spindex];
    m_V0_ptr[m_spindex] = GasConstant * m_temp / m_p0;
    m_hss_RT_ptr[m_spindex] = m_h0_RT_ptr[m_spindex];
    m_cpss_R_ptr[m_spindex] = m_cp0_R_ptr[m_spindex];
    m_sss_R_ptr[m_spindex] = m_s0_R_ptr[m_spindex] + log(m_pres/m_p0);
    m_gss_RT_ptr[m_spindex] = m_hss_RT_ptr[m_spindex] - m_sss_R_ptr[m_spindex];
    m_Vss_ptr[m_spindex] = GasConstant * m_temp / m_pres;
}

void PDSS_IdealGas::setState_TP(doublereal temp, doublereal pres)
{
    m_pres = pres;
    setTemperature(temp);
}

void PDSS_IdealGas::setState_TR(doublereal temp, doublereal rho)
{
    m_pres = GasConstant * temp * rho / m_mw;
    setTemperature(temp);
}

}
