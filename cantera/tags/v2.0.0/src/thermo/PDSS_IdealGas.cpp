/**
 * @file PDSS_IdealGas.cpp
 * Implementation of a pressure dependent standard state
 * virtual function.
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#include "cantera/base/ct_defs.h"
#include "cantera/base/xml.h"
#include "cantera/base/ctml.h"
#include "cantera/thermo/PDSS_IdealGas.h"
#include "cantera/thermo/ThermoFactory.h"

#include "cantera/thermo/VPStandardStateTP.h"

#include <fstream>

using namespace std;

namespace Cantera
{
/**
 * Basic list of constructors and duplicators
 */

PDSS_IdealGas::PDSS_IdealGas(VPStandardStateTP* tp, int spindex) :
    PDSS(tp, spindex)
{
    m_pdssType = cPDSS_IDEALGAS;
}


PDSS_IdealGas::PDSS_IdealGas(VPStandardStateTP* tp, int spindex,
                             std::string inputFile, std::string id) :
    PDSS(tp, spindex)
{
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
    std::string id = "";
    constructPDSSXML(tp, spindex, phaseRoot, id);
}



PDSS_IdealGas::PDSS_IdealGas(const PDSS_IdealGas& b) :
    PDSS(b)
{
    /*
     * Use the assignment operator to do the brunt
     * of the work for the copy constructor.
     */
    *this = b;
}

/*
 * Assignment operator
 */
PDSS_IdealGas& PDSS_IdealGas::operator=(const PDSS_IdealGas& b)
{
    if (&b == this) {
        return *this;
    }
    PDSS::operator=(b);

    m_tmin       =   b.m_tmin;
    m_tmax       =   b.m_tmax;

    return *this;
}

PDSS_IdealGas::~PDSS_IdealGas()
{
}

// Duplicator
PDSS* PDSS_IdealGas::duplMyselfAsPDSS() const
{
    PDSS_IdealGas* idg = new PDSS_IdealGas(*this);
    return (PDSS*) idg;
}



/*
 * constructPDSSXML:
 *
 * Initialization of a PDSS_IdealGas object using an
 * xml file.
 *
 * This routine is a precursor to initThermo(XML_Node*)
 * routine, which does most of the work.
 *
 * @param infile XML file containing the description of the
 *        phase
 *
 * @param id  Optional parameter identifying the name of the
 *            phase. If none is given, the first XML
 *            phase element will be used.
 */
void PDSS_IdealGas::constructPDSSXML(VPStandardStateTP* tp, size_t spindex,
                                     const XML_Node& phaseNode, std::string id)
{
    //initThermo();
    //initThermoXML(phaseNode, id);
}


void PDSS_IdealGas::constructPDSSFile(VPStandardStateTP* tp, size_t spindex,
                                      std::string inputFile, std::string id)
{

    if (inputFile.size() == 0) {
        throw CanteraError("PDSS_IdealGas::constructPDSSFile",
                           "input file is null");
    }
    std::string path = findInputFile(inputFile);
    ifstream fin(path.c_str());
    if (!fin) {
        throw CanteraError("PDSS_IdealGas::constructPDSSFile","could not open "
                           +path+" for reading.");
    }
    /*
     * The phase object automatically constructs an XML object.
     * Use this object to store information.
     */

    XML_Node* fxml = new XML_Node();
    fxml->build(fin);
    XML_Node* fxml_phase = findXMLPhase(fxml, id);
    if (!fxml_phase) {
        throw CanteraError("PDSS_IdealGas::constructPDSSFile",
                           "ERROR: Can not find phase named " +
                           id + " in file named " + inputFile);
    }
    constructPDSSXML(tp, spindex, *fxml_phase, id);
    delete fxml;
}

void PDSS_IdealGas::initThermoXML(const XML_Node& phaseNode, std::string& id)
{
    PDSS::initThermoXML(phaseNode, id);
}

void PDSS_IdealGas::initThermo()
{
    PDSS::initThermo();
    SpeciesThermo& sp = m_tp->speciesThermo();
    m_p0 = sp.refPressure(m_spindex);
    m_minTemp = m_spthermo->minTemp(m_spindex);
    m_maxTemp = m_spthermo->maxTemp(m_spindex);
}

/*
 * Return the molar enthalpy in units of J kmol-1
 */
doublereal
PDSS_IdealGas::enthalpy_mole() const
{
    doublereal val = enthalpy_RT();
    doublereal RT = GasConstant * m_temp;
    return (val * RT);
}

doublereal
PDSS_IdealGas::enthalpy_RT() const
{
    doublereal val = m_h0_RT_ptr[m_spindex];
    return (val);
}


/*
 * Calculate the internal energy in mks units of
 * J kmol-1
 */
doublereal
PDSS_IdealGas::intEnergy_mole() const
{
    doublereal val = m_h0_RT_ptr[m_spindex] - 1.0;
    doublereal RT = GasConstant * m_temp;
    return (val * RT);
}

/*
 * Calculate the entropy in mks units of
 * J kmol-1 K-1
 */
doublereal
PDSS_IdealGas::entropy_mole() const
{
    doublereal val = entropy_R();
    return (val * GasConstant);
}

doublereal
PDSS_IdealGas::entropy_R() const
{
    doublereal val = m_s0_R_ptr[m_spindex] - log(m_pres/m_p0);
    return (val);
}

/*
 * Calculate the Gibbs free energy in mks units of
 * J kmol-1 K-1.
 */
doublereal
PDSS_IdealGas::gibbs_mole() const
{
    doublereal val = gibbs_RT();
    doublereal RT = GasConstant * m_temp;
    return (val * RT);
}

doublereal
PDSS_IdealGas::gibbs_RT() const
{
    doublereal val = m_g0_RT_ptr[m_spindex] + log(m_pres/m_p0);
    return (val);
}

/*
 * Calculate the constant pressure heat capacity
 * in mks units of J kmol-1 K-1
 */
doublereal
PDSS_IdealGas::cp_mole() const
{
    doublereal val = cp_R();
    return (val * GasConstant);
}

doublereal
PDSS_IdealGas::cp_R() const
{
    doublereal val = m_cp0_R_ptr[m_spindex];
    return (val);
}

doublereal
PDSS_IdealGas::molarVolume() const
{
    return (GasConstant * m_temp / m_pres);
}


doublereal
PDSS_IdealGas::density() const
{
    return (m_pres * m_mw / (GasConstant * m_temp));
}

/*
 * Calculate the constant volume heat capacity
 * in mks units of J kmol-1 K-1
 */
doublereal
PDSS_IdealGas::cv_mole() const
{
    return (cp_mole() - GasConstant);
}


doublereal
PDSS_IdealGas::gibbs_RT_ref() const
{
    doublereal val = m_g0_RT_ptr[m_spindex];
    return (val);
}

doublereal PDSS_IdealGas::enthalpy_RT_ref() const
{
    doublereal val = m_h0_RT_ptr[m_spindex];
    return (val);
}

doublereal PDSS_IdealGas::entropy_R_ref() const
{
    doublereal val = m_s0_R_ptr[m_spindex];
    return (val);
}

doublereal PDSS_IdealGas::cp_R_ref() const
{
    return (cp_R());
}

doublereal PDSS_IdealGas::molarVolume_ref() const
{
    return (GasConstant * m_temp / m_p0);
}

/*
 * Calculate the pressure (Pascals), given the temperature and density
 *  Temperature: kelvin
 *  rho: density in kg m-3
 */
doublereal  PDSS_IdealGas::pressure() const
{
    throw CanteraError("PDSS_IdealGas::pressure()", "unimplemented");
    return (0.0);
}

void PDSS_IdealGas::setPressure(doublereal p)
{
    m_sss_R_ptr[m_spindex] = m_s0_R_ptr[m_spindex] + log(m_pres/m_p0);
    m_gss_RT_ptr[m_spindex] = m_hss_RT_ptr[m_spindex] - m_sss_R_ptr[m_spindex];
    m_Vss_ptr[m_spindex] = GasConstant * m_temp / m_pres;
}


// critical temperature
doublereal PDSS_IdealGas::critTemperature() const
{
    throw CanteraError("PDSS_IdealGas::critTemperature()", "unimplemented");
    return (0.0);
}

// critical pressure
doublereal PDSS_IdealGas::critPressure() const
{
    throw CanteraError("PDSS_IdealGas::critPressure()", "unimplemented");
    return (0.0);
}

// critical density
doublereal PDSS_IdealGas::critDensity() const
{
    throw CanteraError("PDSS_IdealGas::critDensity()", "unimplemented");
    return (0.0);
}


/*
 * Return the temperature
 *
 * Obtain the temperature from the owning VPStandardStateTP object
 * if you can.
 */
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
    m_g0_RT_ptr[m_spindex] =  m_h0_RT_ptr[m_spindex] -  m_s0_R_ptr[m_spindex];
    m_V0_ptr[m_spindex] = GasConstant * m_temp / m_p0;

    m_hss_RT_ptr[m_spindex]  = m_h0_RT_ptr[m_spindex];
    m_cpss_R_ptr[m_spindex]  = m_cp0_R_ptr[m_spindex];
    m_sss_R_ptr[m_spindex] = m_s0_R_ptr[m_spindex] + log(m_pres/m_p0);
    m_gss_RT_ptr[m_spindex] = m_hss_RT_ptr[m_spindex] - m_sss_R_ptr[m_spindex];
    m_Vss_ptr[m_spindex] = GasConstant * m_temp / m_pres;
}


void PDSS_IdealGas::setState_TP(doublereal temp, doublereal pres)
{
    m_pres = pres;
    setTemperature(temp);
}

void  PDSS_IdealGas::setState_TR(doublereal temp, doublereal rho)
{
    m_pres = GasConstant * temp * rho / m_mw;
    setTemperature(temp);
}

// saturation pressure
doublereal PDSS_IdealGas::satPressure(doublereal t)
{
    throw CanteraError("PDSS_IdealGas::satPressure()", "unimplemented");
    /*NOTREACHED*/
    return (0.0);
}


}
