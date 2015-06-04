/**
 *  @file FixedChemPotSSTP.cpp
 * Definition file for the FixedChemPotSSTP class, which represents a fixed-composition
 * incompressible substance with a constant chemical potential (see \ref thermoprops and
 * class \link Cantera::FixedChemPotSSTP FixedChemPotSSTP\endlink)
 */

/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/thermo/mix_defs.h"
#include "cantera/thermo/FixedChemPotSSTP.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/thermo/SpeciesThermoInterpType.h"
#include "cantera/base/ctml.h"
#include "cantera/base/stringUtils.h"

namespace Cantera
{
/*
 * ----  Constructors -------
 */

FixedChemPotSSTP::FixedChemPotSSTP() :
    chemPot_(0.0)
{
}

FixedChemPotSSTP::FixedChemPotSSTP(const std::string& infile, std::string id_) :
    chemPot_(0.0)
{
    XML_Node* root = get_XML_File(infile);
    if (id_ == "-") {
        id_ = "";
    }
    XML_Node* xphase = get_XML_NameID("phase", std::string("#")+id_, root);
    if (!xphase) {
        throw CanteraError("FixedChemPotSSTP::FixedChemPotSSTP",
                           "Couldn't find phase name in file:" + id_);
    }
    // Check the model name to ensure we have compatibility
    const XML_Node& th = xphase->child("thermo");
    std::string model = th["model"];
    if (model != "StoichSubstance" && model != "StoichSubstanceSSTP" && model != "FixedChemPot") {
        throw CanteraError("FixedChemPotSSTP::FixedChemPotSSTP",
                           "thermo model attribute must be FixedChemPot or StoichSubstance");
    }
    importPhase(*xphase, this);
}
FixedChemPotSSTP::FixedChemPotSSTP(XML_Node& xmlphase, const std::string& id_) :
    chemPot_(0.0)
{
    if (id_ != "") {
        std::string idxml = xmlphase["id"];
        if (id_ != idxml) {
            throw CanteraError("FixedChemPotSSTP::FixedChemPotSSTP",
                               "id's don't match");
        }
    }
    const XML_Node& th = xmlphase.child("thermo");
    std::string model = th["model"];
    if (model != "StoichSubstance" && model != "StoichSubstanceSSTP" && model != "FixedChemPotSSTP") {
        throw CanteraError("FixedChemPotSSTP::FixedChemPotSSTP",
                           "thermo model attribute must be StoichSubstance or FixedChemPot");
    }
    importPhase(xmlphase, this);

    if (model ==  "StoichSubstance" || model == "StoichSubstanceSSTP") {
        _updateThermo();
        chemPot_ = (m_h0_RT[0] - m_s0_R[0]) * GasConstant * temperature();
    }
}

FixedChemPotSSTP::FixedChemPotSSTP(const std::string& Ename, doublereal val) :
    chemPot_(0.0)
{

    std::string pname = Ename + "Fixed";
    setID(pname);
    setName(pname);
    setNDim(3);
    addElement(Ename);
    shared_ptr<Species> sp(new Species(pname, parseCompString(Ename + ":1.0")));
    double c[4] = {298.15, val, 0.0, 0.0};
    shared_ptr<SpeciesThermoInterpType> stit(
            newSpeciesThermoInterpType("const_cp", 0.1, 1e30, OneAtm, c));
    sp->thermo = stit;
    addSpecies(sp);
    initThermo();
    m_p0 = OneAtm;
    m_tlast = 298.15;
    setChemicalPotential(val);

    // Create an XML_Node entry for this species
    XML_Node s("species", 0);
    s.addAttribute("name", pname);
    std::string aaS = Ename + ":1";
    s.addChild("atomArray", aaS);
    XML_Node& tt = s.addChild("thermo");
    XML_Node& ss = tt.addChild("Simple");
    ss.addAttribute("Pref", "1 bar");
    ss.addAttribute("Tmax", "5000.");
    ss.addAttribute("Tmin", "100.");
    ss.addChild("t0", "298.15");
    ss.addChild("cp0", "0.0");
    std::string sval = fp2str(val);
    ss.addChild("h", sval);
    ss.addChild("s", "0.0");
    saveSpeciesData(0, &s);
}

FixedChemPotSSTP::FixedChemPotSSTP(const FixedChemPotSSTP&  right)
{
    *this = right;
}

FixedChemPotSSTP&
FixedChemPotSSTP::operator=(const FixedChemPotSSTP& right)
{
    if (&right != this) {
        SingleSpeciesTP::operator=(right);

        chemPot_ = right.chemPot_;
    }
    return *this;
}

ThermoPhase* FixedChemPotSSTP::duplMyselfAsThermoPhase() const
{
    return new FixedChemPotSSTP(*this);
}

/*
 * ---- Utilities -----
 */

int FixedChemPotSSTP::eosType() const
{
    return cFixedChemPot;
}

/*
 * ----- Mechanical Equation of State ------
 */

doublereal FixedChemPotSSTP::pressure() const
{
    return m_press;
}

void FixedChemPotSSTP::setPressure(doublereal p)
{
    m_press = p;
}

doublereal FixedChemPotSSTP::isothermalCompressibility() const
{
    return 0.0;
}

doublereal FixedChemPotSSTP::thermalExpansionCoeff() const
{
    return 0.0;
}

/*
 * ---- Chemical Potentials and Activities ----
 */

void FixedChemPotSSTP::getActivityConcentrations(doublereal* c) const
{
    c[0] = 1.0;
}

doublereal FixedChemPotSSTP::standardConcentration(size_t k) const
{
    return 1.0;
}

doublereal FixedChemPotSSTP::logStandardConc(size_t k) const
{
    return 0.0;
}

void FixedChemPotSSTP::getUnitsStandardConc(doublereal* uA, int k,
                                            int sizeUA) const
{
    warn_deprecated("FixedChemPotSSTP::getUnitsStandardConc",
                    "To be removed after Cantera 2.2.");
    for (int i = 0; i < 6; i++) {
        uA[i] = 0;
    }
}

/*
 *  ---- Partial Molar Properties of the Solution ----
 */

void FixedChemPotSSTP::getPartialMolarVolumes(doublereal* vbar) const
{
    vbar[0] = 0.0;
}

/*
 * Properties of the Standard State of the Species in the Solution
 */

void FixedChemPotSSTP::getStandardChemPotentials(doublereal* mu0) const
{
    mu0[0] = chemPot_;
}

void FixedChemPotSSTP::getEnthalpy_RT(doublereal* hrt) const
{
    double rt = _RT();
    hrt[0] = chemPot_ / rt;
}

void FixedChemPotSSTP::getEntropy_R(doublereal* sr) const
{
    sr[0] = 0.0;
}

void FixedChemPotSSTP::getGibbs_RT(doublereal* grt) const
{
    double rt = _RT();
    grt[0] = chemPot_ / rt;
}

void FixedChemPotSSTP::getCp_R(doublereal* cpr) const
{
    cpr[0] = 0.0;
}

void FixedChemPotSSTP::getIntEnergy_RT(doublereal* urt) const
{
    urt[0] = chemPot_;
}

void FixedChemPotSSTP::getStandardVolumes(doublereal* vbar) const
{
    vbar[0] = 0.0;
}

/*
 * ---- Thermodynamic Values for the Species Reference States ----
 */

void FixedChemPotSSTP::getIntEnergy_RT_ref(doublereal* urt) const
{
    urt[0] = chemPot_;
}

void FixedChemPotSSTP::getEnthalpy_RT_ref(doublereal* hrt) const
{
    double rt = _RT();
    hrt[0] = chemPot_ / rt;
}

void FixedChemPotSSTP::getEntropy_R_ref(doublereal* sr) const
{
    sr[0] = 0.0;
}

void FixedChemPotSSTP::getGibbs_RT_ref(doublereal* grt) const
{
    double rt = _RT();
    grt[0] = chemPot_ / rt;
}

void FixedChemPotSSTP::getGibbs_ref(doublereal* g) const
{
    g[0] = chemPot_;
}

void FixedChemPotSSTP::getCp_R_ref(doublereal* cpr) const
{
    cpr[0] = 0.0;
}

/*
 * ---- Initialization and Internal functions
 */

void FixedChemPotSSTP::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    /*
     * Find the Thermo XML node
     */
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError("FixedChemPotSSTP::initThermoXML", "no thermo XML node");
    }
    XML_Node& tnode = phaseNode.child("thermo");
    std::string model = tnode["model"];
    if (model != "StoichSubstance" && model != "FixedChemPot" && model != "StoichSubstanceSSTP") {
        throw CanteraError("FixedChemPotSSTP::initThermoXML()",
                           "thermo model attribute must be FixedChemPot or StoichSubstance or StoichSubstanceSSTP");
    }
    if (model == "FixedChemPot") {
        double val = getFloatDefaultUnits(tnode, "chemicalPotential", "J/kmol");
        chemPot_ = val;
    }
    SingleSpeciesTP::initThermoXML(phaseNode, id_);
}

void FixedChemPotSSTP::setParameters(int n, doublereal* const c)
{
    chemPot_ = c[0];
}

void FixedChemPotSSTP::getParameters(int& n, doublereal* const c) const
{
    n = 1;
    c[0] = chemPot_;
}

void FixedChemPotSSTP::setParametersFromXML(const XML_Node& eosdata)
{
    std::string model = eosdata["model"];
    if (model != "StoichSubstance" && model != "FixedChemPot" && model != "StoichSubstanceSSTP") {
        throw CanteraError("FixedChemPotSSTP::setParametersFromXML",
                           "thermo model attribute must be FixedChemPot or StoichSubstance or StoichSubstanceSSTP");
    }
    if (model == "FixedChemPotSSTP") {
        doublereal val = getFloatDefaultUnits(eosdata, "chemicalPotential", "J/kmol");
        chemPot_ = val;
    }
}

void FixedChemPotSSTP::setChemicalPotential(doublereal chemPot)
{
    chemPot_ = chemPot;
}

}
