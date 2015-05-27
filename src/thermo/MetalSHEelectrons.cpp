/**
 * @file MetalSHEelectrons.cpp
 * Definition file for the MetalSHEElectrons class, which represents the
 * electrons in a metal that are consistent with the
 * SHE electrode (see \ref thermoprops and
 * class \link Cantera::MetalSHEelectrons MetalSHEelectrons\endlink)
 */

/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 *
 */
#include "cantera/thermo/mix_defs.h"
#include "cantera/base/ctml.h"
#include "cantera/thermo/MetalSHEelectrons.h"
#include "cantera/thermo/ThermoFactory.h"

namespace Cantera
{

/*
 * ----  Constructors -------
 */

MetalSHEelectrons::MetalSHEelectrons():
    xdef_(0)
{
}

MetalSHEelectrons::MetalSHEelectrons(const std::string& infile, std::string id_) :
    xdef_(0)
{
    XML_Node* root;
    if (infile == "MetalSHEelectrons_default.xml") {
        xdef_ = MetalSHEelectrons::makeDefaultXMLTree();
        root = xdef_;
    } else {
        root = get_XML_File(infile);
    }
    if (id_ == "-") {
        id_ = "";
    }
    XML_Node* xphase = get_XML_NameID("phase", std::string("#")+id_, root);
    if (!xphase) {
        throw CanteraError("MetalSHEelectrons::MetalSHEelectrons",
                           "Couldn't find phase name in file:" + id_);
    }
    // Check the model name to ensure we have compatibility
    if (xphase->child("thermo")["model"] != "MetalSHEelectrons") {
        throw CanteraError("MetalSHEelectrons::MetalSHEelectrons",
                           "thermo model attribute must be MetalSHEelectrons");
    }
    importPhase(*xphase, this);
}

MetalSHEelectrons::MetalSHEelectrons(XML_Node& xmlphase, const std::string& id_) :
    xdef_(0)
{
    if (id_ != "") {
        if (id_ != xmlphase["id"]) {
            throw CanteraError("MetalSHEelectrons::MetalSHEelectrons",
                               "id's don't match");
        }
    }
    if (xmlphase.child("thermo")["model"] != "MetalSHEelectrons") {
        throw CanteraError("MetalSHEelectrons::MetalSHEelectrons",
                           "thermo model attribute must be MetalSHEelectrons");
    }
    importPhase(xmlphase, this);
}

MetalSHEelectrons::MetalSHEelectrons(const MetalSHEelectrons&  right) :
    xdef_(0)
{
    operator=(right);
}

MetalSHEelectrons::~MetalSHEelectrons()
{
    delete xdef_;
}

MetalSHEelectrons&
MetalSHEelectrons::operator=(const MetalSHEelectrons& right)
{
    if (&right != this) {
        SingleSpeciesTP::operator=(right);
    }

    delete xdef_;
    if(right.xdef_)
    {
      xdef_ = new XML_Node(*right.xdef_);
    }

    return *this;
}

ThermoPhase* MetalSHEelectrons::duplMyselfAsThermoPhase() const
{
    return new MetalSHEelectrons(*this);
}

/*
 * ---- Utilities -----
 */

int MetalSHEelectrons::eosType() const
{
    return cMetalSHEelectrons;
}

/*
 * ----- Mechanical Equation of State ------
 */

doublereal MetalSHEelectrons::pressure() const
{
    return m_press;
}

void MetalSHEelectrons::setPressure(doublereal p)
{
    m_press = p;
}

doublereal MetalSHEelectrons::isothermalCompressibility() const
{
    return 1.0/pressure();
}

doublereal MetalSHEelectrons::thermalExpansionCoeff() const
{
    return 1.0/temperature();

}

/*
 * ---- Chemical Potentials and Activities ----
 */

void MetalSHEelectrons::getActivityConcentrations(doublereal* c) const
{
    c[0] = 1.0;
}

doublereal MetalSHEelectrons::standardConcentration(size_t k) const
{
    return 1.0;
}

doublereal MetalSHEelectrons::logStandardConc(size_t k) const
{
    return 0.0;
}

void MetalSHEelectrons::getUnitsStandardConc(doublereal* uA, int k,
                                             int sizeUA) const
{
    warn_deprecated("MetalSHEelectrons::getUnitsStandardConc",
            "To be removed after Cantera 2.2.");

    for (int i = 0; i < 6; i++) {
        uA[i] = 0;
    }
}

/*
 * Properties of the Standard State of the Species in the Solution
 */

void MetalSHEelectrons::getStandardChemPotentials(doublereal* mu0) const
{
    getGibbs_RT(mu0);
    mu0[0] *= GasConstant * temperature();
}

void MetalSHEelectrons::getEnthalpy_RT(doublereal* hrt) const
{
    getEnthalpy_RT_ref(hrt);
}

void MetalSHEelectrons::getEntropy_R(doublereal* sr) const
{
    getEntropy_R_ref(sr);
    doublereal tmp = log(pressure() / m_p0);
    sr[0] -= tmp;
}

void MetalSHEelectrons::getGibbs_RT(doublereal* grt) const
{
    getGibbs_RT_ref(grt);
    doublereal tmp = log(pressure() / m_p0);
    grt[0] += tmp;
}

void MetalSHEelectrons::getCp_R(doublereal* cpr) const
{
    _updateThermo();
    cpr[0] = m_cp0_R[0];
}
void MetalSHEelectrons::getIntEnergy_RT(doublereal* urt) const
{
    getEnthalpy_RT(urt);
    urt[0] -= 1.0;
}

void MetalSHEelectrons::getIntEnergy_RT_ref(doublereal* urt) const
{
    _updateThermo();
    doublereal RT = GasConstant * temperature();
    urt[0] = m_h0_RT[0] - m_p0 / molarDensity() / RT;
}

/*
 * ---- Initialization and Internal functions
 */

void MetalSHEelectrons::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    /*
     * Find the Thermo XML node
     */
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError("MetalSHEelectrons::initThermoXML",
                           "no thermo XML node");
    }
    XML_Node& tnode = phaseNode.child("thermo");
    doublereal dens = 2.65E3;
    if (tnode.hasChild("density")) {
        dens = getFloatDefaultUnits(tnode, "density", "kg/m3");
    }
    setDensity(dens);
    SingleSpeciesTP::initThermoXML(phaseNode, id_);
}

XML_Node* MetalSHEelectrons::makeDefaultXMLTree()
{
    XML_Node* xtop = new XML_Node("ctml", 0);
    XML_Node& xv = xtop->addChild("validate");
    xv.addAttribute("reactions", "yes");
    xv.addAttribute("species", "yes");

    XML_Node& xp = xtop->addChild("phase");
    xp.addAttribute("dim", "3");
    xp.addAttribute("id", "MetalSHEelectrons");
    XML_Node& xe = xp.addChild("elementArray", "E");
    xe.addAttribute("datasrc", "elements.xml");
    XML_Node& xs = xp.addChild("speciesArray", "she_electron");
    xs.addAttribute("datasrc", "#species_Metal_SHEelectrons");
    XML_Node& xt = xp.addChild("thermo");
    xt.addAttribute("model", "metalSHEelectrons");
    XML_Node& xtr = xp.addChild("transport");
    xtr.addAttribute("model", "none");
    XML_Node& xk = xp.addChild("kinetics");
    xk.addAttribute("model", "none");

    XML_Node& xsd = xtop->addChild("speciesData");
    xsd.addAttribute("id", "species_Metal_SHEelectrons");

    XML_Node& xsp = xsd.addChild("species");
    xsp.addAttribute("name", "she_electron");
    xsp.addChild("atomArray", "E:1");
    xsp.addChild("charge", "-1");
    XML_Node& xspt = xsp.addChild("thermo");

    XML_Node& xN1 = xspt.addChild("NASA");
    xN1.addAttribute("Tmax", "1000.");
    xN1.addAttribute("Tmin", "200.");
    xN1.addAttribute("P0", "100000.0");
    XML_Node& xF1 = xsd.addChild("floatArray",
                                 "1.172165560E+00,   3.990260375E-03,  -9.739075500E-06, "
                                 "1.007860470E-08, -3.688058805E-12, -4.589675865E+02,  3.415051190E-01");
    xF1.addAttribute("name", "coeffs");
    xF1.addAttribute("size", "7");

    XML_Node& xN2 = xspt.addChild("NASA");
    xN2.addAttribute("Tmax", "6000.");
    xN2.addAttribute("Tmin", "1000.");
    xN2.addAttribute("P0", "100000.0");
    XML_Node& xF2 = xsd.addChild("floatArray",
                                 "1.466432895E+00,  4.133039835E-04, -7.320116750E-08, 7.705017950E-12,"
                                 "-3.444022160E-16, -4.065327985E+02, -5.121644350E-01");
    xF2.addAttribute("name", "coeffs");
    xF2.addAttribute("size", "7");

    return xtop;
}

void MetalSHEelectrons::setParameters(int n, doublereal* const c)
{
    setDensity(c[0]);
}

void MetalSHEelectrons::getParameters(int& n, doublereal* const c) const
{
    n = 1;
    c[0] = density();
}

void MetalSHEelectrons::setParametersFromXML(const XML_Node& eosdata)
{
    if ( eosdata["model"] != "MetalSHEelectrons") {
        throw CanteraError("MetalSHEelectrons::setParametersFromXML",
                           "thermo model attribute must be MetalSHEelectrons");
    }
    doublereal rho = 2.65E3;
    if (eosdata.hasChild("density")) {
        rho = getFloat(eosdata, "density", "toSI");
    }
    setDensity(rho);
}

}
