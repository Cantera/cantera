/**
 * @file MetalSHEelectrons.cpp
 * Definition file for the MetalSHEElectrons class, which represents the
 * electrons in a metal that are consistent with the
 * SHE electrode (see \ref thermoprops and
 * class \link Cantera::MetalSHEelectrons MetalSHEelectrons\endlink)
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/base/ctml.h"
#include "cantera/thermo/MetalSHEelectrons.h"
#include "cantera/thermo/ThermoFactory.h"

namespace Cantera
{

// ----  Constructors -------

MetalSHEelectrons::MetalSHEelectrons()
{
    warn_deprecated("Class MetalSHEelectrons", "To be removed after Cantera 2.4");
}

MetalSHEelectrons::MetalSHEelectrons(const std::string& infile, const std::string& id_)
{
    warn_deprecated("Class MetalSHEelectrons", "To be removed after Cantera 2.4");
    initThermoFile(infile, id_);
}

MetalSHEelectrons::MetalSHEelectrons(XML_Node& xmlphase, const std::string& id_)
{
    warn_deprecated("Class MetalSHEelectrons", "To be removed after Cantera 2.4");
    importPhase(xmlphase, this);
}

// ----- Mechanical Equation of State ------

double MetalSHEelectrons::pressure() const
{
    return m_press;
}

void MetalSHEelectrons::setPressure(double p)
{
    m_press = p;
}

double MetalSHEelectrons::isothermalCompressibility() const
{
    return 1.0/pressure();
}

double MetalSHEelectrons::thermalExpansionCoeff() const
{
    return 1.0/temperature();
}

// ---- Chemical Potentials and Activities ----

void MetalSHEelectrons::getActivityConcentrations(double* c) const
{
    c[0] = 1.0;
}

double MetalSHEelectrons::standardConcentration(size_t k) const
{
    return 1.0;
}

double MetalSHEelectrons::logStandardConc(size_t k) const
{
    return 0.0;
}

// Properties of the Standard State of the Species in the Solution

void MetalSHEelectrons::getStandardChemPotentials(double* mu0) const
{
    getGibbs_RT(mu0);
    mu0[0] *= RT();
}

void MetalSHEelectrons::getEnthalpy_RT(double* hrt) const
{
    getEnthalpy_RT_ref(hrt);
}

void MetalSHEelectrons::getEntropy_R(double* sr) const
{
    getEntropy_R_ref(sr);
    double tmp = log(pressure() / m_p0);
    sr[0] -= tmp;
}

void MetalSHEelectrons::getGibbs_RT(double* grt) const
{
    getGibbs_RT_ref(grt);
    double tmp = log(pressure() / m_p0);
    grt[0] += tmp;
}

void MetalSHEelectrons::getCp_R(double* cpr) const
{
    _updateThermo();
    cpr[0] = m_cp0_R;
}
void MetalSHEelectrons::getIntEnergy_RT(double* urt) const
{
    getEnthalpy_RT(urt);
    urt[0] -= 1.0;
}

void MetalSHEelectrons::getIntEnergy_RT_ref(double* urt) const
{
    _updateThermo();
    urt[0] = m_h0_RT - m_p0 / molarDensity() / RT();
}

// ---- Initialization and Internal functions

void MetalSHEelectrons::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    // Find the Thermo XML node
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError("MetalSHEelectrons::initThermoXML",
                           "no thermo XML node");
    }
    XML_Node& tnode = phaseNode.child("thermo");
    double dens = 2.65E3;
    if (tnode.hasChild("density")) {
        dens = getFloat(tnode, "density", "toSI");
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

void MetalSHEelectrons::setParameters(int n, double* const c)
{
    setDensity(c[0]);
}

void MetalSHEelectrons::getParameters(int& n, double* const c) const
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
    double rho = 2.65E3;
    if (eosdata.hasChild("density")) {
        rho = getFloat(eosdata, "density", "toSI");
    }
    setDensity(rho);
}

}
