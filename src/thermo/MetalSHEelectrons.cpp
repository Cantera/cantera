/**
 * @file MetalSHEelectrons.cpp
 * Definition file for the %MetalSHEElectrons class, which represents the
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
#include "cantera/base/ct_defs.h"

#include "cantera/thermo/MetalSHEelectrons.h"
#include "cantera/thermo/SingleSpeciesTP.h"
#include "cantera/thermo/ThermoFactory.h"

#include <string>

namespace Cantera
{

/*
 * ----  Constructors -------
 */
//====================================================================================================================
/*
 * Default Constructor for the MetalSHEelectrons class
 */
MetalSHEelectrons::MetalSHEelectrons():
    SingleSpeciesTP(),
    xdef_(0)
{
}
//====================================================================================================================
// Create and initialize a MetalSHEelectrons ThermoPhase object
// from an ASCII input file
/*
 * @param infile name of the input file
 * @param id     name of the phase id in the file.
 *               If this is blank, the first phase in the file is used.
 */
MetalSHEelectrons::MetalSHEelectrons(const std::string& infile, std::string id) :
    SingleSpeciesTP(),
    xdef_(0)
{
    XML_Node* root;
    if (infile == "MetalSHEelectrons_default.xml") {
        xdef_ = MetalSHEelectrons::makeDefaultXMLTree();
        root = xdef_;
    } else {
        root = get_XML_File(infile);
    }
    if (id == "-") {
        id = "";
    }
    XML_Node* xphase = get_XML_NameID("phase", std::string("#")+id, root);
    if (!xphase) {
        throw CanteraError("MetalSHEelectrons::MetalSHEelectrons",
                           "Couldn't find phase name in file:" + id);
    }
    // Check the model name to ensure we have compatibility
    const XML_Node& th = xphase->child("thermo");
    std::string model = th["model"];
    if (model != "MetalSHEelectrons") {
        throw CanteraError("MetalSHEelectrons::MetalSHEelectrons",
                           "thermo model attribute must be MetalSHEelectrons");
    }
    importPhase(*xphase, this);
}
//====================================================================================================================
// Full Constructor.
/*
 *  @param phaseRef XML node pointing to a MetalSHEelectrons description
 *  @param id       Id of the phase.
 */
MetalSHEelectrons::MetalSHEelectrons(XML_Node& xmlphase, const std::string& id) :
    SingleSpeciesTP(),
    xdef_(0)
{
    if (id != "") {
        std::string idxml = xmlphase["id"];
        if (id != idxml) {
            throw CanteraError("MetalSHEelectrons::MetalSHEelectrons",
                               "id's don't match");
        }
    }
    const XML_Node& th = xmlphase.child("thermo");
    std::string model = th["model"];
    if (model != "MetalSHEelectrons") {
        throw CanteraError("MetalSHEelectrons::MetalSHEelectrons",
                           "thermo model attribute must be MetalSHEelectrons");
    }
    importPhase(xmlphase, this);
}
//====================================================================================================================
// Copy constructor
/*
 * @param right Object to be copied
 */
MetalSHEelectrons::MetalSHEelectrons(const MetalSHEelectrons&  right) :
    SingleSpeciesTP()
{
    operator=(right);
}
//====================================================================================================================
/*
 * Destructor for the routine (virtual)
 *
 */
MetalSHEelectrons::~MetalSHEelectrons()
{
    if (xdef_) {
        delete xdef_;
    }
}
//====================================================================================================================
// Assignment operator
/*
 * @param right Object to be copied
 */
MetalSHEelectrons&
MetalSHEelectrons::operator=(const MetalSHEelectrons& right)
{
    if (&right != this) {
        SingleSpeciesTP::operator=(right);
    }

    if (xdef_) {
        delete xdef_;
    }
    xdef_ = new XML_Node(*right.xdef_);

    return *this;
}
//====================================================================================================================
// Duplication function
/*
 * This virtual function is used to create a duplicate of the
 * current phase. It's used to duplicate the phase when given
 * a ThermoPhase pointer to the phase.
 *
 * @return It returns a ThermoPhase pointer.
 */
ThermoPhase* MetalSHEelectrons::duplMyselfAsThermoPhase() const
{
    MetalSHEelectrons* stp = new MetalSHEelectrons(*this);
    return (ThermoPhase*) stp;
}
//====================================================================================================================

/*
 * ---- Utilities -----
 */

/*
 * Equation of state flag. Returns the value cStoichSubstance,
 * defined in mix_defs.h.
 */
int MetalSHEelectrons::eosType() const
{
    return cMetalSHEelectrons;
}
//====================================================================================================================

/*
 * ---- Molar Thermodynamic properties of the solution ----
 */

/**
 * ----- Mechanical Equation of State ------
 */
//====================================================================================================================
/*
 * Pressure. Units: Pa.
 * For an incompressible substance, the density is independent
 * of pressure. This method simply returns the stored
 * pressure value.
 */
doublereal MetalSHEelectrons::pressure() const
{
    return m_press;
}
//====================================================================================================================
/*
 * Set the pressure at constant temperature. Units: Pa.
 * For an incompressible substance, the density is
 * independent of pressure. Therefore, this method only
 * stores the specified pressure value. It does not
 * modify the density.
 */
void MetalSHEelectrons::setPressure(doublereal p)
{
    m_press = p;
}
//====================================================================================================================
/*
 * The isothermal compressibility. Units: 1/Pa.
 * The isothermal compressibility is defined as
 * \f[
 * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
 * \f]
 *
 *  It's equal to zero for this model, since the molar volume
 *  doesn't change with pressure or temperature.
 */
doublereal MetalSHEelectrons::isothermalCompressibility() const
{
    return -1.0/pressure();
}
//====================================================================================================================
/*
 * The thermal expansion coefficient. Units: 1/K.
 * The thermal expansion coefficient is defined as
 *
 * \f[
 * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
 * \f]
 *
 *  It's equal to zero for this model, since the molar volume
 *  doesn't change with pressure or temperature.
 */
doublereal MetalSHEelectrons::thermalExpansionCoeff() const
{
    return 1.0/temperature();

}
//====================================================================================================================
/*
 * ---- Chemical Potentials and Activities ----
 */
//====================================================================================================================
/*
 * This method returns the array of generalized
 * concentrations.  For a stoichiometric substance, there is
 * only one species, and the generalized concentration is 1.0.
 */
void MetalSHEelectrons::
getActivityConcentrations(doublereal* c) const
{
    c[0] = 1.0;
}
//====================================================================================================================
/*
 * The standard concentration. This is defined as the concentration
 * by which the generalized concentration is normalized to produce
 * the activity.
 */
doublereal MetalSHEelectrons::standardConcentration(size_t k) const
{
    return 1.0;
}
//====================================================================================================================
/*
 * Returns the natural logarithm of the standard
 * concentration of the kth species
 */
doublereal MetalSHEelectrons::logStandardConc(size_t k) const
{
    return 0.0;
}
//====================================================================================================================
/*
 * Returns the units of the standard and generalized
 * concentrations Note they have the same units, as their
 * ratio is defined to be equal to the activity of the kth
 * species in the solution, which is unitless.
 *
 * This routine is used in print out applications where the
 * units are needed. Usually, MKS units are assumed throughout
 * the program and in the XML input files.
 *
 *  uA[0] = kmol units - default  = 1
 *  uA[1] = m    units - default  = -nDim(), the number of spatial
 *                                dimensions in the Phase class.
 *  uA[2] = kg   units - default  = 0;
 *  uA[3] = Pa(pressure) units - default = 0;
 *  uA[4] = Temperature units - default = 0;
 *  uA[5] = time units - default = 0
 */
void MetalSHEelectrons::
getUnitsStandardConc(doublereal* uA, int k, int sizeUA) const
{
    for (int i = 0; i < 6; i++) {
        uA[i] = 0;
    }
}
//====================================================================================================================
/*
 *  ---- Partial Molar Properties of the Solution ----
 */

//====================================================================================================================

/*
 * ---- Properties of the Standard State of the Species in the Solution
 * ----
 */
//====================================================================================================================
/*
 * Get the array of chemical potentials at unit activity
 * \f$ \mu^0_k \f$.
 *
 * For a stoichiometric substance, there is no activity term in
 * the chemical potential expression, and therefore the
 * standard chemical potential and the chemical potential
 * are both equal to the molar Gibbs function.
 */
void MetalSHEelectrons::
getStandardChemPotentials(doublereal* mu0) const
{
    getGibbs_RT(mu0);
    mu0[0] *= GasConstant * temperature();
}
//====================================================================================================================
/*
 * Get the nondimensional Enthalpy functions for the species
 * at their standard states at the current
 * <I>T</I> and <I>P</I> of the solution.
 * Molar enthalpy. Units: J/kmol.  For an incompressible,
 * stoichiometric substance, the internal energy is
 * independent of pressure, and therefore the molar enthalpy
 * is \f[ \hat h(T, P) = \hat u(T) + P \hat v \f], where the
 * molar specific volume is constant.
 */
void MetalSHEelectrons::getEnthalpy_RT(doublereal* hrt) const
{
    getEnthalpy_RT_ref(hrt);
}
//====================================================================================================================
/*
 * Get the array of nondimensional Entropy functions for the
 * standard state species
 * at the current <I>T</I> and <I>P</I> of the solution.
 */
void MetalSHEelectrons::getEntropy_R(doublereal* sr) const
{
    getEntropy_R_ref(sr);
    doublereal tmp = log(pressure() / m_p0);
    sr[0] -= tmp;
}
//====================================================================================================================
/*
 * Get the nondimensional Gibbs functions for the species
 * at their standard states of solution at the current T and P
 * of the solution
 */
void MetalSHEelectrons::getGibbs_RT(doublereal* grt) const
{
    getGibbs_RT_ref(grt);
    doublereal tmp = log(pressure() / m_p0);
    grt[0] += tmp;
}
//====================================================================================================================
/*
 * Get the nondimensional Gibbs functions for the standard
 * state of the species at the current T and P.
 */
void MetalSHEelectrons::getCp_R(doublereal* cpr) const
{
    _updateThermo();
    cpr[0] = m_cp0_R[0];
}
//====================================================================================================================
/*
 * Molar internal energy (J/kmol).
 * For an incompressible,
 * stoichiometric substance, the molar internal energy is
 * independent of pressure. Since the thermodynamic properties
 * are specified by giving the standard-state enthalpy, the
 * term \f$ P_0 \hat v\f$ is subtracted from the specified molar
 * enthalpy to compute the molar internal energy.
 */
void MetalSHEelectrons::getIntEnergy_RT(doublereal* urt) const
{
    getEnthalpy_RT(urt);
    urt[0] -= 1.0;
}
//====================================================================================================================
/*
 * ---- Thermodynamic Values for the Species Reference States ----
 */
/*
 * Molar internal energy or the reference state at the current
 * temperature, T  (J/kmol).
 * For an incompressible,
 * stoichiometric substance, the molar internal energy is
 * independent of pressure. Since the thermodynamic properties
 * are specified by giving the standard-state enthalpy, the
 * term \f$ P_0 \hat v\f$ is subtracted from the specified molar
 * enthalpy to compute the molar internal energy.
 *
 * Note, this is equal to the standard state internal energy
 * evaluated at the reference pressure.
 */
void MetalSHEelectrons::getIntEnergy_RT_ref(doublereal* urt) const
{
    _updateThermo();
    doublereal RT = GasConstant * temperature();
    doublereal PV = m_p0 / molarDensity();
    urt[0] = m_h0_RT[0] - PV / RT;
}

/*
 * ---- Saturation Properties
 */



/*
 * ---- Initialization and Internal functions
 */
//====================================================================================================================
/*
 * @internal Initialize. This method is provided to allow
 * subclasses to perform any initialization required after all
 * species have been added. For example, it might be used to
 * resize internal work arrays that must have an entry for
 * each species.  The base class implementation does nothing,
 * and subclasses that do not require initialization do not
 * need to overload this method.  When importing a CTML phase
 * description, this method is called just prior to returning
 * from function importPhase.
 *
 * @see importCTML.cpp
 */
void MetalSHEelectrons::initThermo()
{
    /*
     * Call the base class thermo initializer
     */
    SingleSpeciesTP::initThermo();
}
//====================================================================================================================

void MetalSHEelectrons::initThermoXML(XML_Node& phaseNode, const std::string& id)
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
        dens = ctml::getFloatDefaultUnits(tnode, "density", "kg/m3");
    }
    setDensity(dens);
    SingleSpeciesTP::initThermoXML(phaseNode, id);
}
//====================================================================================================================
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
//====================================================================================================================
/*
 * setParameters:
 *
 *   Generic routine that is used to set the parameters used
 *   by this model.
 *        C[0] = density of phase [ kg/m3 ]
 */
void MetalSHEelectrons::setParameters(int n, doublereal* const c)
{
    doublereal rho = c[0];
    setDensity(rho);
}
//====================================================================================================================
/*
 * getParameters:
 *
 *   Generic routine that is used to get the parameters used
 *   by this model.
 *        n = 1
 *        C[0] = density of phase [ kg/m3 ]
 */
void MetalSHEelectrons::getParameters(int& n, doublereal* const c) const
{
    doublereal rho = density();
    n = 1;
    c[0] = rho;
}
//====================================================================================================================
/*
 * Reads an xml data block for the parameters needed by this
 * routine. eosdata is a reference to the xml thermo block, and looks
 * like this:
 *
 *   <phase id="stoichsolid" >
 *     <thermo model="StoichSubstance">
 *         <density units="g/cm3">3.52</density>
 *     </thermo>
 *   </phase>
 */
void MetalSHEelectrons::setParametersFromXML(const XML_Node& eosdata)
{
    std::string model = eosdata["model"];
    if (model != "MetalSHEelectrons") {
        throw CanteraError("MetalSHEelectrons::setParametersFromXML",
                           "thermo model attribute must be MetalSHEelectrons");
    }
    doublereal rho = 2.65E3;
    if (eosdata.hasChild("density")) {
        rho = ctml::getFloat(eosdata, "density", "toSI");
    }
    setDensity(rho);
}
//====================================================================================================================

}
