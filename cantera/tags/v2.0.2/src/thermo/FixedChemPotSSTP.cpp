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
 *
 */

#include "cantera/base/ct_defs.h"
#include "cantera/thermo/mix_defs.h"
#include "cantera/thermo/FixedChemPotSSTP.h"
#include "cantera/thermo/SpeciesThermo.h"
#include "cantera/thermo/ThermoFactory.h"


#include <string>
#include "cantera/thermo/SimpleThermo.h"
namespace Cantera
{
//====================================================================================================================
/*
 * ----  Constructors -------
 */
//====================================================================================================================
/*
 * Default Constructor for the FixedChemPotSSTP class
 */
FixedChemPotSSTP::FixedChemPotSSTP() :
    SingleSpeciesTP(),
    chemPot_(0.0)
{
}
//====================================================================================================================
// Create and initialize a FixedChemPotSSTP ThermoPhase object
// from an ASCII input file
/*
 * @param infile name of the input file
 * @param id     name of the phase id in the file.
 *               If this is blank, the first phase in the file is used.
 */
FixedChemPotSSTP::FixedChemPotSSTP(std::string infile, std::string id) :
    SingleSpeciesTP(),
    chemPot_(0.0)
{
    XML_Node* root = get_XML_File(infile);
    if (id == "-") {
        id = "";
    }
    XML_Node* xphase = get_XML_NameID("phase", std::string("#")+id, root);
    if (!xphase) {
        throw CanteraError("FixedChemPotSSTP::FixedChemPotSSTP",
                           "Couldn't find phase name in file:" + id);
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
//====================================================================================================================
// Full Constructor.
/*
 *  @param phaseRef XML node pointing to a FixedChemPotSSTP description
 *  @param id       Id of the phase.
 */
FixedChemPotSSTP::FixedChemPotSSTP(XML_Node& xmlphase, std::string id) :
    SingleSpeciesTP(),
    chemPot_(0.0)
{
    if (id != "") {
        std::string idxml = xmlphase["id"];
        if (id != idxml) {
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
//====================================================================================================================
FixedChemPotSSTP::FixedChemPotSSTP(std::string Ename, doublereal val) :
    SingleSpeciesTP(),
    chemPot_(0.0)
{

    std::string pname = Ename + "Fixed";
    setID(pname);
    setName(pname);
    setNDim(3);
    addUniqueElement(Ename, -12345.);
    freezeElements();
    vector_fp ecomp(nElements(), 0.0);
    ecomp[0] = 1.0;
    double chrg = 0.0;
    SpeciesThermo* spth = new SimpleThermo();
    setSpeciesThermo(spth);
    addUniqueSpecies(pname, &ecomp[0], chrg, 0.0);
    double c[4];
    c[0] = 298.15;
    c[1] = val;
    c[2] = 0.0;
    c[3] = 0.0;
    m_spthermo->install(pname, 0, SIMPLE, c, 0.0, 1.0E30, OneAtm);
    freezeSpecies();
    initThermo();
    m_p0 = OneAtm;
    m_tlast = 298.15;
    setChemicalPotential(val);

    // Create an XML_Node entry for this species
    XML_Node* s = new XML_Node("species", 0);
    s->addAttribute("name", pname);
    std::string aaS = Ename + ":1";
    s->addChild("atomArray", aaS);
    XML_Node& tt = s->addChild("thermo");
    XML_Node& ss = tt.addChild("Simple");
    ss.addAttribute("Pref", "1 bar");
    ss.addAttribute("Tmax", "5000.");
    ss.addAttribute("Tmin", "100.");
    ss.addChild("t0", "298.15");
    ss.addChild("cp0", "0.0");
    std::string sval = fp2str(val);
    ss.addChild("h", sval);
    ss.addChild("s", "0.0");
    saveSpeciesData(0, s);
    delete s;
    s = 0;
}

//====================================================================================================================
// Copy constructor
/*
 * @param right Object to be copied
 */
FixedChemPotSSTP::FixedChemPotSSTP(const FixedChemPotSSTP&  right) :
    SingleSpeciesTP()
{
    *this = operator=(right);
}
//====================================================================================================================
// Assignment operator
/*
 * @param right Object to be copied
 */
FixedChemPotSSTP&
FixedChemPotSSTP::operator=(const FixedChemPotSSTP& right)
{
    if (&right != this) {
        SingleSpeciesTP::operator=(right);

        chemPot_ = right.chemPot_;
    }
    return *this;
}
//====================================================================================================================
/*
 * Destructor for the routine (virtual)
 *
 */
FixedChemPotSSTP::~FixedChemPotSSTP()
{
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
ThermoPhase* FixedChemPotSSTP::duplMyselfAsThermoPhase() const
{
    FixedChemPotSSTP* stp = new FixedChemPotSSTP(*this);
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
int FixedChemPotSSTP::eosType() const
{
    return cFixedChemPot;
}

/*
 * ---- Molar Thermodynamic properties of the solution ----
 */

/*
 * ----- Mechanical Equation of State ------
 */
//====================================================================================================================
/*
 * Pressure. Units: Pa.
 * For an incompressible substance, the density is independent
 * of pressure. This method simply returns the stored
 * pressure value.
 */
doublereal FixedChemPotSSTP::pressure() const
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
void FixedChemPotSSTP::setPressure(doublereal p)
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
doublereal FixedChemPotSSTP::isothermalCompressibility() const
{
    return 0.0;
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
doublereal FixedChemPotSSTP::thermalExpansionCoeff() const
{
    return 0.0;
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
void FixedChemPotSSTP::
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
doublereal FixedChemPotSSTP::standardConcentration(size_t k) const
{
    return 1.0;
}
//====================================================================================================================
/*
 * Returns the natural logarithm of the standard
 * concentration of the kth species
 */
doublereal FixedChemPotSSTP::logStandardConc(size_t k) const
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
void FixedChemPotSSTP::
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
void FixedChemPotSSTP::getPartialMolarVolumes(doublereal* vbar) const
{
    vbar[0] = 0.0;
}
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
void FixedChemPotSSTP::
getStandardChemPotentials(doublereal* mu0) const
{
    mu0[0] = chemPot_;
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
void FixedChemPotSSTP::getEnthalpy_RT(doublereal* hrt) const
{
    double rt = _RT();
    hrt[0] = chemPot_ / rt;
}
//====================================================================================================================
/*
 * Get the array of nondimensional Entropy functions for the
 * standard state species
 * at the current <I>T</I> and <I>P</I> of the solution.
 */
void FixedChemPotSSTP::getEntropy_R(doublereal* sr) const
{
    sr[0] = 0.0;
}
//====================================================================================================================
/*
 * Get the nondimensional Gibbs functions for the species
 * at their standard states of solution at the current T and P
 * of the solution
 */
void FixedChemPotSSTP::getGibbs_RT(doublereal* grt) const
{
    double rt = _RT();
    grt[0] = chemPot_ / rt;
}
//====================================================================================================================
/*
 * Get the nondimensional Gibbs functions for the standard
 * state of the species at the current T and P.
 */
void FixedChemPotSSTP::getCp_R(doublereal* cpr) const
{
    cpr[0] = 0.0;
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
void FixedChemPotSSTP::getIntEnergy_RT(doublereal* urt) const
{
    urt[0] = chemPot_;
}
//====================================================================================================================
// Get the molar volumes of each species in their standard
// states at the current  <I>T</I> and <I>P</I> of the solution.
/*
 *   units = m^3 / kmol
 *
 * We set this to zero
 *
 * @param vbar On output this contains the standard volume of the species
 *             and phase (m^3/kmol). Vector of length 1
 */
void FixedChemPotSSTP::getStandardVolumes(doublereal* vbar) const
{
    vbar[0] = 0.0;
}
//====================================================================================================================
/*
 * ---- Thermodynamic Values for the Species Reference States ----
 */
//====================================================================================================================
void FixedChemPotSSTP::getIntEnergy_RT_ref(doublereal* urt) const
{
    urt[0] = chemPot_;
}
//====================================================================================================================
void FixedChemPotSSTP::getEnthalpy_RT_ref(doublereal* hrt) const
{
    double rt = _RT();
    hrt[0] = chemPot_ / rt;
}
//====================================================================================================================
void FixedChemPotSSTP::getEntropy_R_ref(doublereal* sr) const
{
    sr[0] = 0.0;
}
//====================================================================================================================
void FixedChemPotSSTP::getGibbs_RT_ref(doublereal* grt) const
{
    double rt = _RT();
    grt[0] = chemPot_ / rt;
}
//====================================================================================================================
void FixedChemPotSSTP::getGibbs_ref(doublereal* g) const
{
    g[0] = chemPot_;
}
//====================================================================================================================
void FixedChemPotSSTP::getCp_R_ref(doublereal* cpr) const
{
    cpr[0] = 0.0;
}
//====================================================================================================================
/*
 * ---- Saturation Properties
 */
//====================================================================================================================
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
void FixedChemPotSSTP::initThermo()
{
    /*
     * Call the base class thermo initializer
     */
    SingleSpeciesTP::initThermo();
}
//====================================================================================================================

void FixedChemPotSSTP::initThermoXML(XML_Node& phaseNode, std::string id)
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
        double val = ctml::getFloatDefaultUnits(tnode, "chemicalPotential", "J/kmol");
        chemPot_ = val;
    }
    SingleSpeciesTP::initThermoXML(phaseNode, id);


}
//====================================================================================================================
/*
 * setParameters:
 *
 *   Generic routine that is used to set the parameters used
 *   by this model.
 *        C[0] = density of phase [ kg/m3 ]
 */
void FixedChemPotSSTP::setParameters(int n, doublereal* const c)
{
    chemPot_ = c[0];
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
void FixedChemPotSSTP::getParameters(int& n, doublereal* const c) const
{
    n = 1;
    c[0] = chemPot_;
}
//====================================================================================================================
void FixedChemPotSSTP::setParametersFromXML(const XML_Node& eosdata)
{
    std::string model = eosdata["model"];
    if (model != "StoichSubstance" && model != "FixedChemPot" && model != "StoichSubstanceSSTP") {
        throw CanteraError("FixedChemPotSSTP::setParametersFromXML",
                           "thermo model attribute must be FixedChemPot or StoichSubstance or StoichSubstanceSSTP");
    }
    if (model == "FixedChemPotSSTP") {
        doublereal val = ctml::getFloatDefaultUnits(eosdata, "chemicalPotential", "J/kmol");
        chemPot_ = val;
    }
}
//====================================================================================================================
// Function to set the chemical potential directly
/*
 *  @param chemPot  Value of the chemical potential (units J/kmol)
 */
void FixedChemPotSSTP::setChemicalPotential(doublereal chemPot)
{
    chemPot_ = chemPot;
}
//====================================================================================================================
}
