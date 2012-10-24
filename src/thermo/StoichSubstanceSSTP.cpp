/**
 *  @file StoichSubstanceSSTP.cpp
 * Definition file for the StoichSubstanceSSTP class, which represents a fixed-composition
 * incompressible substance (see \ref thermoprops and
 * class \link Cantera::StoichSubstanceSSTP StoichSubstanceSSTP\endlink)
 */

/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 *
 * Copyright 2001 California Institute of Technology
 */
#include "cantera/base/ct_defs.h"
#include "cantera/thermo/mix_defs.h"
#include "cantera/thermo/StoichSubstanceSSTP.h"
#include "cantera/thermo/SpeciesThermo.h"
#include "cantera/thermo/ThermoFactory.h"

#include <string>

namespace Cantera
{

/*
 * ----  Constructors -------
 */

/*
 * Default Constructor for the StoichSubstanceSSTP class
 */
StoichSubstanceSSTP::StoichSubstanceSSTP():
    SingleSpeciesTP()
{
}

// Create and initialize a StoichSubstanceSSTP ThermoPhase object
// from an ASCII input file
/*
 * @param infile name of the input file
 * @param id     name of the phase id in the file.
 *               If this is blank, the first phase in the file is used.
 */
StoichSubstanceSSTP::StoichSubstanceSSTP(const std::string& infile, std::string id) :
    SingleSpeciesTP()
{
    XML_Node* root = get_XML_File(infile);
    if (id == "-") {
        id = "";
    }
    XML_Node* xphase = get_XML_NameID("phase", std::string("#")+id, root);
    if (!xphase) {
        throw CanteraError("StoichSubstanceSSTP::StoichSubstanceSSTP",
                           "Couldn't find phase name in file:" + id);
    }
    // Check the model name to ensure we have compatibility
    const XML_Node& th = xphase->child("thermo");
    std::string model = th["model"];
    if (model != "StoichSubstance" && model != "StoichSubstanceSSTP") {
        throw CanteraError("StoichSubstanceSSTP::StoichSubstanceSSTP",
                           "thermo model attribute must be StoichSubstance");
    }
    importPhase(*xphase, this);
}

// Full Constructor.
/*
 *  @param phaseRef XML node pointing to a StoichSubstanceSSTP description
 *  @param id       Id of the phase.
 */
StoichSubstanceSSTP::StoichSubstanceSSTP(XML_Node& xmlphase, const std::string& id) :
    SingleSpeciesTP()
{
    if (id != "") {
        std::string idxml = xmlphase["id"];
        if (id != idxml) {
            throw CanteraError("StoichSubstanceSSTP::StoichSubstanceSSTP",
                               "id's don't match");
        }
    }
    const XML_Node& th = xmlphase.child("thermo");
    std::string model = th["model"];
    if (model != "StoichSubstance" && model != "StoichSubstanceSSTP") {
        throw CanteraError("StoichSubstanceSSTP::StoichSubstanceSSTP",
                           "thermo model attribute must be StoichSubstance");
    }
    importPhase(xmlphase, this);
}

//! Copy constructor
/*!
 * @param right Object to be copied
 */
StoichSubstanceSSTP::StoichSubstanceSSTP(const StoichSubstanceSSTP&  right) :
    SingleSpeciesTP()
{
    *this = operator=(right);
}

//! Assignment operator
/*!
 * @param right Object to be copied
 */
StoichSubstanceSSTP&
StoichSubstanceSSTP::operator=(const StoichSubstanceSSTP& right)
{
    if (&right != this) {
        SingleSpeciesTP::operator=(right);
    }
    return *this;
}

/*
 * Destructor for the routine (virtual)
 *
 */
StoichSubstanceSSTP::~StoichSubstanceSSTP()
{
}

// Duplication function
/*
 * This virtual function is used to create a duplicate of the
 * current phase. It's used to duplicate the phase when given
 * a ThermoPhase pointer to the phase.
 *
 * @return It returns a ThermoPhase pointer.
 */
ThermoPhase* StoichSubstanceSSTP::duplMyselfAsThermoPhase() const
{
    return new StoichSubstanceSSTP(*this);
}


/*
 * ---- Utilities -----
 */

/*
 * Equation of state flag. Returns the value cStoichSubstance,
 * defined in mix_defs.h.
 */
int StoichSubstanceSSTP::eosType() const
{
    return cStoichSubstance;
}

/*
 * ---- Molar Thermodynamic properties of the solution ----
 */

/**
 * ----- Mechanical Equation of State ------
 */

/*
 * Pressure. Units: Pa.
 * For an incompressible substance, the density is independent
 * of pressure. This method simply returns the stored
 * pressure value.
 */
doublereal StoichSubstanceSSTP::pressure() const
{
    return m_press;
}

/*
 * Set the pressure at constant temperature. Units: Pa.
 * For an incompressible substance, the density is
 * independent of pressure. Therefore, this method only
 * stores the specified pressure value. It does not
 * modify the density.
 */
void StoichSubstanceSSTP::setPressure(doublereal p)
{
    m_press = p;
}

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
doublereal StoichSubstanceSSTP::isothermalCompressibility() const
{
    return 0.0;
}

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
doublereal StoichSubstanceSSTP::thermalExpansionCoeff() const
{
    return 0.0;
}

/*
 * ---- Chemical Potentials and Activities ----
 */

/*
 * This method returns the array of generalized
 * concentrations.  For a stoichiometric substance, there is
 * only one species, and the generalized concentration is 1.0.
 */
void StoichSubstanceSSTP::
getActivityConcentrations(doublereal* c) const
{
    c[0] = 1.0;
}

/*
 * The standard concentration. This is defined as the concentration
 * by which the generalized concentration is normalized to produce
 * the activity.
 */
doublereal StoichSubstanceSSTP::standardConcentration(size_t k) const
{
    return 1.0;
}

/*
 * Returns the natural logarithm of the standard
 * concentration of the kth species
 */
doublereal StoichSubstanceSSTP::logStandardConc(size_t k) const
{
    return 0.0;
}

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
void StoichSubstanceSSTP::
getUnitsStandardConc(doublereal* uA, int k, int sizeUA) const
{
    for (int i = 0; i < 6; i++) {
        uA[i] = 0;
    }
}

/*
 *  ---- Partial Molar Properties of the Solution ----
 */



/*
 * ---- Properties of the Standard State of the Species in the Solution
 * ----
 */

/*
 * Get the array of chemical potentials at unit activity
 * \f$ \mu^0_k \f$.
 *
 * For a stoichiometric substance, there is no activity term in
 * the chemical potential expression, and therefore the
 * standard chemical potential and the chemical potential
 * are both equal to the molar Gibbs function.
 */
void StoichSubstanceSSTP::
getStandardChemPotentials(doublereal* mu0) const
{
    getGibbs_RT(mu0);
    mu0[0] *= GasConstant * temperature();
}

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
void StoichSubstanceSSTP::getEnthalpy_RT(doublereal* hrt) const
{
    getEnthalpy_RT_ref(hrt);
    doublereal RT = GasConstant * temperature();
    doublereal presCorrect = (m_press - m_p0) /  molarDensity();
    hrt[0] += presCorrect / RT;
}

/*
 * Get the array of nondimensional Entropy functions for the
 * standard state species
 * at the current <I>T</I> and <I>P</I> of the solution.
 */
void StoichSubstanceSSTP::getEntropy_R(doublereal* sr) const
{
    getEntropy_R_ref(sr);
}

/*
 * Get the nondimensional Gibbs functions for the species
 * at their standard states of solution at the current T and P
 * of the solution
 */
void StoichSubstanceSSTP::getGibbs_RT(doublereal* grt) const
{
    getEnthalpy_RT(grt);
    grt[0] -= m_s0_R[0];
}

/*
 * Get the nondimensional Gibbs functions for the standard
 * state of the species at the current T and P.
 */
void StoichSubstanceSSTP::getCp_R(doublereal* cpr) const
{
    _updateThermo();
    cpr[0] = m_cp0_R[0];
}

/*
 * Molar internal energy (J/kmol).
 * For an incompressible,
 * stoichiometric substance, the molar internal energy is
 * independent of pressure. Since the thermodynamic properties
 * are specified by giving the standard-state enthalpy, the
 * term \f$ P_0 \hat v\f$ is subtracted from the specified molar
 * enthalpy to compute the molar internal energy.
 */
void StoichSubstanceSSTP::getIntEnergy_RT(doublereal* urt) const
{
    _updateThermo();
    doublereal RT = GasConstant * temperature();
    doublereal PV = m_p0 / molarDensity();
    urt[0] = m_h0_RT[0] - PV / RT;
}

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
void StoichSubstanceSSTP::getIntEnergy_RT_ref(doublereal* urt) const
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

/**
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
void StoichSubstanceSSTP::initThermo()
{
    /*
     * Make sure there is one and only one species in this phase.
     */
    m_kk = nSpecies();
    if (m_kk != 1) {
        throw CanteraError("initThermo",
                           "stoichiometric substances may only contain one species.");
    }
    doublereal tmin = m_spthermo->minTemp();
    doublereal tmax = m_spthermo->maxTemp();
    if (tmin > 0.0) {
        m_tmin = tmin;
    }
    if (tmax > 0.0) {
        m_tmax = tmax;
    }
    /*
     * Store the reference pressure in the variables for the class.
     */
    m_p0 = refPressure();

    /*
     * Resize temporary arrays.
     */
    int leng = 1;
    m_h0_RT.resize(leng);
    m_cp0_R.resize(leng);
    m_s0_R.resize(leng);
    /*
     * Call the base class thermo initializer
     */
    SingleSpeciesTP::initThermo();
}


void StoichSubstanceSSTP::initThermoXML(XML_Node& phaseNode, const std::string& id)
{
    /*
     * Find the Thermo XML node
     */
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError("StoichSubstanceSSTP::initThermoXML",
                           "no thermo XML node");
    }
    XML_Node& tnode = phaseNode.child("thermo");
    double dens = ctml::getFloatDefaultUnits(tnode, "density", "kg/m3");
    setDensity(dens);
    SingleSpeciesTP::initThermoXML(phaseNode, id);
}

/**
 * setParameters:
 *
 *   Generic routine that is used to set the parameters used
 *   by this model.
 *        C[0] = density of phase [ kg/m3 ]
 */
void StoichSubstanceSSTP::setParameters(int n, doublereal* const c)
{
    doublereal rho = c[0];
    setDensity(rho);
}

/**
 * getParameters:
 *
 *   Generic routine that is used to get the parameters used
 *   by this model.
 *        n = 1
 *        C[0] = density of phase [ kg/m3 ]
 */
void StoichSubstanceSSTP::getParameters(int& n, doublereal* const c) const
{
    doublereal rho = density();
    n = 1;
    c[0] = rho;
}

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
void StoichSubstanceSSTP::setParametersFromXML(const XML_Node& eosdata)
{
    std::string model = eosdata["model"];
    if (model != "StoichSubstance" && model != "StoichSubstanceSSTP") {
        throw CanteraError("StoichSubstanceSSTP::setParametersFromXML",
                           "thermo model attribute must be StoichSubstance");
    }
    doublereal rho = ctml::getFloat(eosdata, "density", "toSI");
    setDensity(rho);
}





/*
 * Default Constructor for the electrodeElectron class
 */
electrodeElectron::electrodeElectron():
    StoichSubstanceSSTP()
{
}

// Create and initialize a electrodeElectron ThermoPhase object
// from an ASCII input file
/*
 * @param infile name of the input file
 * @param id     name of the phase id in the file.
 *               If this is blank, the first phase in the file is used.
 */
electrodeElectron::electrodeElectron(const std::string& infile, std::string id) :
    StoichSubstanceSSTP()
{
    XML_Node* root = get_XML_File(infile);
    if (id == "-") {
        id = "";
    }
    XML_Node* xphase = get_XML_NameID("phase", std::string("#")+id, root);
    if (!xphase) {
        throw CanteraError("electrodeElectron::electrodeElectron",
                           "Couldn't find phase name in file:" + id);
    }
    // Check the model name to ensure we have compatibility
    const XML_Node& th = xphase->child("thermo");
    std::string model = th["model"];
    if (model != "electrodeElectron") {
        throw CanteraError("electrodeElectron::electrodeElectron",
                           "thermo model attribute must be electrodeElectron");
    }
    importPhase(*xphase, this);
}

// Full Constructor.
/*
 *  @param phaseRef XML node pointing to a electrodeElectron description
 *  @param id       Id of the phase.
 */
electrodeElectron::electrodeElectron(XML_Node& xmlphase, const std::string& id) :
    StoichSubstanceSSTP()
{
    if (id != "") {
        std::string idxml = xmlphase["id"];
        if (id != idxml) {
            throw CanteraError("electrodeElectron::electrodeElectron",
                               "id's don't match");
        }
    }
    const XML_Node& th = xmlphase.child("thermo");
    std::string model = th["model"];
    if (model != "electrodeElectron") {
        throw CanteraError("electrodeElectron::electrodeElectron",
                           "thermo model attribute must be electrodeElectron");
    }
    importPhase(xmlphase, this);
}

//! Copy constructor
/*!
 * @param right Object to be copied
 */
electrodeElectron::electrodeElectron(const electrodeElectron&  right) :
    StoichSubstanceSSTP()
{
    *this = operator=(right);
}

//! Assignment operator
/*!
 * @param right Object to be copied
 */
electrodeElectron&
electrodeElectron::operator=(const electrodeElectron& right)
{
    if (&right != this) {
        StoichSubstanceSSTP::operator=(right);
    }
    return *this;
}

/*
 * Destructor for the routine (virtual)
 *
 */
electrodeElectron::~electrodeElectron()
{
}

void electrodeElectron::setParametersFromXML(const XML_Node& eosdata)
{
    std::string model = eosdata["model"];
    if (model != "electrodeElectron") {
        throw CanteraError("electrodeElectron::setParametersFromXML",
                           "thermo model attribute must be electrodeElectron");
    }
}

void electrodeElectron::initThermoXML(XML_Node& phaseNode, const std::string& id)
{
    doublereal rho = 10.0;
    setDensity(rho);
    SingleSpeciesTP::initThermoXML(phaseNode, id);
}

void electrodeElectron::setParameters(int n, doublereal* const c)
{
    doublereal rho = 10.0;
    setDensity(rho);
}

}


