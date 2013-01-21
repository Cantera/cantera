/**
 *  @file MargulesVPSSTP.cpp
 *   Definitions for ThermoPhase object for phases which
 *   employ excess gibbs free energy formulations related to Margules
 *   expansions (see \ref thermoprops
 *    and class \link Cantera::MargulesVPSSTP MargulesVPSSTP\endlink).
 *
 */
/*
 * Copyright (2009) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/thermo/MixedSolventElectrolyte.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"

#include <iomanip>
#include <fstream>

using namespace std;

namespace Cantera
{

/*
 * Default constructor.
 *
 */
MixedSolventElectrolyte::MixedSolventElectrolyte() :
    MolarityIonicVPSSTP(),
    numBinaryInteractions_(0),
    formMargules_(0),
    formTempModel_(0)
{
}

/*
 * Working constructors
 *
 *  The two constructors below are the normal way
 *  the phase initializes itself. They are shells that call
 *  the routine initThermo(), with a reference to the
 *  XML database to get the info for the phase.

 */
MixedSolventElectrolyte::MixedSolventElectrolyte(const std::string& inputFile,
                                                 const std::string& id) :
    MolarityIonicVPSSTP(),
    numBinaryInteractions_(0),
    formMargules_(0),
    formTempModel_(0)
{
    initThermoFile(inputFile, id);
}

MixedSolventElectrolyte::MixedSolventElectrolyte(XML_Node& phaseRoot,
                                                 const std::string& id) :
    MolarityIonicVPSSTP(),
    numBinaryInteractions_(0),
    formMargules_(0),
    formTempModel_(0)
{
    importPhase(*findXMLPhase(&phaseRoot, id), this);
}


/*
 * Copy Constructor:
 *
 *  Note this stuff will not work until the underlying phase
 *  has a working copy constructor
 */
MixedSolventElectrolyte::MixedSolventElectrolyte(const MixedSolventElectrolyte& b) :
    MolarityIonicVPSSTP()
{
    MixedSolventElectrolyte::operator=(b);
}

/*
 * operator=()
 *
 *  Note this stuff will not work until the underlying phase
 *  has a working assignment operator
 */
MixedSolventElectrolyte& MixedSolventElectrolyte::
operator=(const MixedSolventElectrolyte& b)
{
    if (&b == this) {
        return *this;
    }

    MolarityIonicVPSSTP::operator=(b);

    numBinaryInteractions_      = b.numBinaryInteractions_ ;
    m_HE_b_ij                   = b.m_HE_b_ij;
    m_HE_c_ij                   = b.m_HE_c_ij;
    m_HE_d_ij                   = b.m_HE_d_ij;
    m_SE_b_ij                   = b.m_SE_b_ij;
    m_SE_c_ij                   = b.m_SE_c_ij;
    m_SE_d_ij                   = b.m_SE_d_ij;
    m_VHE_b_ij                  = b.m_VHE_b_ij;
    m_VHE_c_ij                  = b.m_VHE_c_ij;
    m_VHE_d_ij                  = b.m_VHE_d_ij;
    m_VSE_b_ij                  = b.m_VSE_b_ij;
    m_VSE_c_ij                  = b.m_VSE_c_ij;
    m_VSE_d_ij                  = b.m_VSE_d_ij;
    m_pSpecies_A_ij             = b.m_pSpecies_A_ij;
    m_pSpecies_B_ij             = b.m_pSpecies_B_ij;
    formMargules_               = b.formMargules_;
    formTempModel_              = b.formTempModel_;

    return *this;
}

/**
 *
 * ~MixedSolventElectrolyte():   (virtual)
 *
 * Destructor: does nothing:
 *
 */
MixedSolventElectrolyte::~MixedSolventElectrolyte()
{
}

/*
 * This routine duplicates the current object and returns
 * a pointer to ThermoPhase.
 */
thermo_t*
MixedSolventElectrolyte::duplMyselfAsThermoPhase() const
{
    return new MixedSolventElectrolyte(*this);
}

// Special constructor for a hard-coded problem
/*
 *
 *   LiKCl treating the PseudoBinary layer as passthrough.
 *   -> test to predict the eutectic and liquidus correctly.
 *
 */
MixedSolventElectrolyte::MixedSolventElectrolyte(int testProb)  :
    MolarityIonicVPSSTP(),
    numBinaryInteractions_(0),
    formMargules_(0),
    formTempModel_(0)
{


    initThermoFile("LiKCl_liquid.xml", "");


    numBinaryInteractions_ = 1;

    m_HE_b_ij.resize(1);
    m_HE_c_ij.resize(1);
    m_HE_d_ij.resize(1);

    m_SE_b_ij.resize(1);
    m_SE_c_ij.resize(1);
    m_SE_d_ij.resize(1);

    m_VHE_b_ij.resize(1);
    m_VHE_c_ij.resize(1);
    m_VHE_d_ij.resize(1);

    m_VSE_b_ij.resize(1);
    m_VSE_c_ij.resize(1);
    m_VSE_d_ij.resize(1);

    m_pSpecies_A_ij.resize(1);
    m_pSpecies_B_ij.resize(1);



    m_HE_b_ij[0] = -17570E3;
    m_HE_c_ij[0] = -377.0E3;
    m_HE_d_ij[0] = 0.0;

    m_SE_b_ij[0] = -7.627E3;
    m_SE_c_ij[0] =  4.958E3;
    m_SE_d_ij[0] =  0.0;


    size_t iLiCl = speciesIndex("LiCl(L)");
    if (iLiCl == npos) {
        throw CanteraError("MixedSolventElectrolyte test1 constructor",
                           "Unable to find LiCl(L)");
    }
    m_pSpecies_B_ij[0] = iLiCl;


    size_t iKCl = speciesIndex("KCl(L)");
    if (iKCl == npos) {
        throw CanteraError("MixedSolventElectrolyte test1 constructor",
                           "Unable to find KCl(L)");
    }
    m_pSpecies_A_ij[0] = iKCl;
}


/*
 *  -------------- Utilities -------------------------------
 */


// Equation of state type flag.
/*
 * The ThermoPhase base class returns
 * zero. Subclasses should define this to return a unique
 * non-zero value. Known constants defined for this purpose are
 * listed in mix_defs.h. The MixedSolventElectrolyte class also returns
 * zero, as it is a non-complete class.
 */
int MixedSolventElectrolyte::eosType() const
{
    return 0;
}

//====================================================================================================================
/*
 * ------------ Molar Thermodynamic Properties ----------------------
 */
//====================================================================================================================
/*
 * - Activities, Standard States, Activity Concentrations -----------
 */
//====================================================================================================================
// Get the array of non-dimensional molar-based activity coefficients at
// the current solution temperature, pressure, and solution concentration.
/*
 * @param ac Output vector of activity coefficients. Length: m_kk.
 */
void MixedSolventElectrolyte::getActivityCoefficients(doublereal* ac) const
{
    /*
     * Update the activity coefficients
     */
    s_update_lnActCoeff();

    /*
     * take the exp of the internally stored coefficients.
     */
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = exp(lnActCoeff_Scaled_[k]);
    }
}
//====================================================================================================================
/*
 * ------------ Partial Molar Properties of the Solution ------------
 */



void MixedSolventElectrolyte::getElectrochemPotentials(doublereal* mu) const
{
    getChemPotentials(mu);
    double ve = Faraday * electricPotential();
    for (size_t k = 0; k < m_kk; k++) {
        mu[k] += ve*charge(k);
    }
}


void MixedSolventElectrolyte::getChemPotentials(doublereal* mu) const
{
    doublereal xx;
    /*
     * First get the standard chemical potentials in
     * molar form.
     *  -> this requires updates of standard state as a function
     *     of T and P
     */
    getStandardChemPotentials(mu);
    /*
     * Update the activity coefficients
     */
    s_update_lnActCoeff();
    /*
     *
     */
    doublereal RT = GasConstant * temperature();
    for (size_t k = 0; k < m_kk; k++) {
        xx = std::max(moleFractions_[k], SmallNumber);
        mu[k] += RT * (log(xx) + lnActCoeff_Scaled_[k]);
    }
}

/// Molar enthalpy. Units: J/kmol.
doublereal MixedSolventElectrolyte::enthalpy_mole() const
{
    size_t kk = nSpecies();
    double h = 0;
    vector_fp hbar(kk);
    getPartialMolarEnthalpies(&hbar[0]);
    for (size_t i = 0; i < kk; i++) {
        h += moleFractions_[i]*hbar[i];
    }
    return h;
}

/// Molar entropy. Units: J/kmol.
doublereal MixedSolventElectrolyte::entropy_mole() const
{
    size_t kk = nSpecies();
    double s = 0;
    vector_fp sbar(kk);
    getPartialMolarEntropies(&sbar[0]);
    for (size_t i = 0; i < kk; i++) {
        s += moleFractions_[i]*sbar[i];
    }
    return s;
}

/// Molar heat capacity at constant pressure. Units: J/kmol/K.
doublereal MixedSolventElectrolyte::cp_mole() const
{
    size_t kk = nSpecies();
    double cp = 0;
    vector_fp cpbar(kk);
    getPartialMolarCp(&cpbar[0]);
    for (size_t i = 0; i < kk; i++) {
        cp += moleFractions_[i]*cpbar[i];
    }
    return cp;
}

/// Molar heat capacity at constant volume. Units: J/kmol/K.
doublereal MixedSolventElectrolyte::cv_mole() const
{
    return cp_mole() - GasConstant;
}

// Returns an array of partial molar enthalpies for the species
// in the mixture.
/*
 * Units (J/kmol)
 *
 * For this phase, the partial molar enthalpies are equal to the
 * standard state enthalpies modified by the derivative of the
 * molality-based activity coefficient wrt temperature
 *
 *  \f[
 * \bar h_k(T,P) = h^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
 * \f]
 *
 */
void MixedSolventElectrolyte::getPartialMolarEnthalpies(doublereal* hbar) const
{
    /*
     * Get the nondimensional standard state enthalpies
     */
    getEnthalpy_RT(hbar);
    /*
     * dimensionalize it.
     */
    double T = temperature();
    double RT = GasConstant * T;
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] *= RT;
    }
    /*
     * Update the activity coefficients, This also update the
     * internally stored molalities.
     */
    s_update_lnActCoeff();
    s_update_dlnActCoeff_dT();
    double RTT = RT * T;
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] -= RTT * dlnActCoeffdT_Scaled_[k];
    }
}

// Returns an array of partial molar heat capacities for the species
// in the mixture.
/*
 * Units (J/kmol)
 *
 * For this phase, the partial molar enthalpies are equal to the
 * standard state enthalpies modified by the derivative of the
 * activity coefficient wrt temperature
 *
 *  \f[
 * ??????????? \bar s_k(T,P) = s^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
 * \f]
 *
 */
void MixedSolventElectrolyte::getPartialMolarCp(doublereal* cpbar) const
{
    /*
     * Get the nondimensional standard state entropies
     */
    getCp_R(cpbar);
    double T = temperature();
    /*
     * Update the activity coefficients, This also update the
     * internally stored molalities.
     */
    s_update_lnActCoeff();
    s_update_dlnActCoeff_dT();

    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] -= 2 * T * dlnActCoeffdT_Scaled_[k] + T * T * d2lnActCoeffdT2_Scaled_[k];
    }
    /*
     * dimensionalize it.
     */
    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] *= GasConstant;
    }
}

// Returns an array of partial molar entropies for the species
// in the mixture.
/*
 * Units (J/kmol)
 *
 * For this phase, the partial molar enthalpies are equal to the
 * standard state enthalpies modified by the derivative of the
 * activity coefficient wrt temperature
 *
 *  \f[
 * \bar s_k(T,P) = s^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
 * \f]
 *
 */
void MixedSolventElectrolyte::getPartialMolarEntropies(doublereal* sbar) const
{
    double xx;
    /*
     * Get the nondimensional standard state entropies
     */
    getEntropy_R(sbar);
    double T = temperature();
    /*
     * Update the activity coefficients, This also update the
     * internally stored molalities.
     */
    s_update_lnActCoeff();
    s_update_dlnActCoeff_dT();

    for (size_t k = 0; k < m_kk; k++) {
        xx = std::max(moleFractions_[k], SmallNumber);
        sbar[k] += - lnActCoeff_Scaled_[k] -log(xx) - T * dlnActCoeffdT_Scaled_[k];
    }
    /*
     * dimensionalize it.
     */
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] *= GasConstant;
    }
}

/*
 * ------------ Partial Molar Properties of the Solution ------------
 */

// Return an array of partial molar volumes for the
// species in the mixture. Units: m^3/kmol.
/*
 *  Frequently, for this class of thermodynamics representations,
 *  the excess Volume due to mixing is zero. Here, we set it as
 *  a default. It may be overridden in derived classes.
 *
 *  @param vbar   Output vector of species partial molar volumes.
 *                Length = m_kk. units are m^3/kmol.
 */
void MixedSolventElectrolyte::getPartialMolarVolumes(doublereal* vbar) const
{
    int delAK, delBK;
    double XA, XB, g0 , g1;
    double T = temperature();

    /*
     * Get the standard state values in m^3 kmol-1
     */
    getStandardVolumes(vbar);

    for (size_t iK = 0; iK < m_kk; iK++) {
        delAK = 0;
        delBK = 0;
        for (size_t i = 0; i <  numBinaryInteractions_; i++) {
            size_t iA =  m_pSpecies_A_ij[i];
            size_t iB =  m_pSpecies_B_ij[i];

            if (iA==iK) {
                delAK = 1;
            } else if (iB==iK) {
                delBK = 1;
            }

            XA = moleFractions_[iA];
            XB = moleFractions_[iB];

            g0 = (m_VHE_b_ij[i] - T * m_VSE_b_ij[i]);
            g1 = (m_VHE_c_ij[i] - T * m_VSE_c_ij[i]);

            vbar[iK] += XA*XB*(g0+g1*XB)+((delAK-XA)*XB+XA*(delBK-XB))*(g0+g1*XB)+XA*XB*(delBK-XB)*g1;
        }
    }
}

doublereal MixedSolventElectrolyte::err(const std::string& msg) const
{
    throw CanteraError("MixedSolventElectrolyte","Base class method "
                       +msg+" called. Equation of state type: "+int2str(eosType()));
    return 0;
}


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
void MixedSolventElectrolyte::initThermo()
{
    initLengths();
    MolarityIonicVPSSTP::initThermo();
}


//   Initialize lengths of local variables after all species have
//   been identified.
void  MixedSolventElectrolyte::initLengths()
{
    m_kk = nSpecies();
    dlnActCoeffdlnN_.resize(m_kk, m_kk);
}

/*
 * initThermoXML()                (virtual from ThermoPhase)
 *   Import and initialize a ThermoPhase object
 *
 * @param phaseNode This object must be the phase node of a
 *             complete XML tree
 *             description of the phase, including all of the
 *             species data. In other words while "phase" must
 *             point to an XML phase object, it must have
 *             sibling nodes "speciesData" that describe
 *             the species in the phase.
 * @param id   ID of the phase. If nonnull, a check is done
 *             to see if phaseNode is pointing to the phase
 *             with the correct id.
 */
void MixedSolventElectrolyte::initThermoXML(XML_Node& phaseNode, const std::string& id)
{
    string subname = "MixedSolventElectrolyte::initThermoXML";
    string stemp;

    if ((int) id.size() > 0) {
        string idp = phaseNode.id();
        if (idp != id) {
            throw CanteraError(subname, "phasenode and Id are incompatible");
        }
    }

    /*
     * Check on the thermo field. Must have:
     * <thermo model="MixedSolventElectrolyte" />
     */
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError(subname, "no thermo XML node");
    }
    XML_Node& thermoNode = phaseNode.child("thermo");
    string mStringa = thermoNode.attrib("model");
    string mString = lowercase(mStringa);
    if (mString != "MixedSolventElectrolyte") {
        throw CanteraError(subname, "Unknown thermo model: " + mStringa);
    }

    /*
     * Go get all of the coefficients and factors in the
     * activityCoefficients XML block
     */
    XML_Node* acNodePtr = 0;
    if (thermoNode.hasChild("activityCoefficients")) {
        XML_Node& acNode = thermoNode.child("activityCoefficients");
        acNodePtr = &acNode;
        string mStringa = acNode.attrib("model");
        string mString = lowercase(mStringa);
        if (mString != "margules") {
            throw CanteraError(subname.c_str(),
                               "Unknown activity coefficient model: " + mStringa);
        }
        size_t n = acNodePtr->nChildren();
        for (size_t i = 0; i < n; i++) {
            XML_Node& xmlACChild = acNodePtr->child(i);
            stemp = xmlACChild.name();
            string nodeName = lowercase(stemp);
            /*
             * Process a binary salt field, or any of the other XML fields
             * that make up the Pitzer Database. Entries will be ignored
             * if any of the species in the entry isn't in the solution.
             */
            if (nodeName == "binaryneutralspeciesparameters") {
                readXMLBinarySpecies(xmlACChild);

            }
        }
    }

    /*
     * Go down the chain
     */
    MolarityIonicVPSSTP::initThermoXML(phaseNode, id);


}
//===================================================================================================================

// Update the activity coefficients
/*
 * This function will be called to update the internally stored
 * natural logarithm of the activity coefficients
 *
 *   he = X_A X_B(B + C X_B)
 */
void MixedSolventElectrolyte::s_update_lnActCoeff() const
{
    int delAK, delBK;
    double XA, XB, g0, g1;
    double T = temperature();
    double RT = GasConstant*T;
    lnActCoeff_Scaled_.assign(m_kk, 0.0);
    for (size_t iK = 0; iK < m_kk; iK++) {
        for (size_t i = 0; i <  numBinaryInteractions_; i++) {
            size_t iA =  m_pSpecies_A_ij[i];
            size_t iB =  m_pSpecies_B_ij[i];
            delAK = 0;
            delBK = 0;
            if (iA==iK) {
                delAK = 1;
            } else if (iB==iK) {
                delBK = 1;
            }
            XA = moleFractions_[iA];
            XB = moleFractions_[iB];
            g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / RT;
            g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / RT;
            lnActCoeff_Scaled_[iK] += (delAK * XB + XA * delBK - XA * XB) * (g0 + g1 * XB) + XA * XB * (delBK - XB) * g1;
        }
    }
}
//===================================================================================================================
// Update the derivative of the log of the activity coefficients wrt T
/*
 * This function will be called to update the internally stored
 * natural logarithm of the activity coefficients
 *
 *   he = X_A X_B(B + C X_B)
 */
void MixedSolventElectrolyte::s_update_dlnActCoeff_dT() const
{
    int delAK, delBK;
    doublereal XA, XB, g0, g1;
    doublereal T = temperature();
    doublereal RTT = GasConstant*T*T;
    dlnActCoeffdT_Scaled_.assign(m_kk, 0.0);
    d2lnActCoeffdT2_Scaled_.assign(m_kk, 0.0);
    for (size_t iK = 0; iK < m_kk; iK++) {
        for (size_t i = 0; i <  numBinaryInteractions_; i++) {
            size_t iA =  m_pSpecies_A_ij[i];
            size_t iB =  m_pSpecies_B_ij[i];
            delAK = 0;
            delBK = 0;
            if (iA==iK) {
                delAK = 1;
            } else if (iB==iK) {
                delBK = 1;
            }
            XA = moleFractions_[iA];
            XB = moleFractions_[iB];
            g0 = -m_HE_b_ij[i] / RTT;
            g1 = -m_HE_c_ij[i] / RTT;
            double temp = (delAK * XB + XA * delBK - XA * XB) * (g0 + g1 * XB) + XA * XB * (delBK - XB) * g1;
            dlnActCoeffdT_Scaled_[iK] += temp;
            d2lnActCoeffdT2_Scaled_[iK] -= 2.0 * temp / T;
        }
    }
}
//====================================================================================================================
void MixedSolventElectrolyte::getdlnActCoeffdT(doublereal* dlnActCoeffdT) const
{
    s_update_dlnActCoeff_dT();
    for (size_t k = 0; k < m_kk; k++) {
        dlnActCoeffdT[k] = dlnActCoeffdT_Scaled_[k];
    }
}
//====================================================================================================================
void MixedSolventElectrolyte::getd2lnActCoeffdT2(doublereal* d2lnActCoeffdT2) const
{
    s_update_dlnActCoeff_dT();
    for (size_t k = 0; k < m_kk; k++) {
        d2lnActCoeffdT2[k] = d2lnActCoeffdT2_Scaled_[k];
    }
}
//====================================================================================================================

// Get the change in activity coefficients w.r.t. change in state (temp, mole fraction, etc.) along
// a line in parameter space or along a line in physical space
/*
 *
 * @param dTds           Input of temperature change along the path
 * @param dXds           Input vector of changes in mole fraction along the path. length = m_kk
 *                       Along the path length it must be the case that the mole fractions sum to one.
 * @param dlnActCoeffds  Output vector of the directional derivatives of the
 *                       log Activity Coefficients along the path. length = m_kk
 *  units are 1/units(s). if s is a physical coordinate then the units are 1/m.
 */
void  MixedSolventElectrolyte::getdlnActCoeffds(const doublereal dTds, const doublereal* const dXds,
        doublereal* dlnActCoeffds) const
{
    int delAK, delBK;
    double XA, XB, g0, g1, dXA, dXB;
    double T = temperature();
    double RT = GasConstant*T;

    //fvo_zero_dbl_1(dlnActCoeff, m_kk);
    s_update_dlnActCoeff_dT();

    for (size_t iK = 0; iK < m_kk; iK++) {
        dlnActCoeffds[iK] = 0.0;
        for (size_t i = 0; i <  numBinaryInteractions_; i++) {

            size_t iA =  m_pSpecies_A_ij[i];
            size_t iB =  m_pSpecies_B_ij[i];

            delAK = 0;
            delBK = 0;

            if (iA==iK) {
                delAK = 1;
            } else if (iB==iK) {
                delBK = 1;
            }

            XA = moleFractions_[iA];
            XB = moleFractions_[iB];

            dXA = dXds[iA];
            dXB = dXds[iB];

            g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / RT;
            g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / RT;

            dlnActCoeffds[iK] += ((delBK-XB)*dXA + (delAK-XA)*dXB)*(g0+2*g1*XB) + (delBK-XB)*2*g1*XA*dXB
                                 + dlnActCoeffdT_Scaled_[iK]*dTds;
        }
    }
}
//====================================================================================================================
// Update the derivative of the log of the activity coefficients wrt dlnN
/*
 * This function will be called to update the internally stored gradients of the
 * logarithm of the activity coefficients.  These are used in the determination
 * of the diffusion coefficients.
 *
 *   he = X_A X_B(B + C X_B)
 */
void MixedSolventElectrolyte::s_update_dlnActCoeff_dlnN_diag() const
{
    int delAK, delBK;
    double XA, XB, XK, g0, g1;
    double T = temperature();
    double RT = GasConstant*T;

    dlnActCoeffdlnN_diag_.assign(m_kk, 0);

    for (size_t iK = 0; iK < m_kk; iK++) {

        XK = moleFractions_[iK];

        for (size_t i = 0; i <  numBinaryInteractions_; i++) {

            size_t iA =  m_pSpecies_A_ij[i];
            size_t iB =  m_pSpecies_B_ij[i];

            delAK = 0;
            delBK = 0;

            if (iA==iK) {
                delAK = 1;
            } else if (iB==iK) {
                delBK = 1;
            }

            XA = moleFractions_[iA];
            XB = moleFractions_[iB];

            g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / RT;
            g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / RT;

            dlnActCoeffdlnN_diag_[iK] += 2*(delBK-XB)*(g0*(delAK-XA)+g1*(2*(delAK-XA)*XB+XA*(delBK-XB)));

            // double gfac = g0 + g1 * XB;
            // double gggg = (delBK - XB) * g1;


            // dlnActCoeffdlnN_diag_[iK] += gfac * delAK * ( - XB + delBK);

            // dlnActCoeffdlnN_diag_[iK] += gfac * delBK * ( - XA + delAK);

            // dlnActCoeffdlnN_diag_[iK] += gfac * (2.0 * XA * XB - delAK * XB - XA * delBK);

            // dlnActCoeffdlnN_diag_[iK] += (delAK * XB + XA * delBK - XA * XB) * g1 * (-XB + delBK);

            // dlnActCoeffdlnN_diag_[iK] += gggg * ( - 2.0 * XA * XB + delAK * XB + XA * delBK);

            // dlnActCoeffdlnN_diag_[iK] += - g1 * XA * XB * (- XB + delBK);
        }
        dlnActCoeffdlnN_diag_[iK] = XK*dlnActCoeffdlnN_diag_[iK];//-XK;
    }
}

//====================================================================================================================
// Update the derivative of the log of the activity coefficients wrt dlnN
/*
 * This function will be called to update the internally stored gradients of the
 * logarithm of the activity coefficients.  These are used in the determination
 * of the diffusion coefficients.
 *
 */
void MixedSolventElectrolyte::s_update_dlnActCoeff_dlnN() const
{
    doublereal delAK, delBK;
    double XA, XB, g0, g1,XM;
    double T = temperature();
    double RT = GasConstant*T;

    doublereal delAM, delBM;
    dlnActCoeffdlnN_.zero();

    /*
     *  Loop over the activity coefficient gamma_k
     */
    for (size_t iK = 0; iK < m_kk; iK++) {
        for (size_t iM = 0; iM < m_kk; iM++) {
            XM = moleFractions_[iM];
            for (size_t i = 0; i <  numBinaryInteractions_; i++) {

                size_t iA =  m_pSpecies_A_ij[i];
                size_t iB =  m_pSpecies_B_ij[i];

                delAK = 0.0;
                delBK = 0.0;
                delAM = 0.0;
                delBM = 0.0;
                if (iA==iK) {
                    delAK = 1.0;
                } else if (iB==iK) {
                    delBK = 1.0;
                }
                if (iA==iM) {
                    delAM = 1.0;
                } else if (iB==iM) {
                    delBM = 1.0;
                }

                XA = moleFractions_[iA];
                XB = moleFractions_[iB];

                g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / RT;
                g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / RT;

                dlnActCoeffdlnN_(iK,iM) += g0*((delAM-XA)*(delBK-XB)+(delAK-XA)*(delBM-XB));
                dlnActCoeffdlnN_(iK,iM) += 2*g1*((delAM-XA)*(delBK-XB)*XB+(delAK-XA)*(delBM-XB)*XB+(delBM-XB)*(delBK-XB)*XA);

                // double gfac = g0 + g1 * XB;
                // double gggg = (delBK - XB) * g1;


                // dlnActCoeffdlnN_(iK, iM) += gfac * delAK * ( - XB + delBM);

                // dlnActCoeffdlnN_(iK, iM) += gfac * delBK * ( - XA + delAM);

                // dlnActCoeffdlnN_(iK, iM) += gfac * (2.0 * XA * XB - delAM * XB - XA * delBM);

                // dlnActCoeffdlnN_(iK, iM) += (delAK * XB + XA * delBK - XA * XB) * g1 * (-XB + delBM);

                // dlnActCoeffdlnN_(iK, iM) += gggg * ( - 2.0 * XA * XB + delAM * XB + XA * delBM);

                // dlnActCoeffdlnN_(iK, iM) += - g1 * XA * XB * (- XB + delBM);
            }
            dlnActCoeffdlnN_(iK,iM) = XM*dlnActCoeffdlnN_(iK,iM);
        }
    }
}
//====================================================================================================================
void MixedSolventElectrolyte::s_update_dlnActCoeff_dlnX_diag() const
{
    doublereal XA, XB, g0 , g1;
    doublereal T = temperature();

    dlnActCoeffdlnX_diag_.assign(m_kk, 0);
    doublereal RT = GasConstant * T;

    for (size_t i = 0; i <  numBinaryInteractions_; i++) {

        size_t iA =  m_pSpecies_A_ij[i];
        size_t iB =  m_pSpecies_B_ij[i];

        XA = moleFractions_[iA];
        XB = moleFractions_[iB];

        g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / RT;
        g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / RT;

        dlnActCoeffdlnX_diag_[iA] += XA*XB*(2*g1*-2*g0-6*g1*XB);
        dlnActCoeffdlnX_diag_[iB] += XA*XB*(2*g1*-2*g0-6*g1*XB);
    }
}

//====================================================================================================================
void MixedSolventElectrolyte::getdlnActCoeffdlnN_diag(doublereal* dlnActCoeffdlnN_diag) const
{
    s_update_dlnActCoeff_dlnN_diag();
    for (size_t k = 0; k < m_kk; k++) {
        dlnActCoeffdlnN_diag[k] = dlnActCoeffdlnN_diag_[k];
    }
}
//====================================================================================================================
void MixedSolventElectrolyte::getdlnActCoeffdlnX_diag(doublereal* dlnActCoeffdlnX_diag) const
{
    s_update_dlnActCoeff_dlnX_diag();
    for (size_t k = 0; k < m_kk; k++) {
        dlnActCoeffdlnX_diag[k] = dlnActCoeffdlnX_diag_[k];
    }
}
//====================================================================================================================
void MixedSolventElectrolyte::getdlnActCoeffdlnN(const size_t ld, doublereal* dlnActCoeffdlnN)
{
    s_update_dlnActCoeff_dlnN();
    double* data =  & dlnActCoeffdlnN_(0,0);
    for (size_t k = 0; k < m_kk; k++) {
        for (size_t m = 0; m < m_kk; m++) {
            dlnActCoeffdlnN[ld * k + m] = data[m_kk * k + m];
        }
    }
}
//====================================================================================================================
void MixedSolventElectrolyte::resizeNumInteractions(const size_t num)
{
    numBinaryInteractions_ = num;
    m_HE_b_ij.resize(num, 0.0);
    m_HE_c_ij.resize(num, 0.0);
    m_HE_d_ij.resize(num, 0.0);
    m_SE_b_ij.resize(num, 0.0);
    m_SE_c_ij.resize(num, 0.0);
    m_SE_d_ij.resize(num, 0.0);
    m_VHE_b_ij.resize(num, 0.0);
    m_VHE_c_ij.resize(num, 0.0);
    m_VHE_d_ij.resize(num, 0.0);
    m_VSE_b_ij.resize(num, 0.0);
    m_VSE_c_ij.resize(num, 0.0);
    m_VSE_d_ij.resize(num, 0.0);

    m_pSpecies_A_ij.resize(num, npos);
    m_pSpecies_B_ij.resize(num, npos);

}
//====================================================================================================================

/*
 * Process an XML node called "binaryNeutralSpeciesParameters"
 * This node contains all of the parameters necessary to describe
 * the Margules Interaction for a single binary interaction
 * This function reads the XML file and writes the coefficients
 * it finds to an internal data structures.
 */
void MixedSolventElectrolyte::readXMLBinarySpecies(XML_Node& xmLBinarySpecies)
{
    string xname = xmLBinarySpecies.name();
    if (xname != "binaryNeutralSpeciesParameters") {
        throw CanteraError("MixedSolventElectrolyte::readXMLBinarySpecies",
                           "Incorrect name for processing this routine: " + xname);
    }
    double* charge = DATA_PTR(m_speciesCharge);
    string stemp;
    size_t nParamsFound;
    vector_fp vParams;
    string iName = xmLBinarySpecies.attrib("speciesA");
    if (iName == "") {
        throw CanteraError("MixedSolventElectrolyte::readXMLBinarySpecies", "no speciesA attrib");
    }
    string jName = xmLBinarySpecies.attrib("speciesB");
    if (jName == "") {
        throw CanteraError("MixedSolventElectrolyte::readXMLBinarySpecies", "no speciesB attrib");
    }
    /*
     * Find the index of the species in the current phase. It's not
     * an error to not find the species
     */
    size_t iSpecies = speciesIndex(iName);
    if (iSpecies == npos) {
        return;
    }
    string ispName = speciesName(iSpecies);
    if (charge[iSpecies] != 0) {
        throw CanteraError("MixedSolventElectrolyte::readXMLBinarySpecies", "speciesA charge problem");
    }
    size_t jSpecies = speciesIndex(jName);
    if (jSpecies == npos) {
        return;
    }
    string jspName = speciesName(jSpecies);
    if (charge[jSpecies] != 0) {
        throw CanteraError("MixedSolventElectrolyte::readXMLBinarySpecies", "speciesB charge problem");
    }

    resizeNumInteractions(numBinaryInteractions_ + 1);
    size_t iSpot = numBinaryInteractions_ - 1;
    m_pSpecies_A_ij[iSpot] = iSpecies;
    m_pSpecies_B_ij[iSpot] = jSpecies;

    size_t num = xmLBinarySpecies.nChildren();
    for (size_t iChild = 0; iChild < num; iChild++) {
        XML_Node& xmlChild = xmLBinarySpecies.child(iChild);
        stemp = xmlChild.name();
        string nodeName = lowercase(stemp);
        /*
         * Process the binary species interaction child elements
         */
        if (nodeName == "excessenthalpy") {
            /*
             * Get the string containing all of the values
             */
            ctml::getFloatArray(xmlChild, vParams, true, "toSI", "excessEnthalpy");
            nParamsFound = vParams.size();

            if (nParamsFound != 2) {
                throw CanteraError("MixedSolventElectrolyte::readXMLBinarySpecies::excessEnthalpy for " + ispName
                                   + "::" + jspName,
                                   "wrong number of params found");
            }
            m_HE_b_ij[iSpot] = vParams[0];
            m_HE_c_ij[iSpot] = vParams[1];
        }

        if (nodeName == "excessentropy") {
            /*
             * Get the string containing all of the values
             */
            ctml::getFloatArray(xmlChild, vParams, true, "toSI", "excessEntropy");
            nParamsFound = vParams.size();

            if (nParamsFound != 2) {
                throw CanteraError("MixedSolventElectrolyte::readXMLBinarySpecies::excessEntropy for " + ispName
                                   + "::" + jspName,
                                   "wrong number of params found");
            }
            m_SE_b_ij[iSpot] = vParams[0];
            m_SE_c_ij[iSpot] = vParams[1];
        }

        if (nodeName == "excessvolume_enthalpy") {
            /*
             * Get the string containing all of the values
             */
            ctml::getFloatArray(xmlChild, vParams, true, "toSI", "excessVolume_Enthalpy");
            nParamsFound = vParams.size();

            if (nParamsFound != 2) {
                throw CanteraError("MixedSolventElectrolyte::readXMLBinarySpecies::excessVolume_Enthalpy for " + ispName
                                   + "::" + jspName,
                                   "wrong number of params found");
            }
            m_VHE_b_ij[iSpot] = vParams[0];
            m_VHE_c_ij[iSpot] = vParams[1];
        }

        if (nodeName == "excessvolume_entropy") {
            /*
             * Get the string containing all of the values
             */
            ctml::getFloatArray(xmlChild, vParams, true, "toSI", "excessVolume_Entropy");
            nParamsFound = vParams.size();

            if (nParamsFound != 2) {
                throw CanteraError("MixedSolventElectrolyte::readXMLBinarySpecies::excessVolume_Entropy for " + ispName
                                   + "::" + jspName,
                                   "wrong number of params found");
            }
            m_VSE_b_ij[iSpot] = vParams[0];
            m_VSE_c_ij[iSpot] = vParams[1];
        }


    }

}

}

