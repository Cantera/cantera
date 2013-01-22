/**
 *  @file RedlichKisterVPSSTP.cpp
 *   Definitions for ThermoPhase object for phases which
 *   employ excess gibbs free energy formulations related to RedlichKister
 *   expansions (see \ref thermoprops
 *    and class \link Cantera::RedlichKisterVPSSTP RedlichKisterVPSSTP\endlink).
 *
 */
/*
 * Copyright (2009) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/thermo/RedlichKisterVPSSTP.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"

#include "cantera/base/ct_defs.h"

#include <iomanip>
#include <fstream>

using namespace std;

namespace Cantera
{

//====================================================================================================================
/*
 * Default constructor.
 *
 */
RedlichKisterVPSSTP::RedlichKisterVPSSTP() :
    GibbsExcessVPSSTP(),
    numBinaryInteractions_(0),
    m_pSpecies_A_ij(0),
    m_pSpecies_B_ij(0),
    m_N_ij(0),
    m_HE_m_ij(0),
    m_SE_m_ij(0),
    formRedlichKister_(0),
    formTempModel_(0),
    dlnActCoeff_dX_()
{
}
//====================================================================================================================
/*
 * Working constructors
 *
 *  The two constructors below are the normal way
 *  the phase initializes itself. They are shells that call
 *  the routine initThermo(), with a reference to the
 *  XML database to get the info for the phase.

 */
RedlichKisterVPSSTP::RedlichKisterVPSSTP(const std::string& inputFile,
                                         const std::string& id) :
    GibbsExcessVPSSTP(),
    numBinaryInteractions_(0),
    m_pSpecies_A_ij(0),
    m_pSpecies_B_ij(0),
    m_N_ij(0),
    m_HE_m_ij(0),
    m_SE_m_ij(0),
    formRedlichKister_(0),
    formTempModel_(0),
    dlnActCoeff_dX_()
{
    initThermoFile(inputFile, id);
}
//====================================================================================================================
RedlichKisterVPSSTP::RedlichKisterVPSSTP(XML_Node& phaseRoot,
                                         const std::string& id) :
    GibbsExcessVPSSTP(),
    numBinaryInteractions_(0),
    m_pSpecies_A_ij(0),
    m_pSpecies_B_ij(0),
    m_N_ij(0),
    m_HE_m_ij(0),
    m_SE_m_ij(0),
    formRedlichKister_(0),
    formTempModel_(0),
    dlnActCoeff_dX_()
{
    importPhase(*findXMLPhase(&phaseRoot, id), this);
}
//====================================================================================================================
// Special constructor for a hard-coded problem
/*
 *
 *   LiKCl treating the PseudoBinary layer as passthrough.
 *   -> test to predict the eutectic and liquidus correctly.
 *
 */
RedlichKisterVPSSTP::RedlichKisterVPSSTP(int testProb)  :
    GibbsExcessVPSSTP(),
    numBinaryInteractions_(0),
    m_pSpecies_A_ij(0),
    m_pSpecies_B_ij(0),
    m_N_ij(0),
    m_HE_m_ij(0),
    m_SE_m_ij(0),
    formRedlichKister_(0),
    formTempModel_(0),
    dlnActCoeff_dX_()
{
    initThermoFile("LiKCl_liquid.xml", "");
    numBinaryInteractions_ = 1;

    m_HE_m_ij.resize(0);
    m_SE_m_ij.resize(0);

    vector_fp he(2);
    he[0] = 0.0;
    he[1] = 0.0;
    vector_fp se(2);
    se[0] = 0.0;
    se[1] = 0.0;

    m_HE_m_ij.push_back(he);
    m_SE_m_ij.push_back(se);
    m_N_ij.push_back(1);
    m_pSpecies_A_ij.resize(1);
    m_pSpecies_B_ij.resize(1);

    size_t iLiLi = speciesIndex("LiLi");
    if (iLiLi == npos) {
        throw CanteraError("RedlichKisterVPSSTP test1 constructor",
                           "Unable to find LiLi");
    }
    m_pSpecies_A_ij[0] = iLiLi;


    size_t iVLi = speciesIndex("VLi");
    if (iVLi == npos) {
        throw CanteraError("RedlichKisterVPSSTP test1 constructor",
                           "Unable to find VLi");
    }
    m_pSpecies_B_ij[0] = iVLi;


}
//====================================================================================================================
/*
 * Copy Constructor:
 *
 *  Note this stuff will not work until the underlying phase
 *  has a working copy constructor
 */
RedlichKisterVPSSTP::RedlichKisterVPSSTP(const RedlichKisterVPSSTP& b) :
    GibbsExcessVPSSTP(),
    numBinaryInteractions_(0),
    m_pSpecies_A_ij(0),
    m_pSpecies_B_ij(0),
    m_N_ij(0),
    m_HE_m_ij(0),
    m_SE_m_ij(0),
    formRedlichKister_(0),
    formTempModel_(0),
    dlnActCoeff_dX_()
{
    RedlichKisterVPSSTP::operator=(b);
}
//====================================================================================================================
/*
 * operator=()
 *
 *  Note this stuff will not work until the underlying phase
 *  has a working assignment operator
 */
RedlichKisterVPSSTP& RedlichKisterVPSSTP::
operator=(const RedlichKisterVPSSTP& b)
{
    if (&b == this) {
        return *this;
    }

    GibbsExcessVPSSTP::operator=(b);

    numBinaryInteractions_      = b.numBinaryInteractions_ ;
    m_pSpecies_A_ij             = b.m_pSpecies_A_ij;
    m_pSpecies_B_ij             = b.m_pSpecies_B_ij;
    m_N_ij                      = b.m_N_ij;
    m_HE_m_ij                   = b.m_HE_m_ij;
    m_SE_m_ij                   = b.m_SE_m_ij;
    formRedlichKister_          = b.formRedlichKister_;
    formTempModel_              = b.formTempModel_;
    dlnActCoeff_dX_             = b.dlnActCoeff_dX_;

    return *this;
}
//====================================================================================================================
/*
 *
 * ~RedlichKisterVPSSTP():   (virtual)
 *
 * Destructor: does nothing:
 *
 */
RedlichKisterVPSSTP::~RedlichKisterVPSSTP()
{
}
//====================================================================================================================
/*
 * This routine duplicates the current object and returns
 * a pointer to ThermoPhase.
 */
thermo_t*
RedlichKisterVPSSTP::duplMyselfAsThermoPhase() const
{
    return new RedlichKisterVPSSTP(*this);
}

//====================================================================================================================
// Equation of state type flag.
/*
 * The ThermoPhase base class returns
 * zero. Subclasses should define this to return a unique
 * non-zero value. Known constants defined for this purpose are
 * listed in mix_defs.h. The RedlichKisterVPSSTP class also returns
 * zero, as it is a non-complete class.
 */
int RedlichKisterVPSSTP::eosType() const
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

void RedlichKisterVPSSTP::getLnActivityCoefficients(doublereal* lnac) const
{
    /*
     * Update the activity coefficients
     */
    s_update_lnActCoeff();

    /*
     * take the exp of the internally stored coefficients.
     */
    for (size_t k = 0; k < m_kk; k++) {
        lnac[k] = lnActCoeff_Scaled_[k];
    }
}
//====================================================================================================================
/*
 * ------------ Partial Molar Properties of the Solution ------------
 */
//====================================================================================================================
void RedlichKisterVPSSTP::getElectrochemPotentials(doublereal* mu) const
{
    getChemPotentials(mu);
    double ve = Faraday * electricPotential();
    for (size_t k = 0; k < m_kk; k++) {
        mu[k] += ve*charge(k);
    }
}
//====================================================================================================================
void RedlichKisterVPSSTP::getChemPotentials(doublereal* mu) const
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
//====================================================================================================================
//Molar enthalpy. Units: J/kmol.
doublereal RedlichKisterVPSSTP::enthalpy_mole() const
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
//====================================================================================================================
/// Molar entropy. Units: J/kmol.
doublereal RedlichKisterVPSSTP::entropy_mole() const
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
//====================================================================================================================
/// Molar heat capacity at constant pressure. Units: J/kmol/K.
doublereal RedlichKisterVPSSTP::cp_mole() const
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
//====================================================================================================================
/// Molar heat capacity at constant volume. Units: J/kmol/K.
doublereal RedlichKisterVPSSTP::cv_mole() const
{
    return cp_mole() - GasConstant;
}
//====================================================================================================================
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
void RedlichKisterVPSSTP::getPartialMolarEnthalpies(doublereal* hbar) const
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
//====================================================================================================================
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
void RedlichKisterVPSSTP::getPartialMolarCp(doublereal* cpbar) const
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
//====================================================================================================================
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
void RedlichKisterVPSSTP::getPartialMolarEntropies(doublereal* sbar) const
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
//====================================================================================================================
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
void RedlichKisterVPSSTP::getPartialMolarVolumes(doublereal* vbar) const
{
    /*
     * Get the standard state values in m^3 kmol-1
     */
    getStandardVolumes(vbar);
    for (size_t iK = 0; iK < m_kk; iK++) {

        vbar[iK] += 0.0;
    }
}
//====================================================================================================================
doublereal RedlichKisterVPSSTP::err(const std::string& msg) const
{
    throw CanteraError("RedlichKisterVPSSTP","Base class method "
                       +msg+" called. Equation of state type: "+int2str(eosType()));
    return 0;
}
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
void RedlichKisterVPSSTP::initThermo()
{
    initLengths();
    GibbsExcessVPSSTP::initThermo();
}
//====================================================================================================================
//   Initialize lengths of local variables after all species have
//   been identified.
void  RedlichKisterVPSSTP::initLengths()
{
    m_kk = nSpecies();
    dlnActCoeffdlnN_.resize(m_kk, m_kk);
}
//====================================================================================================================
/*
 * initThermoXML()                (virtual from ThermoPhase)
 *   Import and initialize a ThermoPhase object
 *
 * @param phaseNode This object must be the phase node of a complete XML tree
 *                  description of the phase, including all of the species data. In other words while "phase" must
 *                  point to an XML phase object, it must have sibling nodes "speciesData" that describe
 *             the species in the phase.
 * @param id   ID of the phase. If nonnull, a check is done to see if phaseNode is pointing to the phase
 *             with the correct id.
 */
void RedlichKisterVPSSTP::initThermoXML(XML_Node& phaseNode, const std::string& id)
{
    std::string subname = "RedlichKisterVPSSTP::initThermoXML";
    std::string stemp;
    if ((int) id.size() > 0) {
        string idp = phaseNode.id();
        if (idp != id) {
            throw CanteraError(subname,
                               "phasenode and Id are incompatible");
        }
    }

    /*
     * Check on the thermo field. Must have:
     * <thermo model="Redlich-Kister" />
     */
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError(subname, "no thermo XML node");
    }
    XML_Node& thermoNode = phaseNode.child("thermo");
    std::string mStringa = thermoNode.attrib("model");
    std::string mString = lowercase(mStringa);
    if (mString != "redlich-kister") {
        throw CanteraError(subname.c_str(),
                           "Unknown thermo model: " + mStringa + " - This object only knows \"Redlich-Kister\" ");
    }

    /*
     * Go get all of the coefficients and factors in the
     * activityCoefficients XML block
     */
    XML_Node* acNodePtr = 0;
    if (thermoNode.hasChild("activityCoefficients")) {
        XML_Node& acNode = thermoNode.child("activityCoefficients");
        acNodePtr = &acNode;
        std::string mStringa = acNode.attrib("model");
        std::string mString = lowercase(mStringa);
        if (mString != "redlich-kister") {
            throw CanteraError(subname.c_str(),
                               "Unknown activity coefficient model: " + mStringa);
        }
        size_t n = acNodePtr->nChildren();
        for (size_t i = 0; i < n; i++) {
            XML_Node& xmlACChild = acNodePtr->child(i);
            stemp = xmlACChild.name();
            std::string nodeName = lowercase(stemp);
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
    GibbsExcessVPSSTP::initThermoXML(phaseNode, id);
}
//===================================================================================================================
// Update the activity coefficients
/*
 * This function will be called to update the internally stored
 * natural logarithm of the activity coefficients
 *
 */
void RedlichKisterVPSSTP::s_update_lnActCoeff() const
{
    doublereal XA, XB;
    doublereal T = temperature();
    doublereal RT = GasConstant * T;

    lnActCoeff_Scaled_.assign(m_kk, 0.0);

    /*
     *  Scaling:  I moved the division of RT higher so that we are always dealing with G/RT dimensionless terms
     *            within the routine. There is a severe problem with roundoff error in these calculations. The
     *            dimensionless terms help.
     */

    for (size_t i = 0; i <  numBinaryInteractions_; i++) {
        size_t iA =  m_pSpecies_A_ij[i];
        size_t iB =  m_pSpecies_B_ij[i];
        XA = moleFractions_[iA];
        XB = moleFractions_[iB];
        doublereal deltaX = XA - XB;
        size_t N = m_N_ij[i];
        vector_fp& he_vec = m_HE_m_ij[i];
        vector_fp& se_vec = m_SE_m_ij[i];
        doublereal poly = 1.0;
        doublereal polyMm1 = 1.0;
        doublereal sum = 0.0;
        doublereal sumMm1 = 0.0;
        doublereal sum2 = 0.0;
        for (size_t m = 0; m < N; m++) {
            doublereal A_ge = (he_vec[m] -  T * se_vec[m]) / RT;
            sum += A_ge * poly;
            sum2 += A_ge * (m + 1) * poly;
            poly *= deltaX;
            if (m >= 1) {
                sumMm1 += (A_ge * polyMm1 * m);
                polyMm1 *= deltaX;
            }
        }
        doublereal oneMXA = 1.0 - XA;
        doublereal oneMXB = 1.0 - XB;
        for (size_t k = 0; k < m_kk; k++) {
            if (iA == k) {
                lnActCoeff_Scaled_[k] += (oneMXA * XB * sum) + (XA * XB * sumMm1 * (oneMXA + XB));
            } else  if (iB == k) {
                lnActCoeff_Scaled_[k] += (oneMXB * XA * sum) + (XA * XB * sumMm1 * (-oneMXB - XA));
            } else {
                lnActCoeff_Scaled_[k] += -(XA * XB * sum2);
            }
        }
        // Debug against formula in literature
#ifdef DEBUG_MODE_NOT
        double lnA = 0.0;
        double lnB = 0.0;
        double polyk = 1.0;
        double fac = 2.0 * XA - 1.0;
        for (int m = 0; m < N; m++) {
            doublereal A_ge = (he_vec[m] - T * se_vec[m]) / RT;
            lnA += A_ge * oneMXA * oneMXA * polyk * (1.0 + 2.0 * XA * m / fac);
            lnB += A_ge * XA * XA * polyk * (1.0 - 2.0 * oneMXA * m / fac);
            polyk *= fac;
        }
        // This gives the same result as above
        //  printf("RT lnActCoeff_Scaled_[iA] = %15.8E   , lnA = %15.8E\n",  lnActCoeff_Scaled_[iA], lnA);
        // printf("RT lnActCoeff_Scaled_[iB] = %15.8E   , lnB = %15.8E\n",  lnActCoeff_Scaled_[iB], lnB);

#endif

    }

}
//===================================================================================================================
// Update the derivative of the log of the activity coefficients wrt T
/*
 * This function will be called to update the internally stored
 * natural logarithm of the activity coefficients
 *

 */
void RedlichKisterVPSSTP::s_update_dlnActCoeff_dT() const
{
    doublereal XA, XB;
    //   doublereal T = temperature();

    dlnActCoeffdT_Scaled_.assign(m_kk, 0.0);
    d2lnActCoeffdT2_Scaled_.assign(m_kk, 0.0);

    for (size_t i = 0; i <  numBinaryInteractions_; i++) {
        size_t iA =  m_pSpecies_A_ij[i];
        size_t iB =  m_pSpecies_B_ij[i];
        XA = moleFractions_[iA];
        XB = moleFractions_[iB];
        doublereal deltaX = XA - XB;
        size_t N = m_N_ij[i];
        doublereal poly = 1.0;
        doublereal sum = 0.0;

        vector_fp& se_vec = m_SE_m_ij[i];
        doublereal sumMm1 = 0.0;
        doublereal polyMm1 = 1.0;
        doublereal sum2 = 0.0;
        for (size_t m = 0; m < N; m++) {
            doublereal A_ge = - se_vec[m];
            sum += A_ge * poly;
            sum2 += A_ge * (m + 1) * poly;
            poly *= deltaX;
            if (m >= 1) {
                sumMm1 += (A_ge * polyMm1 * m);
                polyMm1 *= deltaX;
            }
        }
        doublereal oneMXA = 1.0 - XA;
        doublereal oneMXB = 1.0 - XB;
        for (size_t k = 0; k < m_kk; k++) {
            if (iA == k) {
                dlnActCoeffdT_Scaled_[k] += (oneMXA * XB * sum) + (XA * XB * sumMm1 * (oneMXA + XB));
            } else  if (iB == k) {
                dlnActCoeffdT_Scaled_[k] += (oneMXB * XA * sum) + (XA * XB * sumMm1 * (-oneMXB - XA));
            } else {
                dlnActCoeffdT_Scaled_[k] += -(XA * XB * sum2);
            }
        }
    }
}
//====================================================================================================================
void RedlichKisterVPSSTP::getdlnActCoeffdT(doublereal* dlnActCoeffdT) const
{
    s_update_dlnActCoeff_dT();
    for (size_t k = 0; k < m_kk; k++) {
        dlnActCoeffdT[k] = dlnActCoeffdT_Scaled_[k];
    }
}
//====================================================================================================================
void RedlichKisterVPSSTP::getd2lnActCoeffdT2(doublereal* d2lnActCoeffdT2) const
{
    s_update_dlnActCoeff_dT();
    for (size_t k = 0; k < m_kk; k++) {
        d2lnActCoeffdT2[k] = d2lnActCoeffdT2_Scaled_[k];
    }
}
//====================================================================================================================
void RedlichKisterVPSSTP::s_update_dlnActCoeff_dX_() const
{
    doublereal XA, XB;
    doublereal T = temperature();

    dlnActCoeff_dX_.zero();

    for (size_t i = 0; i <  numBinaryInteractions_; i++) {
        size_t iA =  m_pSpecies_A_ij[i];
        size_t iB =  m_pSpecies_B_ij[i];
        XA = moleFractions_[iA];
        XB = moleFractions_[iB];
        doublereal deltaX = XA - XB;
        size_t N = m_N_ij[i];
        doublereal poly = 1.0;
        doublereal sum = 0.0;
        vector_fp& he_vec = m_HE_m_ij[i];
        vector_fp& se_vec = m_SE_m_ij[i];
        doublereal sumMm1 = 0.0;
        doublereal polyMm1 = 1.0;
        doublereal polyMm2 = 1.0;
        doublereal sum2 = 0.0;
        doublereal sum2Mm1 = 0.0;
        doublereal sumMm2 = 0.0;
        for (size_t m = 0; m < N; m++) {
            doublereal A_ge = he_vec[m] -  T * se_vec[m];
            sum += A_ge * poly;
            sum2 += A_ge * (m + 1) * poly;
            poly *= deltaX;
            if (m >= 1) {
                sumMm1  += (A_ge * polyMm1 * m);
                sum2Mm1 += (A_ge * polyMm1 * m * (1.0 + m));
                polyMm1 *= deltaX;
            }
            if (m >= 2) {
                sumMm2 += (A_ge * polyMm2 * m * (m - 1.0));
                polyMm2 *= deltaX;
            }
        }

        for (size_t k = 0; k < m_kk; k++) {
            if (iA == k) {

                dlnActCoeff_dX_(k, iA) += (- XB * sum + (1.0 - XA) * XB * sumMm1
                                           + XB * sumMm1 * (1.0 - 2.0 * XA + XB)
                                           + XA * XB * sumMm2 * (1.0 - XA + XB));

                dlnActCoeff_dX_(k, iB) += ((1.0 - XA) * sum - (1.0 - XA) * XB * sumMm1
                                           + XA * sumMm1 * (1.0 + 2.0 * XB - XA)
                                           - XA * XB * sumMm2 * (1.0 - XA + XB));

            } else  if (iB == k) {

                dlnActCoeff_dX_(k, iA) += ((1.0 - XB) * sum + (1.0 - XA) * XB * sumMm1
                                           + XB * sumMm1 * (1.0 - 2.0 * XA + XB)
                                           + XA * XB * sumMm2 * (1.0 - XA + XB));

                dlnActCoeff_dX_(k, iB) += (- XA * sum - (1.0 - XB) * XA * sumMm1
                                           + XA * sumMm1 * (XB - XA - (1.0 - XB))
                                           - XA * XB * sumMm2 * (-XA - (1.0 - XB)));
            } else {

                dlnActCoeff_dX_(k, iA) += (- XB * sum2  - XA * XB * sum2Mm1);

                dlnActCoeff_dX_(k, iB) += (- XA * sum2  + XA * XB * sum2Mm1);

            }
        }
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
void RedlichKisterVPSSTP::getdlnActCoeffds(const doublereal dTds, const doublereal* const dXds,
        doublereal* dlnActCoeffds) const
{
    s_update_dlnActCoeff_dT();
    s_update_dlnActCoeff_dX_();
    for (size_t k = 0; k < m_kk; k++) {
        dlnActCoeffds[k] = dlnActCoeffdT_Scaled_[k] * dTds;
        for (size_t l = 0; l < m_kk; l++) {
            dlnActCoeffds[k] += dlnActCoeff_dX_(k, l) * dXds[l];
        }
    }
}

//====================================================================================================================
void RedlichKisterVPSSTP::getdlnActCoeffdlnN_diag(doublereal* dlnActCoeffdlnN_diag) const
{
    s_update_dlnActCoeff_dX_();
    for (size_t l = 0; l < m_kk; l++) {
        dlnActCoeffdlnN_diag[l] = dlnActCoeff_dX_(l, l);
        for (size_t k = 0; k < m_kk; k++) {
            dlnActCoeffdlnN_diag[k] -= dlnActCoeff_dX_(l, k) * moleFractions_[k];
        }
    }
}
//====================================================================================================================
void RedlichKisterVPSSTP::getdlnActCoeffdlnX_diag(doublereal* dlnActCoeffdlnX_diag) const
{
    s_update_dlnActCoeff_dX_();
    for (size_t k = 0; k < m_kk; k++) {
        dlnActCoeffdlnX_diag[k] = dlnActCoeffdlnX_diag_[k];
    }
}
//====================================================================================================================
void RedlichKisterVPSSTP::getdlnActCoeffdlnN(const size_t ld, doublereal* dlnActCoeffdlnN)
{
    s_update_dlnActCoeff_dX_();
    double* data =  & dlnActCoeffdlnN_(0,0);
    for (size_t k = 0; k < m_kk; k++) {
        for (size_t m = 0; m < m_kk; m++) {
            dlnActCoeffdlnN[ld * k + m] = data[m_kk * k + m];
        }
    }
}
//====================================================================================================================
void RedlichKisterVPSSTP::resizeNumInteractions(const size_t num)
{
    numBinaryInteractions_ = num;
    m_pSpecies_A_ij.resize(num, npos);
    m_pSpecies_B_ij.resize(num, npos);
    m_N_ij.resize(num, npos);
    m_HE_m_ij.resize(num);
    m_SE_m_ij.resize(num);
    dlnActCoeff_dX_.resize(num, num);
}
//====================================================================================================================
// Process an XML node called "binaryNeutralSpeciesParameters"
/*
 *  This node contains all of the parameters necessary to describe the RedlichKister Interaction for
 *  a single binary interaction. This function reads the XML file and writes the coefficients
 *  it finds to an internal data structures.
 */
void RedlichKisterVPSSTP::readXMLBinarySpecies(XML_Node& xmLBinarySpecies)
{
    std::string xname = xmLBinarySpecies.name();
    if (xname != "binaryNeutralSpeciesParameters") {
        throw CanteraError("RedlichKisterVPSSTP::readXMLBinarySpecies",
                           "Incorrect name for processing this routine: " + xname);
    }
    double* charge = DATA_PTR(m_speciesCharge);
    std::string stemp;
    size_t Npoly = 0;
    vector_fp hParams, sParams, vParams;
    std::string iName = xmLBinarySpecies.attrib("speciesA");
    if (iName == "") {
        throw CanteraError("RedlichKisterVPSSTP::readXMLBinarySpecies", "no speciesA attrib");
    }
    std::string jName = xmLBinarySpecies.attrib("speciesB");
    if (jName == "") {
        throw CanteraError("RedlichKisterVPSSTP::readXMLBinarySpecies", "no speciesB attrib");
    }
    /*
     * Find the index of the species in the current phase. It's not
     * an error to not find the species. This means that the interaction doesn't occur for the current
     * implementation of the phase.
     */
    size_t iSpecies = speciesIndex(iName);
    if (iSpecies == npos) {
        return;
    }
    string ispName = speciesName(iSpecies);
    if (charge[iSpecies] != 0) {
        throw CanteraError("RedlichKisterVPSSTP::readXMLBinarySpecies", "speciesA charge problem");
    }
    size_t jSpecies = speciesIndex(jName);
    if (jSpecies == npos) {
        return;
    }
    std::string jspName = speciesName(jSpecies);
    if (charge[jSpecies] != 0) {
        throw CanteraError("RedlichKisterVPSSTP::readXMLBinarySpecies", "speciesB charge problem");
    }
    /*
     *  Ok we have found a valid interaction
     */
    numBinaryInteractions_++;
    size_t iSpot = numBinaryInteractions_ - 1;
    m_pSpecies_A_ij.resize(numBinaryInteractions_);
    m_pSpecies_B_ij.resize(numBinaryInteractions_);
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
            ctml::getFloatArray(xmlChild, hParams, true, "toSI", "excessEnthalpy");
            size_t nParamsFound = hParams.size();
            if (nParamsFound > Npoly) {
                Npoly = nParamsFound;
            }

        }

        if (nodeName == "excessentropy") {
            /*
             * Get the string containing all of the values
             */
            ctml::getFloatArray(xmlChild, sParams, true, "toSI", "excessEntropy");
            size_t nParamsFound = sParams.size();
            if (nParamsFound > Npoly) {
                Npoly = nParamsFound;
            }
        }
    }
    hParams.resize(Npoly, 0.0);
    sParams.resize(Npoly, 0.0);
    m_HE_m_ij.push_back(hParams);
    m_SE_m_ij.push_back(sParams);
    m_N_ij.push_back(Npoly);
    resizeNumInteractions(numBinaryInteractions_);
}
//====================================================================================================================
#ifdef DEBUG_MODE
void RedlichKisterVPSSTP::Vint(double& VintOut, double& voltsOut)
{
    int iA, iB, m;
    doublereal XA, XB;
    doublereal T = temperature();
    doublereal RT = GasConstant * T;
    double Volts = 0.0;

    lnActCoeff_Scaled_.assign(m_kk, 0.0);

    for (int i = 0; i < (int) numBinaryInteractions_; i++) {
        iA =  m_pSpecies_A_ij[i];
        iB =  m_pSpecies_B_ij[i];
        XA = moleFractions_[iA];
        XB = moleFractions_[iB];
        if (XA <= 1.0E-14) {
            XA = 1.0E-14;
        }
        if (XA >= (1.0 - 1.0E-14)) {
            XA = 1.0 - 1.0E-14;
        }

        int N = m_N_ij[i];
        vector_fp& he_vec = m_HE_m_ij[i];
        vector_fp& se_vec = m_SE_m_ij[i];
        double fac = 2.0 * XA - 1.0;
        if (fabs(fac) < 1.0E-13) {
            fac = 1.0E-13;
        }
        double polykp1 = fac;
        double poly1mk = fac;

        for (m = 0; m < N; m++) {
            doublereal A_ge = he_vec[m] - T * se_vec[m];
            Volts += A_ge * (polykp1 - (2.0 * XA * m * (1.0-XA)) / poly1mk);
            polykp1 *= fac;
            poly1mk /= fac;
        }
    }
    Volts /= Faraday;

    double termp = RT * log((1.0 - XA)/XA) / Faraday;

    VintOut =  Volts;
    voltsOut = Volts + termp;
}
#endif
//====================================================================================================================
}

