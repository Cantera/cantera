/**
 *  @file DebyeHuckel.cpp
 *    Declarations for the DebyeHuckel ThermoPhase object, which models dilute
 *    electrolyte solutions
 *    (see \ref thermoprops and \link Cantera::DebyeHuckel DebyeHuckel \endlink).
 *
 * Class DebyeHuckel represents a dilute liquid electrolyte phase which
 * obeys the Debye Huckel formulation for nonideality.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/DebyeHuckel.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/PDSS_Water.h"
#include "cantera/thermo/electrolytes.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"

#include <cstdio>

using namespace std;

namespace Cantera
{

DebyeHuckel::DebyeHuckel() :
    m_formDH(DHFORM_DILUTE_LIMIT),
    m_formGC(2),
    m_IionicMolality(0.0),
    m_maxIionicStrength(30.0),
    m_useHelgesonFixedForm(false),
    m_IionicMolalityStoich(0.0),
    m_form_A_Debye(A_DEBYE_CONST),
    m_A_Debye(1.172576), // units = sqrt(kg/gmol)
    m_B_Debye(3.28640E9), // units = sqrt(kg/gmol) / m
    m_waterSS(0),
    m_densWaterSS(1000.)
{
}

DebyeHuckel::DebyeHuckel(const std::string& inputFile,
                         const std::string& id_) :
    m_formDH(DHFORM_DILUTE_LIMIT),
    m_formGC(2),
    m_IionicMolality(0.0),
    m_maxIionicStrength(30.0),
    m_useHelgesonFixedForm(false),
    m_IionicMolalityStoich(0.0),
    m_form_A_Debye(A_DEBYE_CONST),
    m_A_Debye(1.172576), // units = sqrt(kg/gmol)
    m_B_Debye(3.28640E9), // units = sqrt(kg/gmol) / m
    m_waterSS(0),
    m_densWaterSS(1000.)
{
    initThermoFile(inputFile, id_);
}

DebyeHuckel::DebyeHuckel(XML_Node& phaseRoot, const std::string& id_) :
    m_formDH(DHFORM_DILUTE_LIMIT),
    m_formGC(2),
    m_IionicMolality(0.0),
    m_maxIionicStrength(3.0),
    m_useHelgesonFixedForm(false),
    m_IionicMolalityStoich(0.0),
    m_form_A_Debye(A_DEBYE_CONST),
    m_A_Debye(1.172576), // units = sqrt(kg/gmol)
    m_B_Debye(3.28640E9), // units = sqrt(kg/gmol) / m
    m_waterSS(0),
    m_densWaterSS(1000.)
{
    importPhase(phaseRoot, this);
}

DebyeHuckel::DebyeHuckel(const DebyeHuckel& b) :
    m_formDH(DHFORM_DILUTE_LIMIT),
    m_formGC(2),
    m_IionicMolality(0.0),
    m_maxIionicStrength(30.0),
    m_useHelgesonFixedForm(false),
    m_IionicMolalityStoich(0.0),
    m_form_A_Debye(A_DEBYE_CONST),
    m_A_Debye(1.172576), // units = sqrt(kg/gmol)
    m_B_Debye(3.28640E9), // units = sqrt(kg/gmol) / m
    m_waterSS(0),
    m_densWaterSS(1000.)
{
    // Use the assignment operator to do the brunt of the work for the copy
    // constructor.
    *this = b;
}

DebyeHuckel& DebyeHuckel::operator=(const DebyeHuckel& b)
{
    if (&b != this) {
        MolalityVPSSTP::operator=(b);
        m_formDH = b.m_formDH;
        m_formGC = b.m_formGC;
        m_Aionic = b.m_Aionic;
        m_IionicMolality = b.m_IionicMolality;
        m_maxIionicStrength = b.m_maxIionicStrength;
        m_useHelgesonFixedForm= b.m_useHelgesonFixedForm;
        m_IionicMolalityStoich= b.m_IionicMolalityStoich;
        m_form_A_Debye = b.m_form_A_Debye;
        m_A_Debye = b.m_A_Debye;
        m_B_Debye = b.m_B_Debye;
        m_B_Dot = b.m_B_Dot;

        // This is an internal shallow copy of the PDSS_Water pointer
        m_waterSS = dynamic_cast<PDSS_Water*>(providePDSS(0));
        if (!m_waterSS) {
            throw CanteraError("DebyeHuckel::operator=()", "Dynamic cast to waterPDSS failed");
        }

        m_densWaterSS = b.m_densWaterSS;

        if (b.m_waterProps) {
            m_waterProps.reset(new WaterProps(m_waterSS));
        }

        m_tmpV = b.m_tmpV;
        m_speciesCharge_Stoich= b.m_speciesCharge_Stoich;
        m_Beta_ij = b.m_Beta_ij;
        m_lnActCoeffMolal = b.m_lnActCoeffMolal;
        m_d2lnActCoeffMolaldT2= b.m_d2lnActCoeffMolaldT2;
    }
    return *this;
}

DebyeHuckel::~DebyeHuckel()
{
}

ThermoPhase* DebyeHuckel::duplMyselfAsThermoPhase() const
{
    return new DebyeHuckel(*this);
}

int DebyeHuckel::eosType() const
{
    warn_deprecated("DebyeHuckel::eosType",
                    "To be removed after Cantera 2.3.");
    int res;
    switch (m_formGC) {
    case 0:
        res = cDebyeHuckel0;
        break;
    case 1:
        res = cDebyeHuckel1;
        break;
    case 2:
        res = cDebyeHuckel2;
        break;
    default:
        throw CanteraError("eosType", "Unknown type");
        break;
    }
    return res;
}

// -------- Molar Thermodynamic Properties of the Solution ---------------

doublereal DebyeHuckel::enthalpy_mole() const
{
    getPartialMolarEnthalpies(m_tmpV.data());
    return mean_X(m_tmpV);
}

doublereal DebyeHuckel::entropy_mole() const
{
    getPartialMolarEntropies(m_tmpV.data());
    return mean_X(m_tmpV);
}

doublereal DebyeHuckel::gibbs_mole() const
{
    getChemPotentials(m_tmpV.data());
    return mean_X(m_tmpV);
}

doublereal DebyeHuckel::cp_mole() const
{
    getPartialMolarCp(m_tmpV.data());
    return mean_X(m_tmpV);
}

// ------- Mechanical Equation of State Properties ------------------------

void DebyeHuckel::calcDensity()
{
    if (m_waterSS) {
        // Store the internal density of the water SS. Note, we would have to do
        // this for all other species if they had pressure dependent properties.
        m_densWaterSS = m_waterSS->density();
    }
    getPartialMolarVolumes(m_tmpV.data());
    double dd = meanMolecularWeight() / mean_X(m_tmpV);
    Phase::setDensity(dd);
}

void DebyeHuckel::setDensity(doublereal rho)
{
    double dens = density();
    if (rho != dens) {
        throw CanteraError("Idea;MolalSoln::setDensity",
                           "Density is not an independent variable");
    }
}

void DebyeHuckel::setMolarDensity(const doublereal conc)
{
    double concI = molarDensity();
    if (conc != concI) {
        throw CanteraError("Idea;MolalSoln::setMolarDensity",
                           "molarDensity/density is not an independent variable");
    }
}

// ------- Activities and Activity Concentrations

void DebyeHuckel::getActivityConcentrations(doublereal* c) const
{
    double c_solvent = standardConcentration();
    getActivities(c);
    for (size_t k = 0; k < m_kk; k++) {
        c[k] *= c_solvent;
    }
}

doublereal DebyeHuckel::standardConcentration(size_t k) const
{
    double mvSolvent = m_speciesSize[m_indexSolvent];
    return 1.0 / mvSolvent;
}

void DebyeHuckel::getActivities(doublereal* ac) const
{
    _updateStandardStateThermo();

    // Update the molality array, m_molalities(). This requires an update due to
    // mole fractions
    s_update_lnMolalityActCoeff();
    for (size_t k = 0; k < m_kk; k++) {
        if (k != m_indexSolvent) {
            ac[k] = m_molalities[k] * exp(m_lnActCoeffMolal[k]);
        }
    }
    double xmolSolvent = moleFraction(m_indexSolvent);
    ac[m_indexSolvent] =
        exp(m_lnActCoeffMolal[m_indexSolvent]) * xmolSolvent;
}

void DebyeHuckel::getMolalityActivityCoefficients(doublereal* acMolality) const
{
    _updateStandardStateThermo();
    A_Debye_TP(-1.0, -1.0);
    s_update_lnMolalityActCoeff();
    copy(m_lnActCoeffMolal.begin(), m_lnActCoeffMolal.end(), acMolality);
    for (size_t k = 0; k < m_kk; k++) {
        acMolality[k] = exp(acMolality[k]);
    }
}

// ------ Partial Molar Properties of the Solution -----------------

void DebyeHuckel::getChemPotentials(doublereal* mu) const
{
    double xx;

    // First get the standard chemical potentials in molar form. This requires
    // updates of standard state as a function of T and P
    getStandardChemPotentials(mu);

    // Update the activity coefficients. This also updates the internal molality
    // array.
    s_update_lnMolalityActCoeff();
    double xmolSolvent = moleFraction(m_indexSolvent);
    for (size_t k = 0; k < m_kk; k++) {
        if (m_indexSolvent != k) {
            xx = std::max(m_molalities[k], SmallNumber);
            mu[k] += RT() * (log(xx) + m_lnActCoeffMolal[k]);
        }
    }
    xx = std::max(xmolSolvent, SmallNumber);
    mu[m_indexSolvent] +=
        RT() * (log(xx) + m_lnActCoeffMolal[m_indexSolvent]);
}

void DebyeHuckel::getPartialMolarEnthalpies(doublereal* hbar) const
{
    // Get the nondimensional standard state enthalpies
    getEnthalpy_RT(hbar);

    // Dimensionalize it.
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] *= RT();
    }

    // Check to see whether activity coefficients are temperature
    // dependent. If they are, then calculate the their temperature
    // derivatives and add them into the result.
    double dAdT = dA_DebyedT_TP();
    if (dAdT != 0.0) {
        // Update the activity coefficients, This also update the
        // internally stored molalities.
        s_update_lnMolalityActCoeff();
        s_update_dlnMolalityActCoeff_dT();
        for (size_t k = 0; k < m_kk; k++) {
            hbar[k] -= RT() * temperature() * m_dlnActCoeffMolaldT[k];
        }
    }
}

void DebyeHuckel::getPartialMolarEntropies(doublereal* sbar) const
{
    // Get the standard state entropies at the temperature and pressure of the
    // solution.
    getEntropy_R(sbar);

    // Dimensionalize the entropies
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] *= GasConstant;
    }

    // Update the activity coefficients, This also update the internally stored
    // molalities.
    s_update_lnMolalityActCoeff();

    // First we will add in the obvious dependence on the T term out front of
    // the log activity term
    doublereal mm;
    for (size_t k = 0; k < m_kk; k++) {
        if (k != m_indexSolvent) {
            mm = std::max(SmallNumber, m_molalities[k]);
            sbar[k] -= GasConstant * (log(mm) + m_lnActCoeffMolal[k]);
        }
    }
    double xmolSolvent = moleFraction(m_indexSolvent);
    mm = std::max(SmallNumber, xmolSolvent);
    sbar[m_indexSolvent] -= GasConstant *(log(mm) + m_lnActCoeffMolal[m_indexSolvent]);

    // Check to see whether activity coefficients are temperature dependent. If
    // they are, then calculate the their temperature derivatives and add them
    // into the result.
    double dAdT = dA_DebyedT_TP();
    if (dAdT != 0.0) {
        s_update_dlnMolalityActCoeff_dT();
        for (size_t k = 0; k < m_kk; k++) {
            sbar[k] -= RT() * m_dlnActCoeffMolaldT[k];
        }
    }
}

void DebyeHuckel::getPartialMolarVolumes(doublereal* vbar) const
{
    getStandardVolumes(vbar);

    // Update the derivatives wrt the activity coefficients.
    s_update_lnMolalityActCoeff();
    s_update_dlnMolalityActCoeff_dP();
    for (size_t k = 0; k < m_kk; k++) {
        vbar[k] += RT() * m_dlnActCoeffMolaldP[k];
    }
}

void DebyeHuckel::getPartialMolarCp(doublereal* cpbar) const
{
    getCp_R(cpbar);
    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] *= GasConstant;
    }

    // Check to see whether activity coefficients are temperature dependent. If
    // they are, then calculate the their temperature derivatives and add them
    // into the result.
    double dAdT = dA_DebyedT_TP();
    if (dAdT != 0.0) {
        // Update the activity coefficients, This also update the internally
        // stored molalities.
        s_update_lnMolalityActCoeff();
        s_update_dlnMolalityActCoeff_dT();
        s_update_d2lnMolalityActCoeff_dT2();
        for (size_t k = 0; k < m_kk; k++) {
            cpbar[k] -= (2.0 * RT() * m_dlnActCoeffMolaldT[k] +
                         RT() * temperature() * m_d2lnActCoeffMolaldT2[k]);
        }
    }
}

// -------------- Utilities -------------------------------

//! Utility function to assign an integer value from a string for the
//! ElectrolyteSpeciesType field.
/*!
 *  @param estString  input string that will be interpreted
 */
static int interp_est(const std::string& estString)
{
    if (ba::iequals(estString, "solvent")) {
        return cEST_solvent;
    } else if (ba::iequals(estString, "chargedspecies")) {
        return cEST_chargedSpecies;
    } else if (ba::iequals(estString, "weakacidassociated")) {
        return cEST_weakAcidAssociated;
    } else if (ba::iequals(estString, "strongacidassociated")) {
        return cEST_strongAcidAssociated;
    } else if (ba::iequals(estString, "polarneutral")) {
        return cEST_polarNeutral;
    } else if (ba::iequals(estString, "nonpolarneutral")) {
        return cEST_nonpolarNeutral;
    }
    int retn, rval;
    if ((retn = sscanf(estString.c_str(), "%d", &rval)) != 1) {
        return -1;
    }
    return rval;
}

void DebyeHuckel::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    if (id_.size() > 0) {
        std::string idp = phaseNode.id();
        if (idp != id_) {
            throw CanteraError("DebyeHuckel::initThermoXML",
                               "phasenode and Id are incompatible");
        }
    }

    // Find the Thermo XML node
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError("DebyeHuckel::initThermoXML",
                           "no thermo XML node");
    }
    XML_Node& thermoNode = phaseNode.child("thermo");

    // Determine the form of the Debye-Huckel model, m_formDH.  We will use this
    // information to size arrays below.
    if (thermoNode.hasChild("activityCoefficients")) {
        XML_Node& scNode = thermoNode.child("activityCoefficients");
        m_formDH = DHFORM_DILUTE_LIMIT;
        std::string formString = scNode.attrib("model");
        if (formString != "") {
            if (formString == "Dilute_limit") {
                m_formDH = DHFORM_DILUTE_LIMIT;
            } else if (formString == "Bdot_with_variable_a") {
                m_formDH = DHFORM_BDOT_AK;
            } else if (formString == "Bdot_with_common_a") {
                m_formDH = DHFORM_BDOT_ACOMMON;
            } else if (formString == "Beta_ij") {
                m_formDH = DHFORM_BETAIJ;
            } else if (formString == "Pitzer_with_Beta_ij") {
                m_formDH = DHFORM_PITZER_BETAIJ;
            } else {
                throw CanteraError("DebyeHuckel::initThermoXML",
                                   "Unknown standardConc model: " + formString);
            }
        }
    } else {
        // If there is no XML node named "activityCoefficients", assume
        // that we are doing the extreme dilute limit assumption
        m_formDH = DHFORM_DILUTE_LIMIT;
    }

    // Possibly change the form of the standard concentrations
    if (thermoNode.hasChild("standardConc")) {
        XML_Node& scNode = thermoNode.child("standardConc");
        m_formGC = 2;
        std::string formString = scNode.attrib("model");
        if (formString != "") {
            if (formString == "unity") {
                m_formGC = 0;
                throw CanteraError("DebyeHuckel::initThermoXML",
                                   "standardConc = unity not done");
            } else if (formString == "molar_volume") {
                m_formGC = 1;
                throw CanteraError("DebyeHuckel::initThermoXML",
                                   "standardConc = molar_volume not done");
            } else if (formString == "solvent_volume") {
                m_formGC = 2;
            } else {
                throw CanteraError("DebyeHuckel::initThermoXML",
                                   "Unknown standardConc model: " + formString);
            }
        }
    }

    // Reconcile the solvent name and index.

    // Get the Name of the Solvent:
    //      <solvent> solventName </solvent>
    std::string solventName = "";
    if (thermoNode.hasChild("solvent")) {
        XML_Node& scNode = thermoNode.child("solvent");
        vector<std::string> nameSolventa;
        getStringArray(scNode, nameSolventa);
        if (nameSolventa.size() != 1) {
            throw CanteraError("DebyeHuckel::initThermoXML",
                               "badly formed solvent XML node");
        }
        solventName = nameSolventa[0];
    }
    for (size_t k = 0; k < m_kk; k++) {
        std::string sname = speciesName(k);
        if (solventName == sname) {
            m_indexSolvent = k;
            break;
        }
    }
    if (m_indexSolvent == npos) {
        cout << "DebyeHuckel::initThermoXML: Solvent Name not found"
             << endl;
        throw CanteraError("DebyeHuckel::initThermoXML",
                           "Solvent name not found");
    }
    if (m_indexSolvent != 0) {
        throw CanteraError("DebyeHuckel::initThermoXML",
                           "Solvent " + solventName +
                           " should be first species");
    }

    // Now go get the specification of the standard states for species in the
    // solution. This includes the molar volumes data blocks for incompressible
    // species.
    XML_Node& speciesList = phaseNode.child("speciesArray");
    XML_Node* speciesDB =
        get_XML_NameID("speciesData", speciesList["datasrc"],
                       &phaseNode.root());
    const vector<string>&sss = speciesNames();

    for (size_t k = 0; k < m_kk; k++) {
        XML_Node* s = speciesDB->findByAttr("name", sss[k]);
        if (!s) {
            throw CanteraError("DebyeHuckel::initThermoXML",
                               "Species Data Base " + sss[k] + " not found");
        }
        XML_Node* ss = s->findByName("standardState");
        if (!ss) {
            throw CanteraError("DebyeHuckel::initThermoXML",
                               "Species " + sss[k] +
                               " standardState XML block not found");
        }
        std::string modelString = ss->attrib("model");
        if (modelString == "") {
            throw CanteraError("DebyeHuckel::initThermoXML",
                               "Species " + sss[k] +
                               " standardState XML block model attribute not found");
        }

        if (k == 0) {
            if (ba::iequals(modelString, "wateriapws") || ba::iequals(modelString, "real_water") ||
                ba::iequals(modelString, "waterpdss")) {
                // Initialize the water standard state model
                m_waterSS = dynamic_cast<PDSS_Water*>(providePDSS(0));
                if (!m_waterSS) {
                    throw CanteraError("HMWSoln::installThermoXML",
                                       "Dynamic cast to PDSS_Water failed");
                }

                // Fill in the molar volume of water (m3/kmol) at standard
                // conditions to fill in the m_speciesSize entry with something
                // reasonable.
                m_waterSS->setState_TP(300., OneAtm);
                double dens = m_waterSS->density();
                double mw = m_waterSS->molecularWeight();
                m_speciesSize[0] = mw / dens;
            } else if (ba::iequals(modelString, "constant_incompressible")) {
                m_speciesSize[k] = getFloat(*ss, "molarVolume", "toSi");
            } else {
                throw CanteraError("DebyeHuckel::initThermoXML",
                                   "Solvent SS Model \"" + modelString +
                                   "\" is not known");
            }
        } else {
            if (!ba::iequals(modelString, "constant_incompressible")) {
                throw CanteraError("DebyeHuckel::initThermoXML",
                                   "Solute SS Model \"" + modelString +
                                   "\" is not known");
            }
            m_speciesSize[k] = getFloat(*ss, "molarVolume", "toSI");
        }
    }

    // Go get all of the coefficients and factors in the activityCoefficients
    // XML block
    XML_Node* acNodePtr = 0;
    if (thermoNode.hasChild("activityCoefficients")) {
        XML_Node& acNode = thermoNode.child("activityCoefficients");
        acNodePtr = &acNode;

        // Look for parameters for A_Debye
        if (acNode.hasChild("A_Debye")) {
            XML_Node* ss = acNode.findByName("A_Debye");
            string modelString = ss->attrib("model");
            if (modelString != "") {
                if (ba::iequals(modelString, "water")) {
                    m_form_A_Debye = A_DEBYE_WATER;
                } else {
                    throw CanteraError("DebyeHuckel::initThermoXML",
                                       "A_Debye Model \"" + modelString +
                                       "\" is not known");
                }
            } else {
                m_A_Debye = getFloat(acNode, "A_Debye");
            }
        }

        // Initialize the water property calculator. It will share the internal
        // eos water calculator.
        if (m_form_A_Debye == A_DEBYE_WATER) {
            m_waterProps.reset(new WaterProps(m_waterSS));
        }

        // Look for parameters for B_Debye
        if (acNode.hasChild("B_Debye")) {
            m_B_Debye = getFloat(acNode, "B_Debye");
        }

        // Look for parameters for B_dot
        if (acNode.hasChild("B_dot")) {
            if (m_formDH == DHFORM_BETAIJ ||
                    m_formDH == DHFORM_DILUTE_LIMIT ||
                    m_formDH == DHFORM_PITZER_BETAIJ) {
                throw CanteraError("DebyeHuckel:init",
                                   "B_dot entry in the wrong DH form");
            }
            double bdot_common = getFloat(acNode, "B_dot");
            // Set B_dot parameters for charged species
            for (size_t k = 0; k < m_kk; k++) {
                double z_k = charge(k);
                if (fabs(z_k) > 0.0001) {
                    m_B_Dot[k] = bdot_common;
                } else {
                    m_B_Dot[k] = 0.0;
                }
            }
        }

        // Look for Parameters for the Maximum Ionic Strength
        if (acNode.hasChild("maxIonicStrength")) {
            m_maxIionicStrength = getFloat(acNode, "maxIonicStrength");
        }

        // Look for Helgeson Parameters
        if (acNode.hasChild("UseHelgesonFixedForm")) {
            m_useHelgesonFixedForm = true;
        } else {
            m_useHelgesonFixedForm = false;
        }

        // Look for parameters for the Ionic radius
        if (acNode.hasChild("ionicRadius")) {
            XML_Node& irNode = acNode.child("ionicRadius");

            double Afactor = 1.0;
            if (irNode.hasAttrib("units")) {
                std::string Aunits = irNode.attrib("units");
                Afactor = toSI(Aunits);
            }

            if (irNode.hasAttrib("default")) {
                std::string ads = irNode.attrib("default");
                double ad = fpValue(ads);
                for (size_t k = 0; k < m_kk; k++) {
                    m_Aionic[k] = ad * Afactor;
                }
            }

            // If the Debye-Huckel form is BDOT_AK, we can have separate values
            // for the denominator's ionic size. -> That's how the activity
            // coefficient is parameterized. In this case only do we allow the
            // code to read in these parameters.
            if (m_formDH == DHFORM_BDOT_AK) {
                // Define a string-string map, and interpret the value of the
                // XML element as binary pairs separated by colons, e.g.:
                //      Na+:3.0
                //      Cl-:4.0
                //      H+:9.0
                //      OH-:3.5
                // Read them into the map.
                map<string, string> m;
                getMap(irNode, m);

                // Iterate over the map pairs, interpreting the first string as
                // a species in the current phase. If no match is made, silently
                // ignore the lack of agreement (HKM -> may be changed in the
                // future).
                for (const auto& b : m) {
                    size_t kk = speciesIndex(b.first);
                    m_Aionic[kk] = fpValue(b.second) * Afactor;
                }
            }
        }

        // Get the matrix of coefficients for the Beta binary interaction
        // parameters. We assume here that this matrix is symmetric, so that we
        // only have to input 1/2 of the values.
        if (acNode.hasChild("DHBetaMatrix")) {
            if (m_formDH == DHFORM_BETAIJ ||
                    m_formDH == DHFORM_PITZER_BETAIJ) {
                m_Beta_ij.resize(m_kk, m_kk, 0.0);
                XML_Node& irNode = acNode.child("DHBetaMatrix");
                const vector<string>& sn = speciesNames();
                getMatrixValues(irNode, sn, sn, m_Beta_ij, true, true);
            } else {
                throw CanteraError("DebyeHuckel::initThermoXML:",
                                   "DHBetaMatrix found for wrong type");
            }
        }

        // Fill in parameters for the calculation of the stoichiometric Ionic
        // Strength. The default is that stoich charge is the same as the
        // regular charge.
        m_speciesCharge_Stoich.resize(m_kk, 0.0);
        for (size_t k = 0; k < m_kk; k++) {
            m_speciesCharge_Stoich[k] = m_speciesCharge[k];
        }

        // First look at the species database. Look for the subelement
        // "stoichIsMods" in each of the species SS databases.
        std::vector<const XML_Node*> xspecies= speciesData();
        size_t jj = xspecies.size();
        for (size_t k = 0; k < m_kk; k++) {
            size_t jmap = npos;
            std::string kname = speciesName(k);
            for (size_t j = 0; j < jj; j++) {
                const XML_Node& sp = *xspecies[j];
                std::string jname = sp["name"];
                if (jname == kname) {
                    jmap = j;
                    break;
                }
            }
            if (jmap != npos) {
                const XML_Node& sp = *xspecies[jmap];
                if (sp.hasChild("stoichIsMods")) {
                    double val = getFloat(sp, "stoichIsMods");
                    m_speciesCharge_Stoich[k] = val;
                }
            }
        }

        // Now look at the activity coefficient database
        if (acNodePtr && acNodePtr->hasChild("stoichIsMods")) {
            XML_Node& sIsNode = acNodePtr->child("stoichIsMods");
            map<std::string, std::string> msIs;
            getMap(sIsNode, msIs);
            for (const auto& b : msIs) {
                size_t kk = speciesIndex(b.first);
                double val = fpValue(b.second);
                m_speciesCharge_Stoich[kk] = val;
            }
        }
    }

    // Fill in the vector specifying the electrolyte species type. First fill in
    // default values. Everything is either a charge species, a nonpolar
    // neutral, or the solvent.
    for (size_t k = 0; k < m_kk; k++) {
        if (fabs(m_speciesCharge[k]) > 0.0001) {
            m_electrolyteSpeciesType[k] = cEST_chargedSpecies;
            if (fabs(m_speciesCharge_Stoich[k] - m_speciesCharge[k]) > 0.0001) {
                m_electrolyteSpeciesType[k] = cEST_weakAcidAssociated;
            }
        } else if (fabs(m_speciesCharge_Stoich[k]) > 0.0001) {
            m_electrolyteSpeciesType[k] = cEST_weakAcidAssociated;
        } else {
            m_electrolyteSpeciesType[k] = cEST_nonpolarNeutral;
        }
    }
    m_electrolyteSpeciesType[m_indexSolvent] = cEST_solvent;

    // First look at the species database. Look for the subelement
    // "stoichIsMods" in each of the species SS databases.
    std::vector<const XML_Node*> xspecies= speciesData();
    for (size_t k = 0; k < m_kk; k++) {
        std::string kname = speciesName(k);
        const XML_Node* spPtr = xspecies[k];
        if (spPtr && spPtr->hasChild("electrolyteSpeciesType")) {
            std::string est = getChildValue(*spPtr, "electrolyteSpeciesType");
            if ((m_electrolyteSpeciesType[k] = interp_est(est)) == -1) {
                throw CanteraError("DebyeHuckel:initThermoXML",
                                   "Bad electrolyte type: " + est);
            }
        }
    }

    // Then look at the phase thermo specification
    if (acNodePtr && acNodePtr->hasChild("electrolyteSpeciesType")) {
        XML_Node& ESTNode = acNodePtr->child("electrolyteSpeciesType");
        map<std::string, std::string> msEST;
        getMap(ESTNode, msEST);
        for (const auto& b : msEST) {
            size_t kk = speciesIndex(b.first);
            std::string est = b.second;
            if ((m_electrolyteSpeciesType[kk] = interp_est(est))  == -1) {
                throw CanteraError("DebyeHuckel:initThermoXML",
                                   "Bad electrolyte type: " + est);
            }
        }
    }

    // Lastly set the state
    if (phaseNode.hasChild("state")) {
        XML_Node& stateNode = phaseNode.child("state");
        setStateFromXML(stateNode);
    }
}

double DebyeHuckel::A_Debye_TP(double tempArg, double presArg) const
{
    double T = temperature();
    double A;
    if (tempArg != -1.0) {
        T = tempArg;
    }
    double P = pressure();
    if (presArg != -1.0) {
        P = presArg;
    }

    switch (m_form_A_Debye) {
    case A_DEBYE_CONST:
        A = m_A_Debye;
        break;
    case A_DEBYE_WATER:
        A = m_waterProps->ADebye(T, P, 0);
        m_A_Debye = A;
        break;
    default:
        throw CanteraError("DebyeHuckel::A_Debye_TP", "shouldn't be here");
    }
    return A;
}

double DebyeHuckel::dA_DebyedT_TP(double tempArg, double presArg) const
{
    double T = temperature();
    if (tempArg != -1.0) {
        T = tempArg;
    }
    double P = pressure();
    if (presArg != -1.0) {
        P = presArg;
    }
    double dAdT;
    switch (m_form_A_Debye) {
    case A_DEBYE_CONST:
        dAdT = 0.0;
        break;
    case A_DEBYE_WATER:
        dAdT = m_waterProps->ADebye(T, P, 1);
        break;
    default:
        throw CanteraError("DebyeHuckel::dA_DebyedT_TP", "shouldn't be here");
    }
    return dAdT;
}

double DebyeHuckel::d2A_DebyedT2_TP(double tempArg, double presArg) const
{
    double T = temperature();
    if (tempArg != -1.0) {
        T = tempArg;
    }
    double P = pressure();
    if (presArg != -1.0) {
        P = presArg;
    }
    double d2AdT2;
    switch (m_form_A_Debye) {
    case A_DEBYE_CONST:
        d2AdT2 = 0.0;
        break;
    case A_DEBYE_WATER:
        d2AdT2 = m_waterProps->ADebye(T, P, 2);
        break;
    default:
        throw CanteraError("DebyeHuckel::d2A_DebyedT2_TP", "shouldn't be here");
    }
    return d2AdT2;
}

double DebyeHuckel::dA_DebyedP_TP(double tempArg, double presArg) const
{
    double T = temperature();
    if (tempArg != -1.0) {
        T = tempArg;
    }
    double P = pressure();
    if (presArg != -1.0) {
        P = presArg;
    }
    double dAdP;
    switch (m_form_A_Debye) {
    case A_DEBYE_CONST:
        dAdP = 0.0;
        break;
    case A_DEBYE_WATER:
        dAdP = m_waterProps->ADebye(T, P, 3);
        break;
    default:
        throw CanteraError("DebyeHuckel::dA_DebyedP_TP", "shouldn't be here");
    }
    return dAdP;
}

// ---------- Other Property Functions

double DebyeHuckel::AionicRadius(int k) const
{
    return m_Aionic[k];
}

// ------------ Private and Restricted Functions ------------------

bool DebyeHuckel::addSpecies(shared_ptr<Species> spec)
{
    bool added = MolalityVPSSTP::addSpecies(spec);
    if (added) {
        m_electrolyteSpeciesType.push_back(cEST_polarNeutral);
        m_speciesSize.push_back(0.0);
        m_Aionic.push_back(0.0);
        m_lnActCoeffMolal.push_back(0.0);
        m_dlnActCoeffMolaldT.push_back(0.0);
        m_d2lnActCoeffMolaldT2.push_back(0.0);
        m_dlnActCoeffMolaldP.push_back(0.0);
        m_B_Dot.push_back(0.0);
        m_tmpV.push_back(0.0);
    }
    return added;
}

double DebyeHuckel::_nonpolarActCoeff(double IionicMolality)
{
     // These are coefficients to describe the increase in activity coeff for
     // non-polar molecules due to the electrolyte becoming stronger (the so-
     // called salt-out effect)
    const static double npActCoeff[] = {0.1127, -0.01049, 1.545E-3};
    double I2 = IionicMolality * IionicMolality;
    double l10actCoeff =
        npActCoeff[0] * IionicMolality +
        npActCoeff[1] * I2 +
        npActCoeff[2] * I2 * IionicMolality;
    return pow(10.0 , l10actCoeff);
}

double DebyeHuckel::_osmoticCoeffHelgesonFixedForm() const
{
    const double a0 = 1.454;
    const double b0 = 0.02236;
    const double c0 = 9.380E-3;
    const double d0 = -5.362E-4;
    double Is = m_IionicMolalityStoich;
    if (Is <= 0.0) {
        return 0.0;
    }
    double Is2 = Is * Is;
    double bhat = 1.0 + a0 * sqrt(Is);
    double func = bhat - 2.0 * log(bhat) - 1.0/bhat;
    double v1 = m_A_Debye / (a0 * a0 * a0 * Is) * func;
    double oc = 1.0 - v1 + b0 * Is / 2.0 + 2.0 * c0 * Is2 / 3.0
                + 3.0 * d0 * Is2 * Is / 4.0;
    return oc;
}

double DebyeHuckel::_lnactivityWaterHelgesonFixedForm() const
{
    // Update the internally stored vector of molalities
    calcMolalities();
    double oc = _osmoticCoeffHelgesonFixedForm();
    double sum = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        if (k != m_indexSolvent) {
            sum += std::max(m_molalities[k], 0.0);
        }
    }
    if (sum > 2.0 * m_maxIionicStrength) {
        sum = 2.0 * m_maxIionicStrength;
    };
    return - m_Mnaught * sum * oc;
}

void DebyeHuckel::s_update_lnMolalityActCoeff() const
{
    double z_k, zs_k1, zs_k2;

    // Update the internally stored vector of molalities
    calcMolalities();

    // Calculate the apparent (real) ionic strength.
    //
    // Note this is not the stoichiometric ionic strengh, where reactions of
    // ions forming neutral salts are ignorred in calculating the ionic
    // strength.
    m_IionicMolality = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        z_k = m_speciesCharge[k];
        m_IionicMolality += m_molalities[k] * z_k * z_k;
    }
    m_IionicMolality /= 2.0;
    m_IionicMolality = std::min(m_IionicMolality, m_maxIionicStrength);

    // Calculate the stoichiometric ionic charge
    m_IionicMolalityStoich = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        z_k = m_speciesCharge[k];
        zs_k1 = m_speciesCharge_Stoich[k];
        if (z_k == zs_k1) {
            m_IionicMolalityStoich += m_molalities[k] * z_k * z_k;
        } else {
            zs_k2 = z_k - zs_k1;
            m_IionicMolalityStoich += m_molalities[k] * (zs_k1 * zs_k1 + zs_k2 * zs_k2);
        }
    }
    m_IionicMolalityStoich /= 2.0;
    m_IionicMolalityStoich = std::min(m_IionicMolalityStoich, m_maxIionicStrength);

    // Possibly update the stored value of the Debye-Huckel parameter A_Debye
    // This parameter appears on the top of the activity coefficient expression.
    // It depends on T (and P), as it depends explicitly on the temperature.
    // Also, the dielectric constant is usually a fairly strong function of T,
    // also.
    m_A_Debye = A_Debye_TP();

    // Calculate a safe value for the mole fraction of the solvent
    double xmolSolvent = moleFraction(m_indexSolvent);
    xmolSolvent = std::max(8.689E-3, xmolSolvent);

    int est;
    double ac_nonPolar = 1.0;
    double numTmp = m_A_Debye * sqrt(m_IionicMolality);
    double denomTmp = m_B_Debye * sqrt(m_IionicMolality);
    double coeff;
    double lnActivitySolvent = 0.0;
    double tmp;
    double tmpLn;
    double y, yp1, sigma;
    switch (m_formDH) {
    case DHFORM_DILUTE_LIMIT:
        for (size_t k = 0; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            m_lnActCoeffMolal[k] = - z_k * z_k * numTmp;
        }
        lnActivitySolvent =
            (xmolSolvent - 1.0)/xmolSolvent +
            2.0 / 3.0 * m_A_Debye * m_Mnaught *
            m_IionicMolality * sqrt(m_IionicMolality);
        break;

    case DHFORM_BDOT_AK:
        ac_nonPolar = _nonpolarActCoeff(m_IionicMolality);
        for (size_t k = 0; k < m_kk; k++) {
            est = m_electrolyteSpeciesType[k];
            if (est == cEST_nonpolarNeutral) {
                m_lnActCoeffMolal[k] = log(ac_nonPolar);
            } else {
                z_k = m_speciesCharge[k];
                m_lnActCoeffMolal[k] =
                    - z_k * z_k * numTmp / (1.0 + denomTmp * m_Aionic[k])
                    + log(10.0) * m_B_Dot[k] * m_IionicMolality;
            }
        }

        lnActivitySolvent = (xmolSolvent - 1.0)/xmolSolvent;
        coeff = 2.0 / 3.0 * m_A_Debye * m_Mnaught
                * sqrt(m_IionicMolality);
        tmp = 0.0;
        if (denomTmp > 0.0) {
            for (size_t k = 0; k < m_kk; k++) {
                if (k != m_indexSolvent || m_Aionic[k] != 0.0) {
                    y = denomTmp * m_Aionic[k];
                    yp1 = y + 1.0;
                    sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 - 2.0*log(yp1));
                    z_k = m_speciesCharge[k];
                    tmp += m_molalities[k] * z_k * z_k * sigma / 2.0;
                }
            }
        }
        lnActivitySolvent += coeff * tmp;
        tmp = 0.0;
        for (size_t k = 0; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            if ((k != m_indexSolvent) && (z_k != 0.0)) {
                tmp += m_B_Dot[k] * m_molalities[k];
            }
        }
        lnActivitySolvent -=
            m_Mnaught * log(10.0) * m_IionicMolality * tmp / 2.0;

        // Special section to implement the Helgeson fixed form for the water
        // brine activity coefficient.
        if (m_useHelgesonFixedForm) {
            lnActivitySolvent = _lnactivityWaterHelgesonFixedForm();
        }
        break;

    case DHFORM_BDOT_ACOMMON:
        denomTmp *= m_Aionic[0];
        for (size_t k = 0; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            m_lnActCoeffMolal[k] =
                - z_k * z_k * numTmp / (1.0 + denomTmp)
                + log(10.0) * m_B_Dot[k] * m_IionicMolality;
        }
        if (denomTmp > 0.0) {
            y = denomTmp;
            yp1 = y + 1.0;
            sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 - 2.0*log(yp1));
        } else {
            sigma = 0.0;
        }
        lnActivitySolvent =
            (xmolSolvent - 1.0)/xmolSolvent +
            2.0 /3.0 * m_A_Debye * m_Mnaught *
            m_IionicMolality * sqrt(m_IionicMolality) * sigma;
        tmp = 0.0;
        for (size_t k = 0; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            if ((k != m_indexSolvent) && (z_k != 0.0)) {
                tmp += m_B_Dot[k] * m_molalities[k];
            }
        }
        lnActivitySolvent -=
            m_Mnaught * log(10.0) * m_IionicMolality * tmp / 2.0;
        break;

    case DHFORM_BETAIJ:
        denomTmp = m_B_Debye * m_Aionic[0];
        denomTmp *= sqrt(m_IionicMolality);
        lnActivitySolvent =
            (xmolSolvent - 1.0)/xmolSolvent;

        for (size_t k = 0; k < m_kk; k++) {
            if (k != m_indexSolvent) {
                z_k = m_speciesCharge[k];
                m_lnActCoeffMolal[k] =
                    - z_k * z_k * numTmp / (1.0 + denomTmp);
                for (size_t j = 0; j < m_kk; j++) {
                    double beta = m_Beta_ij.value(k, j);
                    m_lnActCoeffMolal[k] += 2.0 * m_molalities[j] * beta;
                }
            }
        }
        if (denomTmp > 0.0) {
            y = denomTmp;
            yp1 = y + 1.0;
            sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 -2.0*log(yp1));
        } else {
            sigma = 0.0;
        }
        lnActivitySolvent =
            (xmolSolvent - 1.0)/xmolSolvent +
            2.0 /3.0 * m_A_Debye * m_Mnaught *
            m_IionicMolality * sqrt(m_IionicMolality) * sigma;
        tmp = 0.0;
        for (size_t k = 0; k < m_kk; k++) {
            for (size_t j = 0; j < m_kk; j++) {
                tmp +=
                    m_Beta_ij.value(k, j) * m_molalities[k] * m_molalities[j];
            }
        }
        lnActivitySolvent -= m_Mnaught * tmp;
        break;

    case DHFORM_PITZER_BETAIJ:
        denomTmp = m_B_Debye * sqrt(m_IionicMolality);
        denomTmp *= m_Aionic[0];
        numTmp = m_A_Debye * sqrt(m_IionicMolality);
        tmpLn = log(1.0 + denomTmp);
        for (size_t k = 0; k < m_kk; k++) {
            if (k != m_indexSolvent) {
                z_k = m_speciesCharge[k];
                m_lnActCoeffMolal[k] =
                    - z_k * z_k * numTmp / 3.0 / (1.0 + denomTmp);
                m_lnActCoeffMolal[k] +=
                    - 2.0 * z_k * z_k * m_A_Debye * tmpLn /
                    (3.0 * m_B_Debye * m_Aionic[0]);
                for (size_t j = 0; j < m_kk; j++) {
                    m_lnActCoeffMolal[k] += 2.0 * m_molalities[j] *
                                            m_Beta_ij.value(k, j);
                }
            }
        }
        sigma = 1.0 / (1.0 + denomTmp);
        lnActivitySolvent =
            (xmolSolvent - 1.0)/xmolSolvent +
            2.0 /3.0 * m_A_Debye * m_Mnaught *
            m_IionicMolality * sqrt(m_IionicMolality) * sigma;
        tmp = 0.0;
        for (size_t k = 0; k < m_kk; k++) {
            for (size_t j = 0; j < m_kk; j++) {
                tmp +=
                    m_Beta_ij.value(k, j) * m_molalities[k] * m_molalities[j];
            }
        }
        lnActivitySolvent -= m_Mnaught * tmp;
        break;

    default:
        throw CanteraError("DebyeHuckel::s_update_lnMolalityActCoeff", "ERROR");
    }

    // Above, we calculated the ln(activitySolvent). Translate that into the
    // molar-based activity coefficient by dividing by the solvent mole
    // fraction. Solvents are not on the molality scale.
    xmolSolvent = moleFraction(m_indexSolvent);
    m_lnActCoeffMolal[m_indexSolvent] =
        lnActivitySolvent - log(xmolSolvent);
}

void DebyeHuckel::s_update_dlnMolalityActCoeff_dT() const
{
    double z_k, coeff, tmp, y, yp1, sigma, tmpLn;
    // First we store dAdT explicitly here
    double dAdT = dA_DebyedT_TP();
    if (dAdT == 0.0) {
        for (size_t k = 0; k < m_kk; k++) {
            m_dlnActCoeffMolaldT[k] = 0.0;
        }
        return;
    }

    // Calculate a safe value for the mole fraction of the solvent
    double xmolSolvent = moleFraction(m_indexSolvent);
    xmolSolvent = std::max(8.689E-3, xmolSolvent);
    double sqrtI = sqrt(m_IionicMolality);
    double numdAdTTmp = dAdT * sqrtI;
    double denomTmp = m_B_Debye * sqrtI;
    double d_lnActivitySolvent_dT = 0;

    switch (m_formDH) {
    case DHFORM_DILUTE_LIMIT:
        for (size_t k = 1; k < m_kk; k++) {
            m_dlnActCoeffMolaldT[k] =
                m_lnActCoeffMolal[k] * dAdT / m_A_Debye;
        }
        d_lnActivitySolvent_dT = 2.0 / 3.0 * dAdT * m_Mnaught *
                                 m_IionicMolality * sqrt(m_IionicMolality);
        m_dlnActCoeffMolaldT[m_indexSolvent] = d_lnActivitySolvent_dT;
        break;

    case DHFORM_BDOT_AK:
        for (size_t k = 0; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            m_dlnActCoeffMolaldT[k] =
                - z_k * z_k * numdAdTTmp / (1.0 + denomTmp * m_Aionic[k]);
        }

        m_dlnActCoeffMolaldT[m_indexSolvent] = 0.0;
        coeff = 2.0 / 3.0 * dAdT * m_Mnaught * sqrtI;
        tmp = 0.0;
        if (denomTmp > 0.0) {
            for (size_t k = 0; k < m_kk; k++) {
                y = denomTmp * m_Aionic[k];
                yp1 = y + 1.0;
                sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 - 2.0*log(yp1));
                z_k = m_speciesCharge[k];
                tmp += m_molalities[k] * z_k * z_k * sigma / 2.0;
            }
        }
        m_dlnActCoeffMolaldT[m_indexSolvent] += coeff * tmp;
        break;

    case DHFORM_BDOT_ACOMMON:
        denomTmp *= m_Aionic[0];
        for (size_t k = 0; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            m_dlnActCoeffMolaldT[k] =
                - z_k * z_k * numdAdTTmp / (1.0 + denomTmp);
        }
        if (denomTmp > 0.0) {
            y = denomTmp;
            yp1 = y + 1.0;
            sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 - 2.0*log(yp1));
        } else {
            sigma = 0.0;
        }
        m_dlnActCoeffMolaldT[m_indexSolvent] =
            2.0 /3.0 * dAdT * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    case DHFORM_BETAIJ:
        denomTmp *= m_Aionic[0];
        for (size_t k = 0; k < m_kk; k++) {
            if (k != m_indexSolvent) {
                z_k = m_speciesCharge[k];
                m_dlnActCoeffMolaldT[k] =
                    - z_k * z_k * numdAdTTmp / (1.0 + denomTmp);
            }
        }
        if (denomTmp > 0.0) {
            y = denomTmp;
            yp1 = y + 1.0;
            sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 - 2.0*log(yp1));
        } else {
            sigma = 0.0;
        }
        m_dlnActCoeffMolaldT[m_indexSolvent] =
            2.0 /3.0 * dAdT * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    case DHFORM_PITZER_BETAIJ:
        denomTmp *= m_Aionic[0];
        tmpLn = log(1.0 + denomTmp);
        for (size_t k = 0; k < m_kk; k++) {
            if (k != m_indexSolvent) {
                z_k = m_speciesCharge[k];
                m_dlnActCoeffMolaldT[k] =
                    - z_k * z_k * numdAdTTmp / (1.0 + denomTmp)
                    - 2.0 * z_k * z_k * dAdT * tmpLn
                    / (m_B_Debye * m_Aionic[0]);
                m_dlnActCoeffMolaldT[k] /= 3.0;
            }
        }

        sigma = 1.0 / (1.0 + denomTmp);
        m_dlnActCoeffMolaldT[m_indexSolvent] =
            2.0 /3.0 * dAdT * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    default:
        throw CanteraError("DebyeHuckel::s_update_dlnMolalityActCoeff_dT",
                           "ERROR");
    }
}

void DebyeHuckel::s_update_d2lnMolalityActCoeff_dT2() const
{
    double z_k, coeff, tmp, y, yp1, sigma, tmpLn;
    double dAdT = dA_DebyedT_TP();
    double d2AdT2 = d2A_DebyedT2_TP();
    if (d2AdT2 == 0.0 && dAdT == 0.0) {
        for (size_t k = 0; k < m_kk; k++) {
            m_d2lnActCoeffMolaldT2[k] = 0.0;
        }
        return;
    }

    // Calculate a safe value for the mole fraction of the solvent
    double xmolSolvent = moleFraction(m_indexSolvent);
    xmolSolvent = std::max(8.689E-3, xmolSolvent);
    double sqrtI = sqrt(m_IionicMolality);
    double numd2AdT2Tmp = d2AdT2 * sqrtI;
    double denomTmp = m_B_Debye * sqrtI;

    switch (m_formDH) {
    case DHFORM_DILUTE_LIMIT:
        for (size_t k = 0; k < m_kk; k++) {
            m_d2lnActCoeffMolaldT2[k] =
                m_lnActCoeffMolal[k] * d2AdT2 / m_A_Debye;
        }
        break;

    case DHFORM_BDOT_AK:
        for (size_t k = 0; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            m_d2lnActCoeffMolaldT2[k] =
                - z_k * z_k * numd2AdT2Tmp / (1.0 + denomTmp * m_Aionic[k]);
        }

        m_d2lnActCoeffMolaldT2[m_indexSolvent] = 0.0;
        coeff = 2.0 / 3.0 * d2AdT2 * m_Mnaught * sqrtI;
        tmp = 0.0;
        if (denomTmp > 0.0) {
            for (size_t k = 0; k < m_kk; k++) {
                y = denomTmp * m_Aionic[k];
                yp1 = y + 1.0;
                sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 - 2.0*log(yp1));
                z_k = m_speciesCharge[k];
                tmp += m_molalities[k] * z_k * z_k * sigma / 2.0;
            }
        }
        m_d2lnActCoeffMolaldT2[m_indexSolvent] += coeff * tmp;
        break;

    case DHFORM_BDOT_ACOMMON:
        denomTmp *= m_Aionic[0];
        for (size_t k = 0; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            m_d2lnActCoeffMolaldT2[k] =
                - z_k * z_k * numd2AdT2Tmp / (1.0 + denomTmp);
        }
        if (denomTmp > 0.0) {
            y = denomTmp;
            yp1 = y + 1.0;
            sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 - 2.0*log(yp1));
        } else {
            sigma = 0.0;
        }
        m_d2lnActCoeffMolaldT2[m_indexSolvent] =
            2.0 /3.0 * d2AdT2 * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    case DHFORM_BETAIJ:
        denomTmp *= m_Aionic[0];
        for (size_t k = 0; k < m_kk; k++) {
            if (k != m_indexSolvent) {
                z_k = m_speciesCharge[k];
                m_d2lnActCoeffMolaldT2[k] =
                    - z_k * z_k * numd2AdT2Tmp / (1.0 + denomTmp);
            }
        }
        if (denomTmp > 0.0) {
            y = denomTmp;
            yp1 = y + 1.0;
            sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 -2.0*log(yp1));
        } else {
            sigma = 0.0;
        }
        m_d2lnActCoeffMolaldT2[m_indexSolvent] =
            2.0 /3.0 * d2AdT2 * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    case DHFORM_PITZER_BETAIJ:
        denomTmp *= m_Aionic[0];
        tmpLn = log(1.0 + denomTmp);
        for (size_t k = 0; k < m_kk; k++) {
            if (k != m_indexSolvent) {
                z_k = m_speciesCharge[k];
                m_d2lnActCoeffMolaldT2[k] =
                    - z_k * z_k * numd2AdT2Tmp / (1.0 + denomTmp)
                    - 2.0 * z_k * z_k * d2AdT2 * tmpLn
                    / (m_B_Debye * m_Aionic[0]);
                m_d2lnActCoeffMolaldT2[k] /= 3.0;
            }
        }

        sigma = 1.0 / (1.0 + denomTmp);
        m_d2lnActCoeffMolaldT2[m_indexSolvent] =
            2.0 /3.0 * d2AdT2 * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    default:
        throw CanteraError("DebyeHuckel::s_update_d2lnMolalityActCoeff_dT2",
                           "ERROR");
    }
}

void DebyeHuckel::s_update_dlnMolalityActCoeff_dP() const
{
    double z_k, coeff, tmp, y, yp1, sigma, tmpLn;
    int est;
    double dAdP = dA_DebyedP_TP();
    if (dAdP == 0.0) {
        for (size_t k = 0; k < m_kk; k++) {
            m_dlnActCoeffMolaldP[k] = 0.0;
        }
        return;
    }

    // Calculate a safe value for the mole fraction of the solvent
    double xmolSolvent = moleFraction(m_indexSolvent);
    xmolSolvent = std::max(8.689E-3, xmolSolvent);
    double sqrtI = sqrt(m_IionicMolality);
    double numdAdPTmp = dAdP * sqrtI;
    double denomTmp = m_B_Debye * sqrtI;

    switch (m_formDH) {
    case DHFORM_DILUTE_LIMIT:
        for (size_t k = 0; k < m_kk; k++) {
            m_dlnActCoeffMolaldP[k] =
                m_lnActCoeffMolal[k] * dAdP / m_A_Debye;
        }
        break;

    case DHFORM_BDOT_AK:
        for (size_t k = 0; k < m_kk; k++) {
            est = m_electrolyteSpeciesType[k];
            if (est == cEST_nonpolarNeutral) {
                m_lnActCoeffMolal[k] = 0.0;
            } else {
                z_k = m_speciesCharge[k];
                m_dlnActCoeffMolaldP[k] =
                    - z_k * z_k * numdAdPTmp / (1.0 + denomTmp * m_Aionic[k]);
            }
        }

        m_dlnActCoeffMolaldP[m_indexSolvent] = 0.0;
        coeff = 2.0 / 3.0 * dAdP * m_Mnaught * sqrtI;
        tmp = 0.0;
        if (denomTmp > 0.0) {
            for (size_t k = 0; k < m_kk; k++) {
                y = denomTmp * m_Aionic[k];
                yp1 = y + 1.0;
                sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 - 2.0*log(yp1));
                z_k = m_speciesCharge[k];
                tmp += m_molalities[k] * z_k * z_k * sigma / 2.0;
            }
        }
        m_dlnActCoeffMolaldP[m_indexSolvent] += coeff * tmp;
        break;

    case DHFORM_BDOT_ACOMMON:
        denomTmp *= m_Aionic[0];
        for (size_t k = 0; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            m_dlnActCoeffMolaldP[k] =
                - z_k * z_k * numdAdPTmp / (1.0 + denomTmp);
        }
        if (denomTmp > 0.0) {
            y = denomTmp;
            yp1 = y + 1.0;
            sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 - 2.0*log(yp1));
        } else {
            sigma = 0.0;
        }
        m_dlnActCoeffMolaldP[m_indexSolvent] =
            2.0 /3.0 * dAdP * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    case DHFORM_BETAIJ:
        denomTmp *= m_Aionic[0];
        for (size_t k = 0; k < m_kk; k++) {
            if (k != m_indexSolvent) {
                z_k = m_speciesCharge[k];
                m_dlnActCoeffMolaldP[k] =
                    - z_k * z_k * numdAdPTmp / (1.0 + denomTmp);
            }
        }
        if (denomTmp > 0.0) {
            y = denomTmp;
            yp1 = y + 1.0;
            sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 - 2.0*log(yp1));
        } else {
            sigma = 0.0;
        }
        m_dlnActCoeffMolaldP[m_indexSolvent] =
            2.0 /3.0 * dAdP * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    case DHFORM_PITZER_BETAIJ:
        denomTmp *= m_Aionic[0];
        tmpLn = log(1.0 + denomTmp);
        for (size_t k = 0; k < m_kk; k++) {
            if (k != m_indexSolvent) {
                z_k = m_speciesCharge[k];
                m_dlnActCoeffMolaldP[k] =
                    - z_k * z_k * numdAdPTmp / (1.0 + denomTmp)
                    - 2.0 * z_k * z_k * dAdP * tmpLn
                    / (m_B_Debye * m_Aionic[0]);
                m_dlnActCoeffMolaldP[k] /= 3.0;
            }
        }

        sigma = 1.0 / (1.0 + denomTmp);
        m_dlnActCoeffMolaldP[m_indexSolvent] =
            2.0 /3.0 * dAdP * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    default:
        throw CanteraError("DebyeHuckel::s_update_dlnMolalityActCoeff_dP",
                           "ERROR");
    }
}

}
