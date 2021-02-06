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
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/DebyeHuckel.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/PDSS_Water.h"
#include "cantera/thermo/PDSS_ConstVol.h"
#include "cantera/thermo/electrolytes.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"

#include <cstdio>

using namespace std;

namespace Cantera
{

DebyeHuckel::DebyeHuckel() :
    m_formDH(DHFORM_DILUTE_LIMIT),
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

DebyeHuckel::~DebyeHuckel()
{
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
    Phase::assignDensity(dd);
}

void DebyeHuckel::setDensity(doublereal rho)
{
    warn_deprecated("DebyeHuckel::setDensity",
        "Overloaded function to be removed after Cantera 2.5. "
        "Error will be thrown by Phase::setDensity instead");
    double dens = density();
    if (rho != dens) {
        throw CanteraError("DebyeHuckel::setDensity",
                           "Density is not an independent variable");
    }
}

void DebyeHuckel::setMolarDensity(const doublereal conc)
{
    warn_deprecated("DebyeHuckel::setMolarDensity",
        "Overloaded function to be removed after Cantera 2.5. "
        "Error will be thrown by Phase::setMolarDensity instead");
    double concI = molarDensity();
    if (conc != concI) {
        throw CanteraError("DebyeHuckel::setMolarDensity",
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
    double mvSolvent = providePDSS(0)->molarVolume();
    return 1.0 / mvSolvent;
}

void DebyeHuckel::getActivities(doublereal* ac) const
{
    _updateStandardStateThermo();

    // Update the molality array, m_molalities(). This requires an update due to
    // mole fractions
    s_update_lnMolalityActCoeff();
    for (size_t k = 1; k < m_kk; k++) {
        ac[k] = m_molalities[k] * exp(m_lnActCoeffMolal[k]);
    }
    double xmolSolvent = moleFraction(0);
    ac[0] = exp(m_lnActCoeffMolal[0]) * xmolSolvent;
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
    double xmolSolvent = moleFraction(0);
    for (size_t k = 1; k < m_kk; k++) {
        xx = std::max(m_molalities[k], SmallNumber);
        mu[k] += RT() * (log(xx) + m_lnActCoeffMolal[k]);
    }
    xx = std::max(xmolSolvent, SmallNumber);
    mu[0] += RT() * (log(xx) + m_lnActCoeffMolal[0]);
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
    for (size_t k = 1; k < m_kk; k++) {
        mm = std::max(SmallNumber, m_molalities[k]);
        sbar[k] -= GasConstant * (log(mm) + m_lnActCoeffMolal[k]);
    }
    double xmolSolvent = moleFraction(0);
    mm = std::max(SmallNumber, xmolSolvent);
    sbar[0] -= GasConstant *(log(mm) + m_lnActCoeffMolal[0]);

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
    if (caseInsensitiveEquals(estString, "solvent")) {
        return cEST_solvent;
    } else if (estString == "charged-species"
               || caseInsensitiveEquals(estString, "chargedspecies")) {
        return cEST_chargedSpecies;
    } else if (estString == "weak-acid-associated"
               || caseInsensitiveEquals(estString, "weakacidassociated")) {
        return cEST_weakAcidAssociated;
    } else if (estString == "strong-acid-associated"
               || caseInsensitiveEquals(estString, "strongacidassociated")) {
        return cEST_strongAcidAssociated;
    } else if (estString == "polar-neutral"
               || caseInsensitiveEquals(estString, "polarneutral")) {
        return cEST_polarNeutral;
    } else if (estString == "nonpolar-neutral"
               || caseInsensitiveEquals(estString, "nonpolarneutral")) {
        return cEST_nonpolarNeutral;
    } else {
        throw CanteraError("interp_est (DebyeHuckel)",
            "Invalid electrolyte species type '{}'", estString);
    }
}

void DebyeHuckel::setDebyeHuckelModel(const std::string& model) {
    if (model == ""
        || model == "dilute-limit"
        || caseInsensitiveEquals(model, "Dilute_limit")) {
        m_formDH = DHFORM_DILUTE_LIMIT;
    } else if (model == "B-dot-with-variable-a"
               || caseInsensitiveEquals(model, "Bdot_with_variable_a")) {
        m_formDH = DHFORM_BDOT_AK;
    } else if (model == "B-dot-with-common-a"
               || caseInsensitiveEquals(model, "Bdot_with_common_a")) {
        m_formDH = DHFORM_BDOT_ACOMMON;
    } else if (caseInsensitiveEquals(model, "beta_ij")) {
        m_formDH = DHFORM_BETAIJ;
        m_Beta_ij.resize(m_kk, m_kk, 0.0);
    } else if (model == "Pitzer-with-beta_ij"
               || caseInsensitiveEquals(model, "Pitzer_with_Beta_ij")) {
        m_formDH = DHFORM_PITZER_BETAIJ;
        m_Beta_ij.resize(m_kk, m_kk, 0.0);
    } else {
        throw CanteraError("DebyeHuckel::setDebyeHuckelModel",
                           "Unknown model '{}'", model);
    }
}

void DebyeHuckel::setA_Debye(double A)
{
    if (A < 0) {
        m_form_A_Debye = A_DEBYE_WATER;
    } else {
        m_form_A_Debye = A_DEBYE_CONST;
        m_A_Debye = A;
    }
}

void DebyeHuckel::setB_dot(double bdot)
{
    if (m_formDH == DHFORM_BETAIJ || m_formDH == DHFORM_DILUTE_LIMIT ||
            m_formDH == DHFORM_PITZER_BETAIJ) {
        throw CanteraError("DebyeHuckel::setB_dot",
                           "B_dot entry in the wrong DH form");
    }
    // Set B_dot parameters for charged species
    for (size_t k = 0; k < nSpecies(); k++) {
        if (fabs(charge(k)) > 0.0001) {
            m_B_Dot[k] = bdot;
        } else {
            m_B_Dot[k] = 0.0;
        }
    }
}

void DebyeHuckel::setDefaultIonicRadius(double value)
{
    for (size_t k = 0; k < m_kk; k++) {
        if (std::isnan(m_Aionic[k])) {
            m_Aionic[k] = value;
        }
    }
}

void DebyeHuckel::setBeta(const std::string& sp1, const std::string& sp2,
                          double value)
{
    size_t k1 = speciesIndex(sp1);
    if (k1 == npos) {
        throw CanteraError("DebyeHuckel::setBeta", "Species '{}' not found", sp1);
    }
    size_t k2 = speciesIndex(sp2);
    if (k2 == npos) {
        throw CanteraError("DebyeHuckel::setBeta", "Species '{}' not found", sp2);
    }
    m_Beta_ij(k1, k2) = value;
    m_Beta_ij(k2, k1) = value;
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

    // Determine the form of the Debye-Huckel model, m_formDH. We will use
    // this information to size arrays below. If there is no XML node named
    // "activityCoefficients", assume that we are doing the extreme dilute
    // limit assumption
    if (thermoNode.hasChild("activityCoefficients")) {
        setDebyeHuckelModel(thermoNode.child("activityCoefficients")["model"]);
    } else {
        setDebyeHuckelModel("Dilute_limit");
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
                if (caseInsensitiveEquals(modelString, "water")) {
                    setA_Debye(-1);
                } else {
                    throw CanteraError("DebyeHuckel::initThermoXML",
                                       "A_Debye Model \"" + modelString +
                                       "\" is not known");
                }
            } else {
                setA_Debye(getFloat(acNode, "A_Debye"));
            }
        }

        // Look for parameters for B_Debye
        if (acNode.hasChild("B_Debye")) {
            setB_Debye(getFloat(acNode, "B_Debye"));
        }

        // Look for parameters for B_dot
        if (acNode.hasChild("B_dot")) {
            setB_dot(getFloat(acNode, "B_dot"));
        }

        // Look for Parameters for the Maximum Ionic Strength
        if (acNode.hasChild("maxIonicStrength")) {
            setMaxIonicStrength(getFloat(acNode, "maxIonicStrength"));
        }

        // Look for Helgeson Parameters
        useHelgesonFixedForm(acNode.hasChild("UseHelgesonFixedForm"));

        // Look for parameters for the Ionic radius
        if (acNode.hasChild("ionicRadius")) {
            XML_Node& irNode = acNode.child("ionicRadius");

            double Afactor = 1.0;
            if (irNode.hasAttrib("units")) {
                std::string Aunits = irNode.attrib("units");
                Afactor = toSI(Aunits);
            }

            if (irNode.hasAttrib("default")) {
                setDefaultIonicRadius(Afactor * fpValue(irNode.attrib("default")));
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
                    size_t k = speciesIndex(b.first);
                    if (k != npos) {
                        m_Aionic[k] = fpValue(b.second) * Afactor;
                    }
                }
            }
        }

        // Get the matrix of coefficients for the Beta binary interaction
        // parameters. We assume here that this matrix is symmetric, so that we
        // only have to input 1/2 of the values.
        if (acNode.hasChild("DHBetaMatrix")) {
            if (m_formDH == DHFORM_BETAIJ ||
                    m_formDH == DHFORM_PITZER_BETAIJ) {
                XML_Node& irNode = acNode.child("DHBetaMatrix");
                const vector<string>& sn = speciesNames();
                getMatrixValues(irNode, sn, sn, m_Beta_ij, true, true);
            } else {
                throw CanteraError("DebyeHuckel::initThermoXML:",
                                   "DHBetaMatrix found for wrong type");
            }
        }

        // Override stoichiometric Ionic Strength based on the phase definition
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

    // Override electrolyte species type based on the phase definition
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

void DebyeHuckel::initThermo()
{
    MolalityVPSSTP::initThermo();
    if (m_input.hasKey("activity-data")) {
        auto& node = m_input["activity-data"].as<AnyMap>();
        setDebyeHuckelModel(node["model"].asString());
        if (node.hasKey("A_Debye")) {
            if (node["A_Debye"] == "variable") {
                setA_Debye(-1);
            } else {
                setA_Debye(node.convert("A_Debye", "kg^0.5/gmol^0.5"));
            }
        }
        if (node.hasKey("B_Debye")) {
            setB_Debye(node.convert("B_Debye", "kg^0.5/gmol^0.5/m"));
        }
        if (node.hasKey("max-ionic-strength")) {
            setMaxIonicStrength(node["max-ionic-strength"].asDouble());
        }
        if (node.hasKey("use-Helgeson-fixed-form")) {
            useHelgesonFixedForm(node["use-Helgeson-fixed-form"].asBool());
        }
        if (node.hasKey("default-ionic-radius")) {
            setDefaultIonicRadius(node.convert("default-ionic-radius", "m"));
        }
        if (node.hasKey("B-dot")) {
            setB_dot(node["B-dot"].asDouble());
        }
        if (node.hasKey("beta")) {
            for (auto& item : node["beta"].asVector<AnyMap>()) {
                auto& species = item["species"].asVector<string>(2);
                setBeta(species[0], species[1], item["beta"].asDouble());
            }
        }
    }

    // Solvent
    m_waterSS = dynamic_cast<PDSS_Water*>(providePDSS(0));
    if (m_waterSS) {
        // Initialize the water property calculator. It will share the internal
        // eos water calculator.
        if (m_form_A_Debye == A_DEBYE_WATER) {
            m_waterProps.reset(new WaterProps(m_waterSS));
        }
    } else if (dynamic_cast<PDSS_ConstVol*>(providePDSS(0)) == 0) {
        throw CanteraError("DebyeHuckel::initThermo", "Solvent standard state"
            " model must be WaterIAPWS or constant_incompressible.");
    }

    // Solutes
    for (size_t k = 1; k < nSpecies(); k++) {
        if (dynamic_cast<PDSS_ConstVol*>(providePDSS(k)) == 0) {
            throw CanteraError("DebyeHuckel::initThermo", "Solute standard"
                " state model must be constant_incompressible.");
        }
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
        m_lnActCoeffMolal.push_back(0.0);
        m_dlnActCoeffMolaldT.push_back(0.0);
        m_d2lnActCoeffMolaldT2.push_back(0.0);
        m_dlnActCoeffMolaldP.push_back(0.0);
        m_B_Dot.push_back(0.0);
        m_tmpV.push_back(0.0);

        // NAN will be replaced with default value
        double Aionic = NAN;

        // Guess electrolyte species type based on charge properties
        int est = cEST_nonpolarNeutral;
        double stoichCharge = spec->charge;
        if (fabs(spec->charge) > 0.0001) {
            est = cEST_chargedSpecies;
        }

        if (spec->input.hasKey("Debye-Huckel")) {
            auto& dhNode = spec->input["Debye-Huckel"].as<AnyMap>();
            Aionic = dhNode.convert("ionic-radius", "m", NAN);
            if (dhNode.hasKey("weak-acid-charge")) {
                stoichCharge = dhNode["weak-acid-charge"].asDouble();
                if (fabs(stoichCharge - spec->charge) > 0.0001) {
                    est = cEST_weakAcidAssociated;
                }
            }
            // Apply override of the electrolyte species type
            if (dhNode.hasKey("electrolyte-species-type")) {
                est = interp_est(dhNode["electrolyte-species-type"].asString());
            }
        }

        if (m_electrolyteSpeciesType.size() == 0) {
            est = cEST_solvent; // species 0 is the solvent
        }

        m_Aionic.push_back(Aionic);
        m_speciesCharge_Stoich.push_back(stoichCharge);
        m_electrolyteSpeciesType.push_back(est);
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
    for (size_t k = 1; k < m_kk; k++) {
        sum += std::max(m_molalities[k], 0.0);
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
    double xmolSolvent = moleFraction(0);
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
                if (k != 0 || m_Aionic[k] != 0.0) {
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
        for (size_t k = 1; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            if (z_k != 0.0) {
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
        for (size_t k = 1; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            if (z_k != 0.0) {
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

        for (size_t k = 1; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            m_lnActCoeffMolal[k] =
                - z_k * z_k * numTmp / (1.0 + denomTmp);
            for (size_t j = 0; j < m_kk; j++) {
                double beta = m_Beta_ij.value(k, j);
                m_lnActCoeffMolal[k] += 2.0 * m_molalities[j] * beta;
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
        for (size_t k = 1; k < m_kk; k++) {
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
    xmolSolvent = moleFraction(0);
    m_lnActCoeffMolal[0] = lnActivitySolvent - log(xmolSolvent);
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
    double xmolSolvent = moleFraction(0);
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
        m_dlnActCoeffMolaldT[0] = d_lnActivitySolvent_dT;
        break;

    case DHFORM_BDOT_AK:
        for (size_t k = 0; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            m_dlnActCoeffMolaldT[k] =
                - z_k * z_k * numdAdTTmp / (1.0 + denomTmp * m_Aionic[k]);
        }

        m_dlnActCoeffMolaldT[0] = 0.0;
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
        m_dlnActCoeffMolaldT[0] += coeff * tmp;
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
        m_dlnActCoeffMolaldT[0] = 2.0 /3.0 * dAdT * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    case DHFORM_BETAIJ:
        denomTmp *= m_Aionic[0];
        for (size_t k = 1; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            m_dlnActCoeffMolaldT[k] = -z_k*z_k * numdAdTTmp / (1.0 + denomTmp);
        }
        if (denomTmp > 0.0) {
            y = denomTmp;
            yp1 = y + 1.0;
            sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 - 2.0*log(yp1));
        } else {
            sigma = 0.0;
        }
        m_dlnActCoeffMolaldT[0] = 2.0 /3.0 * dAdT * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    case DHFORM_PITZER_BETAIJ:
        denomTmp *= m_Aionic[0];
        tmpLn = log(1.0 + denomTmp);
        for (size_t k = 1; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            m_dlnActCoeffMolaldT[k] =
                - z_k * z_k * numdAdTTmp / (1.0 + denomTmp)
                - 2.0 * z_k * z_k * dAdT * tmpLn / (m_B_Debye * m_Aionic[0]);
            m_dlnActCoeffMolaldT[k] /= 3.0;
        }

        sigma = 1.0 / (1.0 + denomTmp);
        m_dlnActCoeffMolaldT[0] = 2.0 /3.0 * dAdT * m_Mnaught *
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
    double xmolSolvent = moleFraction(0);
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

        m_d2lnActCoeffMolaldT2[0] = 0.0;
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
        m_d2lnActCoeffMolaldT2[0] += coeff * tmp;
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
        m_d2lnActCoeffMolaldT2[0] = 2.0 /3.0 * d2AdT2 * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    case DHFORM_BETAIJ:
        denomTmp *= m_Aionic[0];
        for (size_t k = 1; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            m_d2lnActCoeffMolaldT2[k] = -z_k*z_k * numd2AdT2Tmp / (1.0 + denomTmp);
        }
        if (denomTmp > 0.0) {
            y = denomTmp;
            yp1 = y + 1.0;
            sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 -2.0*log(yp1));
        } else {
            sigma = 0.0;
        }
        m_d2lnActCoeffMolaldT2[0] = 2.0 /3.0 * d2AdT2 * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    case DHFORM_PITZER_BETAIJ:
        denomTmp *= m_Aionic[0];
        tmpLn = log(1.0 + denomTmp);
        for (size_t k = 1; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            m_d2lnActCoeffMolaldT2[k] =
                - z_k * z_k * numd2AdT2Tmp / (1.0 + denomTmp)
                - 2.0 * z_k * z_k * d2AdT2 * tmpLn / (m_B_Debye * m_Aionic[0]);
            m_d2lnActCoeffMolaldT2[k] /= 3.0;
        }

        sigma = 1.0 / (1.0 + denomTmp);
        m_d2lnActCoeffMolaldT2[0] = 2.0 /3.0 * d2AdT2 * m_Mnaught *
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
    double xmolSolvent = moleFraction(0);
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

        m_dlnActCoeffMolaldP[0] = 0.0;
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
        m_dlnActCoeffMolaldP[0] += coeff * tmp;
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
        m_dlnActCoeffMolaldP[0] =
            2.0 /3.0 * dAdP * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    case DHFORM_BETAIJ:
        denomTmp *= m_Aionic[0];
        for (size_t k = 1; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            m_dlnActCoeffMolaldP[k] = - z_k*z_k * numdAdPTmp / (1.0 + denomTmp);
        }
        if (denomTmp > 0.0) {
            y = denomTmp;
            yp1 = y + 1.0;
            sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 - 2.0*log(yp1));
        } else {
            sigma = 0.0;
        }
        m_dlnActCoeffMolaldP[0] = 2.0 /3.0 * dAdP * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    case DHFORM_PITZER_BETAIJ:
        denomTmp *= m_Aionic[0];
        tmpLn = log(1.0 + denomTmp);
        for (size_t k = 1; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            m_dlnActCoeffMolaldP[k] =
                - z_k * z_k * numdAdPTmp / (1.0 + denomTmp)
                - 2.0 * z_k * z_k * dAdP * tmpLn
                / (m_B_Debye * m_Aionic[0]);
            m_dlnActCoeffMolaldP[k] /= 3.0;
        }

        sigma = 1.0 / (1.0 + denomTmp);
        m_dlnActCoeffMolaldP[0] = 2.0 /3.0 * dAdP * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    default:
        throw CanteraError("DebyeHuckel::s_update_dlnMolalityActCoeff_dP",
                           "ERROR");
    }
}

}
