/**
 *  @file MolarityIonicVPSSTP.cpp
 *   Definitions for intermediate ThermoPhase object for phases which
 *   employ excess Gibbs free energy formulations
 *  (see \ref thermoprops
 * and class \link Cantera::MolarityIonicVPSSTP MolarityIonicVPSSTP\endlink).
 *
 * Header file for a derived class of ThermoPhase that handles variable pressure
 * standard state methods for calculating thermodynamic properties that are
 * further based upon expressions for the excess Gibbs free energy expressed as
 * a function of the mole fractions.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/MolarityIonicVPSSTP.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"

#include <cstdio>

using namespace std;

namespace Cantera
{

MolarityIonicVPSSTP::MolarityIonicVPSSTP() :
    PBType_(PBTYPE_PASSTHROUGH),
    numPBSpecies_(m_kk),
    neutralPBindexStart(0)
{
    warn_deprecated("Class MolarityIonicVPSSTP", "To be removed after Cantera 2.4");
}

MolarityIonicVPSSTP::MolarityIonicVPSSTP(const std::string& inputFile,
        const std::string& id_) :
    PBType_(PBTYPE_PASSTHROUGH),
    numPBSpecies_(m_kk),
    neutralPBindexStart(0)
{
    warn_deprecated("Class MolarityIonicVPSSTP", "To be removed after Cantera 2.4");
    initThermoFile(inputFile, id_);
}

MolarityIonicVPSSTP::MolarityIonicVPSSTP(XML_Node& phaseRoot,
        const std::string& id_) :
    PBType_(PBTYPE_PASSTHROUGH),
    numPBSpecies_(m_kk),
    neutralPBindexStart(0)
{
    warn_deprecated("Class MolarityIonicVPSSTP", "To be removed after Cantera 2.4");
    importPhase(phaseRoot, this);
}

// - Activities, Standard States, Activity Concentrations -----------

void MolarityIonicVPSSTP::getLnActivityCoefficients(doublereal* lnac) const
{
    // Update the activity coefficients
    s_update_lnActCoeff();

    // take the exp of the internally stored coefficients.
    for (size_t k = 0; k < m_kk; k++) {
        lnac[k] = lnActCoeff_Scaled_[k];
    }
}

void MolarityIonicVPSSTP::getChemPotentials(doublereal* mu) const
{
    // First get the standard chemical potentials in molar form. This requires
    // updates of standard state as a function of T and P
    getStandardChemPotentials(mu);

    // Update the activity coefficients
    s_update_lnActCoeff();
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(moleFractions_[k], SmallNumber);
        mu[k] += RT() * (log(xx) + lnActCoeff_Scaled_[k]);
    }
}

void MolarityIonicVPSSTP::getPartialMolarEnthalpies(doublereal* hbar) const
{
    // Get the nondimensional standard state enthalpies
    getEnthalpy_RT(hbar);

    // dimensionalize it.
    double T = temperature();
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] *= GasConstant * T;
    }

    // Update the activity coefficients, This also update the internally stored
    // molalities.
    s_update_lnActCoeff();
    s_update_dlnActCoeff_dT();
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] -= GasConstant * T * T * dlnActCoeffdT_Scaled_[k];
    }
}

void MolarityIonicVPSSTP::getPartialMolarCp(doublereal* cpbar) const
{
    // Get the nondimensional standard state entropies
    getCp_R(cpbar);
    double T = temperature();

    // Update the activity coefficients, This also update the internally stored
    // molalities.
    s_update_lnActCoeff();
    s_update_dlnActCoeff_dT();

    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] -= 2 * T * dlnActCoeffdT_Scaled_[k] + T * T * d2lnActCoeffdT2_Scaled_[k];
    }

    // dimensionalize it.
    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] *= GasConstant;
    }
}

void MolarityIonicVPSSTP::getPartialMolarEntropies(doublereal* sbar) const
{
    // Get the nondimensional standard state entropies
    getEntropy_R(sbar);
    double T = temperature();

    // Update the activity coefficients, This also update the internally stored
    // molalities.
    s_update_lnActCoeff();
    s_update_dlnActCoeff_dT();

    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(moleFractions_[k], SmallNumber);
        sbar[k] += - lnActCoeff_Scaled_[k] -log(xx) - T * dlnActCoeffdT_Scaled_[k];
    }

    // dimensionalize it.
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] *= GasConstant;
    }
}

void MolarityIonicVPSSTP::getPartialMolarVolumes(doublereal* vbar) const
{
    // Get the standard state values in m^3 kmol-1
    getStandardVolumes(vbar);
    for (size_t iK = 0; iK < m_kk; iK++) {
        vbar[iK] += 0.0;
    }
}

void MolarityIonicVPSSTP::calcPseudoBinaryMoleFractions() const
{
    switch (PBType_) {
    case PBTYPE_PASSTHROUGH:
        for (size_t k = 0; k < m_kk; k++) {
            PBMoleFractions_[k] = moleFractions_[k];
        }
        break;
    case PBTYPE_SINGLEANION:
    {
        double sumCat = 0.0;
        double sumAnion = 0.0;
        for (size_t k = 0; k < m_kk; k++) {
            moleFractionsTmp_[k] = moleFractions_[k];
        }
        size_t kMax = npos;
        double sumMax = 0.0;
        for (size_t k = 0; k < cationList_.size(); k++) {
            size_t kCat = cationList_[k];
            double chP = m_speciesCharge[kCat];
            if (moleFractions_[kCat] > sumMax) {
                kMax = k;
                sumMax = moleFractions_[kCat];
            }
            sumCat += chP * moleFractions_[kCat];
        }
        size_t ka = anionList_[0];
        sumAnion = moleFractions_[ka] * m_speciesCharge[ka];
        double sum = sumCat - sumAnion;
        if (fabs(sum) > 1.0E-16) {
            moleFractionsTmp_[cationList_[kMax]] -= sum / m_speciesCharge[kMax];
            sum = 0.0;
            for (size_t k = 0; k < cationList_.size(); k++) {
                sum += moleFractionsTmp_[k];
            }
            for (size_t k = 0; k < cationList_.size(); k++) {
                moleFractionsTmp_[k]/= sum;
            }
        }

        for (size_t k = 0; k < cationList_.size(); k++) {
            PBMoleFractions_[k] = moleFractionsTmp_[cationList_[k]];
        }
        for (size_t k = 0; k < passThroughList_.size(); k++) {
            PBMoleFractions_[neutralPBindexStart + k] = moleFractions_[passThroughList_[k]];
        }

        sum = std::max(0.0, PBMoleFractions_[0]);
        for (size_t k = 1; k < numPBSpecies_; k++) {
            sum += PBMoleFractions_[k];
        }
        for (size_t k = 0; k < numPBSpecies_; k++) {
            PBMoleFractions_[k] /= sum;
        }
        break;
    }
    case PBTYPE_SINGLECATION:
        throw CanteraError("eosType", "Unknown type");
    case PBTYPE_MULTICATIONANION:
        throw CanteraError("eosType", "Unknown type");
    default:
        throw CanteraError("eosType", "Unknown type");
    }
}

void MolarityIonicVPSSTP::s_update_lnActCoeff() const
{
    for (size_t k = 0; k < m_kk; k++) {
        lnActCoeff_Scaled_[k] = 0.0;
    }
}

void MolarityIonicVPSSTP::s_update_dlnActCoeff_dT() const
{
}

void MolarityIonicVPSSTP::s_update_dlnActCoeff_dX_() const
{
}

void MolarityIonicVPSSTP::initThermo()
{
    GibbsExcessVPSSTP::initThermo();
    initLengths();

    // Go find the list of cations and anions
    cationList_.clear();
    anionList_.clear();
    passThroughList_.clear();
    for (size_t k = 0; k < m_kk; k++) {
        double ch = m_speciesCharge[k];
        if (ch > 0.0) {
            cationList_.push_back(k);
        } else if (ch < 0.0) {
            anionList_.push_back(k);
        } else {
            passThroughList_.push_back(k);
        }
    }
    numPBSpecies_ = cationList_.size() + anionList_.size() - 1;
    neutralPBindexStart = numPBSpecies_;
    PBType_ = PBTYPE_MULTICATIONANION;
    if (anionList_.size() == 1) {
        PBType_ = PBTYPE_SINGLEANION;
    } else if (cationList_.size() == 1) {
        PBType_ = PBTYPE_SINGLECATION;
    }
    if (anionList_.size() == 0 && cationList_.size() == 0) {
        PBType_ = PBTYPE_PASSTHROUGH;
    }
}

void MolarityIonicVPSSTP::initLengths()
{
    moleFractionsTmp_.resize(m_kk);
}

void MolarityIonicVPSSTP::initThermoXML(XML_Node& phaseNode, const std::string& id)
{
    if ((int) id.size() > 0 && phaseNode.id() != id) {
        throw CanteraError("MolarityIonicVPSSTP::initThermoXML",
                           "phasenode and Id are incompatible");
    }

    // Check on the thermo field. Must have one of:
    //     <thermo model="MolarityIonicVPSS" />
    //     <thermo model="MolarityIonicVPSSTP" />
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError("MolarityIonicVPSSTP::initThermoXML",
                           "no thermo XML node");
    }
    XML_Node& thermoNode = phaseNode.child("thermo");
    if (!caseInsensitiveEquals(thermoNode["model"], "molarityionicvpss")
        && !caseInsensitiveEquals(thermoNode["model"], "molarityionicvpsstp")) {
        throw CanteraError("MolarityIonicVPSSTP::initThermoXML",
                           "Unknown thermo model: " + thermoNode["model"]
                           + " - This object only knows \"MolarityIonicVPSSTP\" ");
    }

    // Go get all of the coefficients and factors in the activityCoefficients
    // XML block
    if (thermoNode.hasChild("activityCoefficients")) {
        XML_Node& acNode = thermoNode.child("activityCoefficients");
        for (size_t i = 0; i < acNode.nChildren(); i++) {
            XML_Node& xmlACChild = acNode.child(i);
            // Process a binary interaction
            if (caseInsensitiveEquals(xmlACChild.name(), "binaryneutralspeciesparameters")) {
                readXMLBinarySpecies(xmlACChild);
            }
        }
    }

    // Go down the chain
    GibbsExcessVPSSTP::initThermoXML(phaseNode, id);
}

void MolarityIonicVPSSTP::readXMLBinarySpecies(XML_Node& xmLBinarySpecies)
{
    std::string xname = xmLBinarySpecies.name();
}

std::string MolarityIonicVPSSTP::report(bool show_thermo, doublereal threshold) const
{
    fmt::memory_buffer b;
    try {
        if (name() != "") {
            format_to(b, "\n  {}:\n", name());
        }
        format_to(b, "\n");
        format_to(b, "       temperature    {:12.6g}  K\n", temperature());
        format_to(b, "          pressure    {:12.6g}  Pa\n", pressure());
        format_to(b, "           density    {:12.6g}  kg/m^3\n", density());
        format_to(b, "  mean mol. weight    {:12.6g}  amu\n", meanMolecularWeight());

        doublereal phi = electricPotential();
        format_to(b, "         potential    {:12.6g}  V\n", phi);

        vector_fp x(m_kk);
        vector_fp molal(m_kk);
        vector_fp mu(m_kk);
        vector_fp muss(m_kk);
        vector_fp acMolal(m_kk);
        vector_fp actMolal(m_kk);
        getMoleFractions(&x[0]);

        getChemPotentials(&mu[0]);
        getStandardChemPotentials(&muss[0]);
        getActivities(&actMolal[0]);

        if (show_thermo) {
            format_to(b, "\n");
            format_to(b, "                          1 kg            1 kmol\n");
            format_to(b, "                       -----------      ------------\n");
            format_to(b, "          enthalpy    {:12.6g}     {:12.4g}     J\n",
                    enthalpy_mass(), enthalpy_mole());
            format_to(b, "   internal energy    {:12.6g}     {:12.4g}     J\n",
                    intEnergy_mass(), intEnergy_mole());
            format_to(b, "           entropy    {:12.6g}     {:12.4g}     J/K\n",
                    entropy_mass(), entropy_mole());
            format_to(b, "    Gibbs function    {:12.6g}     {:12.4g}     J\n",
                    gibbs_mass(), gibbs_mole());
            format_to(b, " heat capacity c_p    {:12.6g}     {:12.4g}     J/K\n",
                    cp_mass(), cp_mole());
            try {
                format_to(b, " heat capacity c_v    {:12.6g}     {:12.4g}     J/K\n",
                        cv_mass(), cv_mole());
            } catch (NotImplementedError&) {
                format_to(b, " heat capacity c_v    <not implemented>\n");
            }
        }
    } catch (CanteraError& e) {
        return to_string(b) + e.what();
    }
    return to_string(b);
}

}
