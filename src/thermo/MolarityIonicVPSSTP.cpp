/**
 *  @file MolarityIonicVPSSTP.cpp
 *   Definitions for intermediate ThermoPhase object for phases which
 *   employ excess Gibbs free energy formulations
 *  (see \ref thermoprops
 * and class \link Cantera::MolarityIonicVPSSTP MolarityIonicVPSSTP\endlink).
 *
 * Header file for a derived class of ThermoPhase that handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties that are further based upon expressions
 * for the excess Gibbs free energy expressed as a function of
 * the mole fractions.
 */
/*
 * Copyright (2009) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
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
    indexSpecialSpecies_(npos),
    neutralPBindexStart(0)
{
}

MolarityIonicVPSSTP::MolarityIonicVPSSTP(const std::string& inputFile,
        const std::string& id_) :
    PBType_(PBTYPE_PASSTHROUGH),
    numPBSpecies_(m_kk),
    indexSpecialSpecies_(npos),
    neutralPBindexStart(0)
{
    initThermoFile(inputFile, id_);
}

MolarityIonicVPSSTP::MolarityIonicVPSSTP(XML_Node& phaseRoot,
        const std::string& id_) :
    PBType_(PBTYPE_PASSTHROUGH),
    numPBSpecies_(m_kk),
    indexSpecialSpecies_(npos),
    neutralPBindexStart(0)
{
    importPhase(*findXMLPhase(&phaseRoot, id_), this);
}

MolarityIonicVPSSTP::MolarityIonicVPSSTP(const MolarityIonicVPSSTP& b) :
    PBType_(PBTYPE_PASSTHROUGH),
    numPBSpecies_(m_kk),
    indexSpecialSpecies_(npos),
    neutralPBindexStart(0)
{
    *this = b;
}

MolarityIonicVPSSTP& MolarityIonicVPSSTP::operator=(const MolarityIonicVPSSTP& b)
{
    if (&b != this) {
        GibbsExcessVPSSTP::operator=(b);
    }

    PBType_                     = b.PBType_;
    numPBSpecies_               = b.numPBSpecies_;
    indexSpecialSpecies_        = b.indexSpecialSpecies_;
    PBMoleFractions_            = b.PBMoleFractions_;
    cationList_                 = b.cationList_;
    anionList_                  = b.anionList_;
    passThroughList_            = b.passThroughList_;
    neutralPBindexStart         = b.neutralPBindexStart;
    moleFractionsTmp_           = b.moleFractionsTmp_;

    return *this;
}

ThermoPhase*
MolarityIonicVPSSTP::duplMyselfAsThermoPhase() const
{
    return new MolarityIonicVPSSTP(*this);
}

/*
 * - Activities, Standard States, Activity Concentrations -----------
 */

void MolarityIonicVPSSTP::getLnActivityCoefficients(doublereal* lnac) const
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

void MolarityIonicVPSSTP::getChemPotentials(doublereal* mu) const
{
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
    doublereal RT = GasConstant * temperature();
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(moleFractions_[k], SmallNumber);
        mu[k] += RT * (log(xx) + lnActCoeff_Scaled_[k]);
    }
}

void MolarityIonicVPSSTP::getElectrochemPotentials(doublereal* mu) const
{
    getChemPotentials(mu);
    double ve = Faraday * electricPotential();
    for (size_t k = 0; k < m_kk; k++) {
        mu[k] += ve*charge(k);
    }
}

void MolarityIonicVPSSTP::getPartialMolarEnthalpies(doublereal* hbar) const
{
    /*
     * Get the nondimensional standard state enthalpies
     */
    getEnthalpy_RT(hbar);
    /*
     * dimensionalize it.
     */
    double T = temperature();
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] *= GasConstant * T;
    }
    /*
     * Update the activity coefficients, This also update the
     * internally stored molalities.
     */
    s_update_lnActCoeff();
    s_update_dlnActCoeff_dT();
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] -= GasConstant * T * T * dlnActCoeffdT_Scaled_[k];
    }
}

void MolarityIonicVPSSTP::getPartialMolarCp(doublereal* cpbar) const
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

void MolarityIonicVPSSTP::getPartialMolarEntropies(doublereal* sbar) const
{
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
        double xx = std::max(moleFractions_[k], SmallNumber);
        sbar[k] += - lnActCoeff_Scaled_[k] -log(xx) - T * dlnActCoeffdT_Scaled_[k];
    }
    /*
     * dimensionalize it.
     */
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] *= GasConstant;
    }
}

void MolarityIonicVPSSTP::getPartialMolarVolumes(doublereal* vbar) const
{
    /*
     * Get the standard state values in m^3 kmol-1
     */
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
                sum +=  moleFractionsTmp_[k];
            }
            for (size_t k = 0; k < cationList_.size(); k++) {
                moleFractionsTmp_[k]/= sum;
            }
        }

        for (size_t k = 0; k < cationList_.size(); k++) {
            PBMoleFractions_[k] = moleFractionsTmp_[cationList_[k]];
        }
        for (size_t k = 0; k <  passThroughList_.size(); k++) {
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

        break;

    case PBTYPE_MULTICATIONANION:
        throw CanteraError("eosType", "Unknown type");

        break;
    default:
        throw CanteraError("eosType", "Unknown type");
        break;

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

void  MolarityIonicVPSSTP::s_update_dlnActCoeff_dX_() const
{
}

void MolarityIonicVPSSTP::initThermo()
{
    GibbsExcessVPSSTP::initThermo();
    initLengths();
    /*
     *  Go find the list of cations and anions
     */
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

void  MolarityIonicVPSSTP::initLengths()
{
    moleFractionsTmp_.resize(m_kk);
}

void MolarityIonicVPSSTP::initThermoXML(XML_Node& phaseNode, const std::string& id)
{
    if ((int) id.size() > 0 && phaseNode.id() != id) {
        throw CanteraError("MolarityIonicVPSSTP::initThermoXML",
                           "phasenode and Id are incompatible");
    }

    /*
     * Check on the thermo field. Must have one of:
     * <thermo model="MolarityIonicVPSS" />
     * <thermo model="MolarityIonicVPSSTP" />
     */
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError("MolarityIonicVPSSTP::initThermoXML",
                           "no thermo XML node");
    }
    XML_Node& thermoNode = phaseNode.child("thermo");
    std::string mStringa = thermoNode.attrib("model");
    std::string mString = lowercase(mStringa);
    if (mString != "molarityionicvpss" && mString != "molarityionicvpsstp") {
        throw CanteraError("MolarityIonicVPSSTP::initThermoXML",
                           "Unknown thermo model: " + mStringa + " - This object only knows \"MolarityIonicVPSSTP\" ");
    }

    /*
     * Go get all of the coefficients and factors in the
     * activityCoefficients XML block
     */
    if (thermoNode.hasChild("activityCoefficients")) {
        XML_Node& acNode = thermoNode.child("activityCoefficients");
        for (size_t i = 0; i < acNode.nChildren(); i++) {
            XML_Node& xmlACChild = acNode.child(i);
            /*
             * Process a binary interaction
             */
            if (lowercase(xmlACChild.name()) == "binaryneutralspeciesparameters") {
                readXMLBinarySpecies(xmlACChild);
            }
        }
    }


    /*
     * Go down the chain
     */
    GibbsExcessVPSSTP::initThermoXML(phaseNode, id);
}

void MolarityIonicVPSSTP::readXMLBinarySpecies(XML_Node& xmLBinarySpecies)
{
    std::string xname = xmLBinarySpecies.name();
}

std::string MolarityIonicVPSSTP::report(bool show_thermo, doublereal threshold) const
{
    char p[800];
    string s = "";
    try {
        if (name() != "") {
            sprintf(p, " \n  %s:\n", name().c_str());
            s += p;
        }
        sprintf(p, " \n       temperature    %12.6g  K\n", temperature());
        s += p;
        sprintf(p, "          pressure    %12.6g  Pa\n", pressure());
        s += p;
        sprintf(p, "           density    %12.6g  kg/m^3\n", density());
        s += p;
        sprintf(p, "  mean mol. weight    %12.6g  amu\n", meanMolecularWeight());
        s += p;

        doublereal phi = electricPotential();
        sprintf(p, "         potential    %12.6g  V\n", phi);
        s += p;

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
            sprintf(p, " \n");
            s += p;
            sprintf(p, "                          1 kg            1 kmol\n");
            s += p;
            sprintf(p, "                       -----------      ------------\n");
            s += p;
            sprintf(p, "          enthalpy    %12.6g     %12.4g     J\n",
                    enthalpy_mass(), enthalpy_mole());
            s += p;
            sprintf(p, "   internal energy    %12.6g     %12.4g     J\n",
                    intEnergy_mass(), intEnergy_mole());
            s += p;
            sprintf(p, "           entropy    %12.6g     %12.4g     J/K\n",
                    entropy_mass(), entropy_mole());
            s += p;
            sprintf(p, "    Gibbs function    %12.6g     %12.4g     J\n",
                    gibbs_mass(), gibbs_mole());
            s += p;
            sprintf(p, " heat capacity c_p    %12.6g     %12.4g     J/K\n",
                    cp_mass(), cp_mole());
            s += p;
            try {
                sprintf(p, " heat capacity c_v    %12.6g     %12.4g     J/K\n",
                        cv_mass(), cv_mole());
                s += p;
            } catch (CanteraError& e) {
                e.save();
                sprintf(p, " heat capacity c_v    <not implemented>       \n");
                s += p;
            }
        }

    } catch (CanteraError& e) {
        e.save();
    }
    return s;
}

}
