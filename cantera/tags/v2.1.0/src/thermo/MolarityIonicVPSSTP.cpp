/**
 *  @file MolarityIonicVPSSTP.cpp
 *   Definitions for intermediate ThermoPhase object for phases which
 *   employ excess gibbs free energy formulations
 *  (see \ref thermoprops
 * and class \link Cantera::MolarityIonicVPSSTP MolarityIonicVPSSTP\endlink).
 *
 * Header file for a derived class of ThermoPhase that handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties that are further based upon expressions
 * for the excess gibbs free energy expressed as a function of
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
    GibbsExcessVPSSTP(),
    PBType_(PBTYPE_PASSTHROUGH),
    numPBSpecies_(m_kk),
    indexSpecialSpecies_(npos),
    numCationSpecies_(0),
    numAnionSpecies_(0),
    numPassThroughSpecies_(0),
    neutralPBindexStart(0)
{
}

MolarityIonicVPSSTP::MolarityIonicVPSSTP(const std::string& inputFile,
        const std::string& id_) :
    GibbsExcessVPSSTP(),
    PBType_(PBTYPE_PASSTHROUGH),
    numPBSpecies_(m_kk),
    indexSpecialSpecies_(npos),
    numCationSpecies_(0),
    numAnionSpecies_(0),
    numPassThroughSpecies_(0),
    neutralPBindexStart(0)
{
    initThermoFile(inputFile, id_);
}

MolarityIonicVPSSTP::MolarityIonicVPSSTP(XML_Node& phaseRoot,
        const std::string& id_) :
    GibbsExcessVPSSTP(),
    PBType_(PBTYPE_PASSTHROUGH),
    numPBSpecies_(m_kk),
    indexSpecialSpecies_(npos),
    numCationSpecies_(0),
    numAnionSpecies_(0),
    numPassThroughSpecies_(0),
    neutralPBindexStart(0)
{
    importPhase(*findXMLPhase(&phaseRoot, id_), this);
}

MolarityIonicVPSSTP::MolarityIonicVPSSTP(const MolarityIonicVPSSTP& b) :
    GibbsExcessVPSSTP(),
    PBType_(PBTYPE_PASSTHROUGH),
    numPBSpecies_(m_kk),
    indexSpecialSpecies_(npos),
    numCationSpecies_(0),
    numAnionSpecies_(0),
    numPassThroughSpecies_(0),
    neutralPBindexStart(0)
{
    *this = operator=(b);
}

MolarityIonicVPSSTP& MolarityIonicVPSSTP::
operator=(const MolarityIonicVPSSTP& b)
{
    if (&b != this) {
        GibbsExcessVPSSTP::operator=(b);
    }

    PBType_                     = b.PBType_;
    numPBSpecies_               = b.numPBSpecies_;
    indexSpecialSpecies_        = b.indexSpecialSpecies_;
    PBMoleFractions_            = b.PBMoleFractions_;
    cationList_                 = b.cationList_;
    numCationSpecies_           = b.numCationSpecies_;
    anionList_                  = b.anionList_;
    numAnionSpecies_            = b.numAnionSpecies_;
    passThroughList_            = b.passThroughList_;
    numPassThroughSpecies_      = b.numPassThroughSpecies_;
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
 *  -------------- Utilities -------------------------------
 */

int MolarityIonicVPSSTP::eosType() const
{
    return 0;
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
    size_t k;
    size_t kCat;
    size_t kMax;
    doublereal sumCat;
    doublereal sumAnion;
    doublereal chP, chM;
    doublereal sum = 0.0;
    doublereal sumMax;
    switch (PBType_) {
    case PBTYPE_PASSTHROUGH:
        for (k = 0; k < m_kk; k++) {
            PBMoleFractions_[k] = moleFractions_[k];
        }
        break;
    case PBTYPE_SINGLEANION:
        sumCat = 0.0;
        sumAnion = 0.0;
        for (k = 0; k < m_kk; k++) {
            moleFractionsTmp_[k] = moleFractions_[k];
        }
        kMax = npos;
        sumMax = 0.0;
        for (k = 0; k < cationList_.size(); k++) {
            kCat = cationList_[k];
            chP = m_speciesCharge[kCat];
            if (moleFractions_[kCat] > sumMax) {
                kMax = k;
                sumMax = moleFractions_[kCat];
            }
            sumCat += chP * moleFractions_[kCat];
        }
        k = anionList_[0];
        chM = m_speciesCharge[k];
        sumAnion = moleFractions_[k] * chM;
        sum = sumCat - sumAnion;
        if (fabs(sum) > 1.0E-16) {
            moleFractionsTmp_[cationList_[kMax]] -= sum / m_speciesCharge[kMax];
            sum = 0.0;
            for (k = 0; k < numCationSpecies_; k++) {
                sum +=  moleFractionsTmp_[k];
            }
            for (k = 0; k < numCationSpecies_; k++) {
                moleFractionsTmp_[k]/= sum;
            }
        }

        for (k = 0; k < numCationSpecies_; k++) {
            PBMoleFractions_[k] = moleFractionsTmp_[cationList_[k]];
        }
        for (k = 0; k <  numPassThroughSpecies_; k++) {
            PBMoleFractions_[neutralPBindexStart + k] = moleFractions_[passThroughList_[k]];
        }

        sum = std::max(0.0, PBMoleFractions_[0]);
        for (k = 1; k < numPBSpecies_; k++) {
            sum += PBMoleFractions_[k];

        }
        for (k = 0; k < numPBSpecies_; k++) {
            PBMoleFractions_[k] /= sum;
        }

        break;
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

/*
 * ------------ Partial Molar Properties of the Solution ------------
 */

doublereal MolarityIonicVPSSTP::err(const std::string& msg) const
{
    throw CanteraError("MolarityIonicVPSSTP","Base class method "
                       +msg+" called. Equation of state type: "+int2str(eosType()));
    return 0;
}

void MolarityIonicVPSSTP::initThermo()
{
    GibbsExcessVPSSTP::initThermo();
    initLengths();
    /*
     *  Go find the list of cations and anions
     */
    double ch;
    numCationSpecies_ = 0;
    cationList_.clear();
    anionList_.clear();
    passThroughList_.clear();
    for (size_t k = 0; k < m_kk; k++) {
        ch = m_speciesCharge[k];
        if (ch > 0.0) {
            cationList_.push_back(k);
            numCationSpecies_++;
        } else if (ch < 0.0) {
            anionList_.push_back(k);
            numAnionSpecies_++;
        } else {
            passThroughList_.push_back(k);
            numPassThroughSpecies_++;
        }
    }
    numPBSpecies_ = numCationSpecies_ + numAnionSpecies_ - 1;
    neutralPBindexStart = numPBSpecies_;
    PBType_ = PBTYPE_MULTICATIONANION;
    if (numAnionSpecies_ == 1) {
        PBType_ = PBTYPE_SINGLEANION;
    } else if (numCationSpecies_ == 1) {
        PBType_ = PBTYPE_SINGLECATION;
    }
    if (numAnionSpecies_ == 0 && numCationSpecies_ == 0) {
        PBType_ = PBTYPE_PASSTHROUGH;
    }
}

void  MolarityIonicVPSSTP::initLengths()
{
    m_kk = nSpecies();
    moleFractionsTmp_.resize(m_kk);
}

void MolarityIonicVPSSTP::initThermoXML(XML_Node& phaseNode, const std::string& id)
{
    std::string subname = "MolarityIonicVPSSTP::initThermoXML";
    std::string stemp;

    if ((int) id.size() > 0) {
        string idp = phaseNode.id();
        if (idp != id) {
            throw CanteraError(subname, "phasenode and Id are incompatible");
        }
    }

    /*
     * Check on the thermo field. Must have one of:
     * <thermo model="MolarityIonicVPSS" />
     * <thermo model="MolarityIonicVPSSTP" />
     */
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError(subname, "no thermo XML node");
    }
    XML_Node& thermoNode = phaseNode.child("thermo");
    std::string mStringa = thermoNode.attrib("model");
    std::string mString = lowercase(mStringa);
    if (mString != "molarityionicvpss" && mString != "molarityionicvpsstp") {
        throw CanteraError(subname.c_str(),
                           "Unknown thermo model: " + mStringa + " - This object only knows \"MolarityIonicVPSSTP\" ");
    }

    /*
     * Go get all of the coefficients and factors in the
     * activityCoefficients XML block
     */
    XML_Node* acNodePtr = 0;
    if (thermoNode.hasChild("activityCoefficients")) {
        XML_Node& acNode = thermoNode.child("activityCoefficients");
        acNodePtr = &acNode;
        mStringa = acNode.attrib("model");
        mString = lowercase(mStringa);
        // if (mString != "redlich-kister") {
        //   throw CanteraError(subname.c_str(),
        //        "Unknown activity coefficient model: " + mStringa);
        //}
        size_t n = acNodePtr->nChildren();
        for (size_t i = 0; i < n; i++) {
            XML_Node& xmlACChild = acNodePtr->child(i);
            stemp = xmlACChild.name();
            std::string nodeName = lowercase(stemp);
            /*
             * Process a binary interaction
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

void MolarityIonicVPSSTP::readXMLBinarySpecies(XML_Node& xmLBinarySpecies)
{
    std::string xname = xmLBinarySpecies.name();

}

std::string MolarityIonicVPSSTP::report(bool show_thermo) const
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

        size_t kk = nSpecies();
        vector_fp x(kk);
        vector_fp molal(kk);
        vector_fp mu(kk);
        vector_fp muss(kk);
        vector_fp acMolal(kk);
        vector_fp actMolal(kk);
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
