/**
 *  @file IonsFromNeutralVPSSTP.cpp
 *   Definitions for the object which treats ionic liquids as made of ions as species
 *   even though the thermodynamics is obtained from the neutral molecule representation.
 *  (see \ref thermoprops
 *   and class \link Cantera::IonsFromNeutralVPSSTP IonsFromNeutralVPSSTP\endlink).
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
#include "cantera/thermo/IonsFromNeutralVPSSTP.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/PDSS_IonsFromNeutral.h"
#include "cantera/base/stringUtils.h"

#include <fstream>

using namespace std;

namespace Cantera
{

IonsFromNeutralVPSSTP::IonsFromNeutralVPSSTP() :
    ionSolnType_(cIonSolnType_SINGLEANION),
    numNeutralMoleculeSpecies_(0),
    indexSpecialSpecies_(npos),
    indexSecondSpecialSpecies_(npos),
    neutralMoleculePhase_(0),
    geThermo(0),
    IOwnNThermoPhase_(true)
{
}

IonsFromNeutralVPSSTP::IonsFromNeutralVPSSTP(const std::string& inputFile,
        const std::string& id_,
        ThermoPhase* neutralPhase) :
    ionSolnType_(cIonSolnType_SINGLEANION),
    numNeutralMoleculeSpecies_(0),
    indexSpecialSpecies_(npos),
    indexSecondSpecialSpecies_(npos),
    neutralMoleculePhase_(neutralPhase),
    IOwnNThermoPhase_(true)
{
    if (neutralPhase) {
        IOwnNThermoPhase_ = false;
    }
    constructPhaseFile(inputFile, id_);
    geThermo = dynamic_cast<GibbsExcessVPSSTP*>(neutralMoleculePhase_);
}

IonsFromNeutralVPSSTP::IonsFromNeutralVPSSTP(XML_Node& phaseRoot,
        const std::string& id_, ThermoPhase* neutralPhase) :
    ionSolnType_(cIonSolnType_SINGLEANION),
    numNeutralMoleculeSpecies_(0),
    indexSpecialSpecies_(npos),
    indexSecondSpecialSpecies_(npos),
    neutralMoleculePhase_(neutralPhase),
    IOwnNThermoPhase_(true)
{
    if (neutralPhase) {
        IOwnNThermoPhase_ = false;
    }
    constructPhaseXML(phaseRoot, id_);
    geThermo = dynamic_cast<GibbsExcessVPSSTP*>(neutralMoleculePhase_);
    y_.resize(numNeutralMoleculeSpecies_,0.0);
    size_t numNeutMolSpec = geThermo->nSpecies();
    dlnActCoeff_NeutralMolecule_.resize(numNeutMolSpec);
    dX_NeutralMolecule_.resize(numNeutMolSpec);
}

IonsFromNeutralVPSSTP::IonsFromNeutralVPSSTP(const IonsFromNeutralVPSSTP& b) :
    ionSolnType_(cIonSolnType_SINGLEANION),
    numNeutralMoleculeSpecies_(0),
    indexSpecialSpecies_(npos),
    indexSecondSpecialSpecies_(npos),
    neutralMoleculePhase_(0),
    geThermo(0),
    IOwnNThermoPhase_(true)
{
    IonsFromNeutralVPSSTP::operator=(b);
}

IonsFromNeutralVPSSTP&
IonsFromNeutralVPSSTP::operator=(const IonsFromNeutralVPSSTP& b)
{
    if (&b == this) {
        return *this;
    }

    /*
     *  If we own the underlying neutral molecule phase, then we do a deep
     *  copy. If not, we do a shallow copy. We get a valid pointer for
     *  neutralMoleculePhase_ first, because we need it to assign the pointers
     *  within the PDSS_IonsFromNeutral object. which is done in the
     *  GibbsExcessVPSSTP::operator=(b) step.
     */
    if (IOwnNThermoPhase_) {
        if (b.neutralMoleculePhase_) {
            delete neutralMoleculePhase_;
            neutralMoleculePhase_   = (b.neutralMoleculePhase_)->duplMyselfAsThermoPhase();
        } else {
            neutralMoleculePhase_   = 0;
        }
    } else {
        neutralMoleculePhase_     = b.neutralMoleculePhase_;
    }
    geThermo = dynamic_cast<GibbsExcessVPSSTP*>(neutralMoleculePhase_);

    GibbsExcessVPSSTP::operator=(b);


    ionSolnType_                = b.ionSolnType_;
    numNeutralMoleculeSpecies_  = b.numNeutralMoleculeSpecies_;
    indexSpecialSpecies_        = b.indexSpecialSpecies_;
    indexSecondSpecialSpecies_  = b.indexSecondSpecialSpecies_;
    fm_neutralMolec_ions_       = b.fm_neutralMolec_ions_;
    fm_invert_ionForNeutral     = b.fm_invert_ionForNeutral;
    NeutralMolecMoleFractions_  = b.NeutralMolecMoleFractions_;
    cationList_                 = b.cationList_;
    anionList_                  = b.anionList_;
    passThroughList_            = b.passThroughList_;

    y_                          = b.y_;
    dlnActCoeff_NeutralMolecule_ = b.dlnActCoeff_NeutralMolecule_;
    dX_NeutralMolecule_         = b.dX_NeutralMolecule_;

    IOwnNThermoPhase_           = b.IOwnNThermoPhase_;
    moleFractionsTmp_           = b.moleFractionsTmp_;
    muNeutralMolecule_          = b.muNeutralMolecule_;
    //  gammaNeutralMolecule_       = b.gammaNeutralMolecule_;
    lnActCoeff_NeutralMolecule_ = b.lnActCoeff_NeutralMolecule_;
    dlnActCoeffdT_NeutralMolecule_ = b.dlnActCoeffdT_NeutralMolecule_;
    dlnActCoeffdlnX_diag_NeutralMolecule_ = b.dlnActCoeffdlnX_diag_NeutralMolecule_;
    dlnActCoeffdlnN_diag_NeutralMolecule_ = b.dlnActCoeffdlnN_diag_NeutralMolecule_;
    dlnActCoeffdlnN_NeutralMolecule_ = b.dlnActCoeffdlnN_NeutralMolecule_;

    return *this;
}

IonsFromNeutralVPSSTP::~IonsFromNeutralVPSSTP()
{
    if (IOwnNThermoPhase_) {
        delete neutralMoleculePhase_;
        neutralMoleculePhase_ = 0;
    }
}

ThermoPhase*
IonsFromNeutralVPSSTP::duplMyselfAsThermoPhase() const
{
    return new IonsFromNeutralVPSSTP(*this);
}

void IonsFromNeutralVPSSTP::constructPhaseFile(std::string inputFile, std::string id_)
{

    if (inputFile.size() == 0) {
        throw CanteraError("MargulesVPSSTP:constructPhaseFile",
                           "input file is null");
    }
    string path = findInputFile(inputFile);
    std::ifstream fin(path.c_str());
    if (!fin) {
        throw CanteraError("MargulesVPSSTP:constructPhaseFile","could not open "
                           +path+" for reading.");
    }
    /*
     * The phase object automatically constructs an XML object.
     * Use this object to store information.
     */
    XML_Node fxml;
    fxml.build(fin);
    XML_Node* fxml_phase = findXMLPhase(&fxml, id_);
    if (!fxml_phase) {
        throw CanteraError("MargulesVPSSTP:constructPhaseFile",
                           "ERROR: Can not find phase named " +
                           id_ + " in file named " + inputFile);
    }
    setXMLdata(*fxml_phase);
    constructPhaseXML(*fxml_phase, id_);
}

void IonsFromNeutralVPSSTP::constructPhaseXML(XML_Node& phaseNode, std::string id_)
{
    if (id_.size() > 0) {
        if (phaseNode.id() != id_) {
            throw CanteraError("IonsFromNeutralVPSSTP::constructPhaseXML",
                               "phasenode and Id are incompatible");
        }
    }

    /*
     * Find the thermo XML node
     */
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError("IonsFromNeutralVPSSTP::constructPhaseXML",
                           "no thermo XML node");
    }
    XML_Node& thermoNode = phaseNode.child("thermo");



    /*
     * Make sure that the thermo model is IonsFromNeutralMolecule
     */
    string formString = lowercase(thermoNode.attrib("model"));
    if (formString != "ionsfromneutralmolecule") {
        throw CanteraError("IonsFromNeutralVPSSTP::constructPhaseXML",
                           "model name isn't IonsFromNeutralMolecule: " + formString);
    }

    /*
     * Find the Neutral Molecule Phase
     */
    if (!thermoNode.hasChild("neutralMoleculePhase")) {
        throw CanteraError("IonsFromNeutralVPSSTP::constructPhaseXML",
                           "no neutralMoleculePhase XML node");
    }
    XML_Node& neutralMoleculeNode = thermoNode.child("neutralMoleculePhase");

    XML_Node* neut_ptr = get_XML_Node(neutralMoleculeNode["datasrc"], 0);
    if (!neut_ptr) {
        throw CanteraError("IonsFromNeutralVPSSTP::constructPhaseXML",
                           "neut_ptr = 0");
    }

    /*
     *  Create the neutralMolecule ThermoPhase if we haven't already
     */
    if (!neutralMoleculePhase_) {
        neutralMoleculePhase_  = newPhase(*neut_ptr);
    }

    /*
     * Call the Cantera importPhase() function. This will import
     * all of the species into the phase. This will also handle
     * all of the solvent and solute standard states
     */
    importPhase(phaseNode, this);
}

/*
 *  -------------- Utilities -------------------------------
 */

int IonsFromNeutralVPSSTP::eosType() const
{
    return cIonsFromNeutral;
}

/*
 * ------------ Molar Thermodynamic Properties ----------------------
 */

doublereal IonsFromNeutralVPSSTP::enthalpy_mole() const
{
    getPartialMolarEnthalpies(DATA_PTR(m_pp));
    return mean_X(m_pp);
}

doublereal IonsFromNeutralVPSSTP::entropy_mole() const
{
    getPartialMolarEntropies(DATA_PTR(m_pp));
    return mean_X(m_pp);
}

doublereal IonsFromNeutralVPSSTP::gibbs_mole() const
{
    getChemPotentials(DATA_PTR(m_pp));
    return mean_X(m_pp);
}

doublereal IonsFromNeutralVPSSTP::cp_mole() const
{
    getPartialMolarCp(DATA_PTR(m_pp));
    return mean_X(m_pp);
}

doublereal IonsFromNeutralVPSSTP::cv_mole() const
{
    // Need to revisit this, as it is wrong
    getPartialMolarCp(DATA_PTR(m_pp));
    return mean_X(m_pp);
}

/*
 * - Activities, Standard States, Activity Concentrations -----------
 */

void IonsFromNeutralVPSSTP::getDissociationCoeffs(vector_fp& coeffs,
        vector_fp& charges, std::vector<size_t>& neutMolIndex) const
{
    coeffs = fm_neutralMolec_ions_;
    charges = m_speciesCharge;
    neutMolIndex = fm_invert_ionForNeutral;
}

void IonsFromNeutralVPSSTP::getActivityCoefficients(doublereal* ac) const
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

/*
 * ---------  Partial Molar Properties of the Solution -------------
 */

void IonsFromNeutralVPSSTP::getChemPotentials(doublereal* mu) const
{
    size_t icat, jNeut;
    doublereal xx, fact2;

    /*
     * Get the standard chemical potentials of netural molecules
     */
    neutralMoleculePhase_->getStandardChemPotentials(DATA_PTR(muNeutralMolecule_));

    doublereal RT_ = GasConstant * temperature();

    switch (ionSolnType_) {
    case cIonSolnType_PASSTHROUGH:
        neutralMoleculePhase_->getChemPotentials(mu);
        break;
    case cIonSolnType_SINGLEANION:
        neutralMoleculePhase_->getLnActivityCoefficients(DATA_PTR(lnActCoeff_NeutralMolecule_));

        fact2 = 2.0 * RT_ * log(2.0);

        // Do the cation list
        for (size_t k = 0; k < cationList_.size(); k++) {
            //! Get the id for the next cation
            icat = cationList_[k];
            jNeut = fm_invert_ionForNeutral[icat];
            xx = std::max(SmallNumber, moleFractions_[icat]);
            mu[icat] = muNeutralMolecule_[jNeut] + fact2 + RT_ * (lnActCoeff_NeutralMolecule_[jNeut] + log(xx));
        }

        // Do the anion list
        icat = anionList_[0];
        jNeut = fm_invert_ionForNeutral[icat];
        xx = std::max(SmallNumber, moleFractions_[icat]);
        mu[icat] = RT_ * log(xx);

        // Do the list of neutral molecules
        for (size_t k = 0; k < passThroughList_.size(); k++) {
            icat = passThroughList_[k];
            jNeut = fm_invert_ionForNeutral[icat];
            xx = std::max(SmallNumber, moleFractions_[icat]);
            mu[icat] =  muNeutralMolecule_[jNeut]  + RT_ * (lnActCoeff_NeutralMolecule_[jNeut] + log(xx));
        }
        break;

    case cIonSolnType_SINGLECATION:
        throw CanteraError("eosType", "Unknown type");
        break;
    case cIonSolnType_MULTICATIONANION:
        throw CanteraError("eosType", "Unknown type");
        break;
    default:
        throw CanteraError("eosType", "Unknown type");
        break;
    }
}

void IonsFromNeutralVPSSTP::getPartialMolarEnthalpies(doublereal* hbar) const
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
    s_update_dlnActCoeffdT();
    double RTT = RT * T;
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] -= RTT * dlnActCoeffdT_Scaled_[k];
    }
}

void IonsFromNeutralVPSSTP::getPartialMolarEntropies(doublereal* sbar) const
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
    s_update_dlnActCoeffdT();

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

void IonsFromNeutralVPSSTP::getdlnActCoeffdlnX_diag(doublereal* dlnActCoeffdlnX_diag) const
{
    s_update_lnActCoeff();
    s_update_dlnActCoeff_dlnX_diag();

    for (size_t k = 0; k < m_kk; k++) {
        dlnActCoeffdlnX_diag[k] = dlnActCoeffdlnX_diag_[k];
    }
}

void IonsFromNeutralVPSSTP::getdlnActCoeffdlnN_diag(doublereal* dlnActCoeffdlnN_diag) const
{
    s_update_lnActCoeff();
    s_update_dlnActCoeff_dlnN_diag();

    for (size_t k = 0; k < m_kk; k++) {
        dlnActCoeffdlnN_diag[k] = dlnActCoeffdlnN_diag_[k];
    }
}

void IonsFromNeutralVPSSTP::getdlnActCoeffdlnN(const size_t ld, doublereal* dlnActCoeffdlnN)
{
    s_update_lnActCoeff();
    s_update_dlnActCoeff_dlnN();
    double* data =  & dlnActCoeffdlnN_(0,0);
    for (size_t k = 0; k < m_kk; k++) {
        for (size_t m = 0; m < m_kk; m++) {
            dlnActCoeffdlnN[ld * k + m] = data[m_kk * k + m];
        }
    }
}

void IonsFromNeutralVPSSTP::setTemperature(const doublereal temp)
{
    IonsFromNeutralVPSSTP::setState_TP(temp, pressure());
}

void IonsFromNeutralVPSSTP::setPressure(doublereal p)
{
    IonsFromNeutralVPSSTP::setState_TP(temperature(), p);
}

void IonsFromNeutralVPSSTP::setState_TP(doublereal t, doublereal p)
{
    /*
     *  This is a two phase process. First, we calculate the standard states
     *  within the neutral molecule phase.
     */
    neutralMoleculePhase_->setState_TP(t, p);
    VPStandardStateTP::setState_TP(t,p);

    /*
     * Calculate the partial molar volumes, and then the density of the fluid
     */
    Phase::setDensity(neutralMoleculePhase_->density());
}

void IonsFromNeutralVPSSTP::calcIonMoleFractions(doublereal* const mf) const
{
    /*
     * Download the neutral mole fraction vector into the
     * vector, NeutralMolecMoleFractions_[]
     */
    neutralMoleculePhase_->getMoleFractions(DATA_PTR(NeutralMolecMoleFractions_));

    // Zero the mole fractions
    for (size_t k = 0; k < m_kk; k++) {
        mf[k] = 0.0;
    }

    /*
     *  Use the formula matrix to calculate the relative mole numbers.
     */
    for (size_t jNeut = 0; jNeut <  numNeutralMoleculeSpecies_; jNeut++) {
        for (size_t k = 0; k < m_kk; k++) {
            double fmij = fm_neutralMolec_ions_[k + jNeut * m_kk];
            mf[k] += fmij * NeutralMolecMoleFractions_[jNeut];
        }
    }

    /*
     * Normalize the new mole fractions
     */
    doublereal sum = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        sum += mf[k];
    }
    for (size_t k = 0; k < m_kk; k++) {
        mf[k] /= sum;
    }

}

void IonsFromNeutralVPSSTP::calcNeutralMoleculeMoleFractions() const
{
    size_t icat, jNeut;
    doublereal fmij;
    doublereal sum = 0.0;

    //! Zero the vector we are trying to find.
    for (size_t k = 0; k < numNeutralMoleculeSpecies_; k++) {
        NeutralMolecMoleFractions_[k] = 0.0;
    }
    if (DEBUG_MODE_ENABLED) {
        sum = -1.0;
        for (size_t k = 0; k < m_kk; k++) {
            sum += moleFractions_[k];
        }
        if (fabs(sum) > 1.0E-11)  {
            throw CanteraError("IonsFromNeutralVPSSTP::calcNeutralMoleculeMoleFractions",
                               "molefracts don't sum to one: " + fp2str(sum));
        }
    }

    switch (ionSolnType_) {

    case cIonSolnType_PASSTHROUGH:
        for (size_t k = 0; k < m_kk; k++) {
            NeutralMolecMoleFractions_[k] = moleFractions_[k];
        }
        break;

    case cIonSolnType_SINGLEANION:
        for (size_t k = 0; k < numNeutralMoleculeSpecies_; k++) {
            NeutralMolecMoleFractions_[k] = 0.0;
        }

        for (size_t k = 0; k < cationList_.size(); k++) {
            //! Get the id for the next cation
            icat = cationList_[k];
            jNeut = fm_invert_ionForNeutral[icat];
            if (jNeut != npos) {
                fmij =  fm_neutralMolec_ions_[icat + jNeut * m_kk];
                AssertTrace(fmij != 0.0);
                NeutralMolecMoleFractions_[jNeut] += moleFractions_[icat] / fmij;
            }
        }

        for (size_t k = 0; k <  passThroughList_.size(); k++) {
            icat = passThroughList_[k];
            jNeut = fm_invert_ionForNeutral[icat];
            fmij = fm_neutralMolec_ions_[ icat + jNeut * m_kk];
            NeutralMolecMoleFractions_[jNeut] += moleFractions_[icat] / fmij;
        }

        if (DEBUG_MODE_ENABLED) {
            for (size_t k = 0; k < m_kk; k++) {
                moleFractionsTmp_[k] = moleFractions_[k];
            }
            for (jNeut = 0; jNeut <  numNeutralMoleculeSpecies_; jNeut++) {
                for (size_t k = 0; k < m_kk; k++) {
                    fmij =  fm_neutralMolec_ions_[k + jNeut * m_kk];
                    moleFractionsTmp_[k] -= fmij * NeutralMolecMoleFractions_[jNeut];
                }
            }
            for (size_t k = 0; k < m_kk; k++) {
                if (fabs(moleFractionsTmp_[k]) > 1.0E-13) {
                    //! Check to see if we have in fact found the inverse.
                    if (anionList_[0] != k) {
                        throw CanteraError("", "neutral molecule calc error");
                    } else {
                        //! For the single anion case, we will allow some slippage
                        if (fabs(moleFractionsTmp_[k]) > 1.0E-5) {
                            throw CanteraError("", "neutral molecule calc error - anion");
                        }
                    }
                }
            }
        }

        // Normalize the Neutral Molecule mole fractions
        sum = 0.0;
        for (size_t k = 0; k < numNeutralMoleculeSpecies_; k++) {
            sum += NeutralMolecMoleFractions_[k];
        }
        for (size_t k = 0; k < numNeutralMoleculeSpecies_; k++) {
            NeutralMolecMoleFractions_[k] /= sum;
        }

        break;

    case  cIonSolnType_SINGLECATION:

        throw CanteraError("eosType", "Unknown type");

        break;

    case  cIonSolnType_MULTICATIONANION:

        throw CanteraError("eosType", "Unknown type");
        break;

    default:

        throw CanteraError("eosType", "Unknown type");
        break;

    }
}

void IonsFromNeutralVPSSTP::getNeutralMoleculeMoleGrads(const doublereal* const dx, doublereal* const dy) const
{
    doublereal sumy, sumdy;

    //check sum dx = 0

    //! Zero the vector we are trying to find.
    for (size_t k = 0; k < numNeutralMoleculeSpecies_; k++) {
        y_[k] = 0.0;
        dy[k] = 0.0;
    }

    switch (ionSolnType_) {

    case cIonSolnType_PASSTHROUGH:
        for (size_t k = 0; k < m_kk; k++) {
            dy[k] = dx[k];
        }
        break;

    case cIonSolnType_SINGLEANION:
        for (size_t k = 0; k < cationList_.size(); k++) {
            //! Get the id for the next cation
            size_t icat = cationList_[k];
            size_t jNeut = fm_invert_ionForNeutral[icat];
            if (jNeut != npos) {
                double fmij =  fm_neutralMolec_ions_[icat + jNeut * m_kk];
                AssertTrace(fmij != 0.0);
                const doublereal temp = 1.0/fmij;
                dy[jNeut] += dx[icat] * temp;
                y_[jNeut] += moleFractions_[icat] * temp;
            }
        }

        for (size_t k = 0; k <  passThroughList_.size(); k++) {
            size_t icat = passThroughList_[k];
            size_t jNeut = fm_invert_ionForNeutral[icat];
            double fmij = fm_neutralMolec_ions_[ icat + jNeut * m_kk];
            const doublereal temp = 1.0/fmij;
            dy[jNeut] += dx[icat] * temp;
            y_[jNeut] += moleFractions_[icat] * temp;
        }
#ifdef DEBUG_MODE_NOT
        //check dy sum to zero
        for (size_t k = 0; k < m_kk; k++) {
            moleFractionsTmp_[k] = dx[k];
        }
        for (jNeut = 0; jNeut <  numNeutralMoleculeSpecies_; jNeut++) {
            for (size_t k = 0; k < m_kk; k++) {
                fmij =  fm_neutralMolec_ions_[k + jNeut * m_kk];
                moleFractionsTmp_[k] -= fmij * dy[jNeut];
            }
        }
        for (size_t k = 0; k < m_kk; k++) {
            if (fabs(moleFractionsTmp_[k]) > 1.0E-13) {
                //! Check to see if we have in fact found the inverse.
                if (anionList_[0] != k) {
                    throw CanteraError("", "neutral molecule calc error");
                } else {
                    //! For the single anion case, we will allow some slippage
                    if (fabs(moleFractionsTmp_[k]) > 1.0E-5) {
                        throw CanteraError("", "neutral molecule calc error - anion");
                    }
                }
            }
        }
#endif
        // Normalize the Neutral Molecule mole fractions
        sumy = 0.0;
        sumdy = 0.0;
        for (size_t k = 0; k < numNeutralMoleculeSpecies_; k++) {
            sumy += y_[k];
            sumdy += dy[k];
        }
        sumy = 1.0 / sumy;
        for (size_t k = 0; k < numNeutralMoleculeSpecies_; k++) {
            dy[k] = dy[k] * sumy - y_[k]*sumdy*sumy*sumy;
        }

        break;

    case  cIonSolnType_SINGLECATION:

        throw CanteraError("eosType", "Unknown type");

        break;

    case  cIonSolnType_MULTICATIONANION:

        throw CanteraError("eosType", "Unknown type");
        break;

    default:

        throw CanteraError("eosType", "Unknown type");
        break;

    }
}

void IonsFromNeutralVPSSTP::setMassFractions(const doublereal* const y)
{
    GibbsExcessVPSSTP::setMassFractions(y);
    calcNeutralMoleculeMoleFractions();
    neutralMoleculePhase_->setMoleFractions(DATA_PTR(NeutralMolecMoleFractions_));
}

void IonsFromNeutralVPSSTP::setMassFractions_NoNorm(const doublereal* const y)
{
    GibbsExcessVPSSTP::setMassFractions_NoNorm(y);
    calcNeutralMoleculeMoleFractions();
    neutralMoleculePhase_->setMoleFractions(DATA_PTR(NeutralMolecMoleFractions_));
}

void IonsFromNeutralVPSSTP::setMoleFractions(const doublereal* const x)
{
    GibbsExcessVPSSTP::setMoleFractions(x);
    calcNeutralMoleculeMoleFractions();
    neutralMoleculePhase_->setMoleFractions(DATA_PTR(NeutralMolecMoleFractions_));
}

void IonsFromNeutralVPSSTP::setMoleFractions_NoNorm(const doublereal* const x)
{
    GibbsExcessVPSSTP::setMoleFractions_NoNorm(x);
    calcNeutralMoleculeMoleFractions();
    neutralMoleculePhase_->setMoleFractions_NoNorm(DATA_PTR(NeutralMolecMoleFractions_));
}

void IonsFromNeutralVPSSTP::setConcentrations(const doublereal* const c)
{
    GibbsExcessVPSSTP::setConcentrations(c);
    calcNeutralMoleculeMoleFractions();
    neutralMoleculePhase_->setMoleFractions(DATA_PTR(NeutralMolecMoleFractions_));
}

/*
 * ------------ Partial Molar Properties of the Solution ------------
 */

void IonsFromNeutralVPSSTP::initThermo()
{
    initLengths();
    GibbsExcessVPSSTP::initThermo();
}

void  IonsFromNeutralVPSSTP::initLengths()
{
    numNeutralMoleculeSpecies_ =  neutralMoleculePhase_->nSpecies();
    moleFractions_.resize(m_kk);
    fm_neutralMolec_ions_.resize(numNeutralMoleculeSpecies_ * m_kk);
    fm_invert_ionForNeutral.resize(m_kk);
    NeutralMolecMoleFractions_.resize(numNeutralMoleculeSpecies_);
    cationList_.resize(m_kk);
    anionList_.resize(m_kk);
    passThroughList_.resize(m_kk);
    moleFractionsTmp_.resize(m_kk);
    muNeutralMolecule_.resize(numNeutralMoleculeSpecies_);
    lnActCoeff_NeutralMolecule_.resize(numNeutralMoleculeSpecies_);
    dlnActCoeffdT_NeutralMolecule_.resize(numNeutralMoleculeSpecies_);
    dlnActCoeffdlnX_diag_NeutralMolecule_.resize(numNeutralMoleculeSpecies_);
    dlnActCoeffdlnN_diag_NeutralMolecule_.resize(numNeutralMoleculeSpecies_);
    dlnActCoeffdlnN_NeutralMolecule_.resize(numNeutralMoleculeSpecies_, numNeutralMoleculeSpecies_, 0.0);

    y_.resize(numNeutralMoleculeSpecies_, 0.0);
    dlnActCoeff_NeutralMolecule_.resize(numNeutralMoleculeSpecies_, 0.0);
    dX_NeutralMolecule_.resize(numNeutralMoleculeSpecies_, 0.0);

}

//!  Return the factor overlap
/*!
 *     @param elnamesVN
 *     @param elemVectorN
 *     @param nElementsN
 *     @param elnamesVI
 *     @param elemVectorI
 *     @param nElementsI
 */
static double factorOverlap(const std::vector<std::string>&  elnamesVN ,
                            const std::vector<double>& elemVectorN,
                            const size_t nElementsN,
                            const std::vector<std::string>&  elnamesVI ,
                            const std::vector<double>& elemVectorI,
                            const size_t nElementsI)
{
    double fMax = 1.0E100;
    for (size_t mi = 0; mi < nElementsI; mi++) {
        if (elnamesVI[mi] != "E") {
            if (elemVectorI[mi] > 1.0E-13) {
                double eiNum = elemVectorI[mi];
                for (size_t mn = 0; mn < nElementsN; mn++) {
                    if (elnamesVI[mi] == elnamesVN[mn]) {
                        if (elemVectorN[mn] <= 1.0E-13) {
                            return 0.0;
                        }
                        fMax = std::min(fMax, elemVectorN[mn]/eiNum);
                    }
                }
            }
        }
    }
    return fMax;
}
void IonsFromNeutralVPSSTP::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    if (id_.size() > 0) {
        if (phaseNode.id() != id_) {
            throw CanteraError("IonsFromNeutralVPSSTP::initThermoXML",
                               "phasenode and Id are incompatible");
        }
    }

    /*
     * Find the Thermo XML node
     */
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError("IonsFromNeutralVPSSTP::initThermoXML",
                           "no thermo XML node");
    }
    XML_Node& thermoNode = phaseNode.child("thermo");



    /*
     * Make sure that the thermo model is IonsFromNeutralMolecule
     */
    string formString = lowercase(thermoNode.attrib("model"));
    if (formString != "ionsfromneutralmolecule") {
        throw CanteraError("IonsFromNeutralVPSSTP::initThermoXML",
                           "model name isn't IonsFromNeutralMolecule: " + formString);
    }

    /*
     * Find the Neutral Molecule Phase
     */
    if (!thermoNode.hasChild("neutralMoleculePhase")) {
        throw CanteraError("IonsFromNeutralVPSSTP::initThermoXML",
                           "no neutralMoleculePhase XML node");
    }
    XML_Node& neutralMoleculeNode = thermoNode.child("neutralMoleculePhase");

    XML_Node* neut_ptr = get_XML_Node(neutralMoleculeNode["datasrc"], 0);
    if (!neut_ptr) {
        throw CanteraError("IonsFromNeutralVPSSTP::initThermoXML",
                           "neut_ptr = 0");
    }

    /*
     *  Create the neutralMolecule ThermoPhase if we haven't already
     */
    if (!neutralMoleculePhase_) {
        neutralMoleculePhase_  = newPhase(*neut_ptr);
    }

    cationList_.clear();
    for (size_t k = 0; k < m_kk; k++) {
        if (charge(k) > 0) {
            cationList_.push_back(k);
        }
    }

    anionList_.clear();
    for (size_t k = 0; k < m_kk; k++) {
        if (charge(k) < 0) {
            anionList_.push_back(k);
        }
    }

    passThroughList_.clear();
    for (size_t k = 0; k < m_kk; k++) {
        if (charge(k) == 0) {
            passThroughList_.push_back(k);
        }
    }

    indexSpecialSpecies_ = npos;
    for (size_t k = 0; k < m_kk; k++) {
        PDSS_IonsFromNeutral* speciesSS =
            dynamic_cast<PDSS_IonsFromNeutral*>(providePDSS(k));
        if (!speciesSS) {
            throw CanteraError("initThermoXML", "Dynamic cast failed");
        }
        if (speciesSS->specialSpecies_ == 1) {
            indexSpecialSpecies_ = k;
        }
        if (speciesSS->specialSpecies_ == 2) {
            indexSecondSpecialSpecies_ = k;
        }
    }


    size_t nElementsN =  neutralMoleculePhase_->nElements();
    const std::vector<std::string>&  elnamesVN = neutralMoleculePhase_->elementNames();
    std::vector<double> elemVectorN(nElementsN);
    std::vector<double> elemVectorN_orig(nElementsN);

    size_t nElementsI =  nElements();
    const std::vector<std::string>&  elnamesVI = elementNames();
    std::vector<double> elemVectorI(nElementsI);

    vector<doublereal> fm_tmp(m_kk);
    for (size_t k = 0; k <  m_kk; k++) {
        fm_invert_ionForNeutral[k] = npos;
    }
    for (size_t jNeut = 0; jNeut <  numNeutralMoleculeSpecies_; jNeut++) {
        for (size_t m = 0; m < nElementsN; m++) {
            elemVectorN[m] = neutralMoleculePhase_->nAtoms(jNeut, m);
        }
        elemVectorN_orig = elemVectorN;
        fm_tmp.assign(m_kk, 0.0);

        for (size_t m = 0; m < nElementsI; m++) {
            elemVectorI[m] = nAtoms(indexSpecialSpecies_, m);
        }
        double fac = factorOverlap(elnamesVN, elemVectorN, nElementsN,
                                   elnamesVI ,elemVectorI, nElementsI);
        if (fac > 0.0) {
            for (size_t m = 0; m < nElementsN; m++) {
                std::string mName = elnamesVN[m];
                for (size_t mi = 0; mi < nElementsI; mi++) {
                    std::string eName = elnamesVI[mi];
                    if (mName == eName) {
                        elemVectorN[m] -= fac * elemVectorI[mi];
                    }

                }
            }
        }
        fm_neutralMolec_ions_[indexSpecialSpecies_  + jNeut * m_kk ] += fac;


        for (size_t k = 0; k < m_kk; k++) {
            for (size_t m = 0; m < nElementsI; m++) {
                elemVectorI[m] = nAtoms(k, m);
            }
            fac = factorOverlap(elnamesVN, elemVectorN, nElementsN,
                                elnamesVI ,elemVectorI, nElementsI);
            if (fac > 0.0) {
                for (size_t m = 0; m < nElementsN; m++) {
                    std::string mName = elnamesVN[m];
                    for (size_t mi = 0; mi < nElementsI; mi++) {
                        std::string eName = elnamesVI[mi];
                        if (mName == eName) {
                            elemVectorN[m] -= fac * elemVectorI[mi];
                        }

                    }
                }
                bool notTaken = true;
                for (size_t iNeut = 0; iNeut < jNeut; iNeut++) {
                    if (fm_invert_ionForNeutral[k] == iNeut) {
                        notTaken = false;
                    }
                }
                if (notTaken) {
                    fm_invert_ionForNeutral[k] = jNeut;
                } else {
                    throw CanteraError("IonsFromNeutralVPSSTP::initThermoXML",
                                       "Simple formula matrix generation failed, one cation is shared between two salts");
                }
            }
            fm_neutralMolec_ions_[k  + jNeut * m_kk] += fac;
        }

        // Ok check the work
        for (size_t m = 0; m < nElementsN; m++) {
            if (fabs(elemVectorN[m]) > 1.0E-13) {
                throw CanteraError("IonsFromNeutralVPSSTP::initThermoXML",
                                   "Simple formula matrix generation failed");
            }
        }


    }
    /*
     * This includes the setStateFromXML calls
     */
    GibbsExcessVPSSTP::initThermoXML(phaseNode, id_);

    /*
     * There is one extra step here. We assure ourselves that we
     * have charge conservation.
     */
}

void IonsFromNeutralVPSSTP::s_update_lnActCoeff() const
{
    size_t icat, jNeut;
    /*
     * Get the activity coefficiens of the neutral molecules
     */
    neutralMoleculePhase_->getLnActivityCoefficients(DATA_PTR(lnActCoeff_NeutralMolecule_));

    switch (ionSolnType_) {
    case cIonSolnType_PASSTHROUGH:
        break;
    case cIonSolnType_SINGLEANION:

        // Do the cation list
        for (size_t k = 0; k < cationList_.size(); k++) {
            //! Get the id for the next cation
            icat = cationList_[k];
            jNeut = fm_invert_ionForNeutral[icat];
            double fmij = fm_neutralMolec_ions_[icat + jNeut * m_kk];
            lnActCoeff_Scaled_[icat] = lnActCoeff_NeutralMolecule_[jNeut] / fmij;
        }

        // Do the anion list
        icat = anionList_[0];
        jNeut = fm_invert_ionForNeutral[icat];
        lnActCoeff_Scaled_[icat]= 0.0;

        // Do the list of neutral molecules
        for (size_t k = 0; k <  passThroughList_.size(); k++) {
            icat = passThroughList_[k];
            jNeut = fm_invert_ionForNeutral[icat];
            lnActCoeff_Scaled_[icat] = lnActCoeff_NeutralMolecule_[jNeut];
        }
        break;

    case cIonSolnType_SINGLECATION:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeff", "Unimplemented type");
        break;
    case cIonSolnType_MULTICATIONANION:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeff", "Unimplemented type");
        break;
    default:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeff", "Unimplemented type");
        break;
    }

}

void IonsFromNeutralVPSSTP::getdlnActCoeffds(const doublereal dTds, const doublereal* const dXds,
        doublereal* dlnActCoeffds) const
{
    size_t icat, jNeut;
    /*
     * Get the activity coefficients of the neutral molecules
     */
    if (!geThermo) {
        for (size_t k = 0; k < m_kk; k++) {
            dlnActCoeffds[k] = dXds[k] / moleFractions_[k];
        }
        return;
    }

    //    static vector_fp dlnActCoeff_NeutralMolecule(numNeutMolSpec);
    //    static vector_fp dX_NeutralMolecule(numNeutMolSpec);


    getNeutralMoleculeMoleGrads(DATA_PTR(dXds),DATA_PTR(dX_NeutralMolecule_));

    // All mole fractions returned to normal

    geThermo->getdlnActCoeffds(dTds, DATA_PTR(dX_NeutralMolecule_), DATA_PTR(dlnActCoeff_NeutralMolecule_));

    switch (ionSolnType_) {
    case cIonSolnType_PASSTHROUGH:
        break;
    case cIonSolnType_SINGLEANION:

        // Do the cation list
        for (size_t k = 0; k < cationList_.size(); k++) {
            //! Get the id for the next cation
            icat = cationList_[k];
            jNeut = fm_invert_ionForNeutral[icat];
            double fmij = fm_neutralMolec_ions_[icat + jNeut * m_kk];
            dlnActCoeffds[icat] = dlnActCoeff_NeutralMolecule_[jNeut]/fmij;
        }

        // Do the anion list
        icat = anionList_[0];
        jNeut = fm_invert_ionForNeutral[icat];
        dlnActCoeffds[icat]= 0.0;

        // Do the list of neutral molecules
        for (size_t k = 0; k < passThroughList_.size(); k++) {
            icat = passThroughList_[k];
            jNeut = fm_invert_ionForNeutral[icat];
            dlnActCoeffds[icat] = dlnActCoeff_NeutralMolecule_[jNeut];
        }
        break;

    case cIonSolnType_SINGLECATION:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeffds", "Unimplemented type");
        break;
    case cIonSolnType_MULTICATIONANION:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeffds", "Unimplemented type");
        break;
    default:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeffds", "Unimplemented type");
        break;
    }

}

void IonsFromNeutralVPSSTP::s_update_dlnActCoeffdT() const
{
    size_t icat, jNeut;
    /*
     * Get the activity coefficients of the neutral molecules
     */
    if (!geThermo) {
        dlnActCoeffdT_Scaled_.assign(m_kk, 0.0);
        return;
    }

    geThermo->getdlnActCoeffdT(DATA_PTR(dlnActCoeffdT_NeutralMolecule_));

    switch (ionSolnType_) {
    case cIonSolnType_PASSTHROUGH:
        break;
    case cIonSolnType_SINGLEANION:

        // Do the cation list
        for (size_t k = 0; k < cationList_.size(); k++) {
            //! Get the id for the next cation
            icat = cationList_[k];
            jNeut = fm_invert_ionForNeutral[icat];
            double fmij = fm_neutralMolec_ions_[icat + jNeut * m_kk];
            dlnActCoeffdT_Scaled_[icat] = dlnActCoeffdT_NeutralMolecule_[jNeut]/fmij;
        }

        // Do the anion list
        icat = anionList_[0];
        jNeut = fm_invert_ionForNeutral[icat];
        dlnActCoeffdT_Scaled_[icat]= 0.0;

        // Do the list of neutral molecules
        for (size_t k = 0; k <  passThroughList_.size(); k++) {
            icat = passThroughList_[k];
            jNeut = fm_invert_ionForNeutral[icat];
            dlnActCoeffdT_Scaled_[icat] = dlnActCoeffdT_NeutralMolecule_[jNeut];
        }
        break;

    case cIonSolnType_SINGLECATION:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeffdT", "Unimplemented type");
        break;
    case cIonSolnType_MULTICATIONANION:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeffdT", "Unimplemented type");
        break;
    default:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeffdT", "Unimplemented type");
        break;
    }

}

void IonsFromNeutralVPSSTP::s_update_dlnActCoeff_dlnX_diag() const
{
    size_t icat, jNeut;
    /*
     * Get the activity coefficients of the neutral molecules
     */
    if (!geThermo) {
        dlnActCoeffdlnX_diag_.assign(m_kk, 0.0);
        return;
    }

    geThermo->getdlnActCoeffdlnX_diag(DATA_PTR(dlnActCoeffdlnX_diag_NeutralMolecule_));

    switch (ionSolnType_) {
    case cIonSolnType_PASSTHROUGH:
        break;
    case cIonSolnType_SINGLEANION:

        // Do the cation list
        for (size_t k = 0; k < cationList_.size(); k++) {
            //! Get the id for the next cation
            icat = cationList_[k];
            jNeut = fm_invert_ionForNeutral[icat];
            double fmij = fm_neutralMolec_ions_[icat + jNeut * m_kk];
            dlnActCoeffdlnX_diag_[icat] = dlnActCoeffdlnX_diag_NeutralMolecule_[jNeut]/fmij;
        }

        // Do the anion list
        icat = anionList_[0];
        jNeut = fm_invert_ionForNeutral[icat];
        dlnActCoeffdlnX_diag_[icat]= 0.0;

        // Do the list of neutral molecules
        for (size_t k = 0; k <  passThroughList_.size(); k++) {
            icat = passThroughList_[k];
            jNeut = fm_invert_ionForNeutral[icat];
            dlnActCoeffdlnX_diag_[icat] = dlnActCoeffdlnX_diag_NeutralMolecule_[jNeut];
        }
        break;

    case cIonSolnType_SINGLECATION:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeff_dlnX_diag()", "Unimplemented type");
        break;
    case cIonSolnType_MULTICATIONANION:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeff_dlnX_diag()", "Unimplemented type");
        break;
    default:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeff_dlnX_diag()", "Unimplemented type");
        break;
    }

}

void IonsFromNeutralVPSSTP::s_update_dlnActCoeff_dlnN_diag() const
{
    size_t icat, jNeut;
    /*
     * Get the activity coefficients of the neutral molecules
     */
    if (!geThermo) {
        dlnActCoeffdlnN_diag_.assign(m_kk, 0.0);
        return;
    }

    geThermo->getdlnActCoeffdlnN_diag(DATA_PTR(dlnActCoeffdlnN_diag_NeutralMolecule_));

    switch (ionSolnType_) {
    case cIonSolnType_PASSTHROUGH:
        break;
    case cIonSolnType_SINGLEANION:

        // Do the cation list
        for (size_t k = 0; k < cationList_.size(); k++) {
            //! Get the id for the next cation
            icat = cationList_[k];
            jNeut = fm_invert_ionForNeutral[icat];
            double fmij = fm_neutralMolec_ions_[icat + jNeut * m_kk];
            dlnActCoeffdlnN_diag_[icat] = dlnActCoeffdlnN_diag_NeutralMolecule_[jNeut]/fmij;
        }

        // Do the anion list
        icat = anionList_[0];
        jNeut = fm_invert_ionForNeutral[icat];
        dlnActCoeffdlnN_diag_[icat]= 0.0;

        // Do the list of neutral molecules
        for (size_t k = 0; k < passThroughList_.size(); k++) {
            icat = passThroughList_[k];
            jNeut = fm_invert_ionForNeutral[icat];
            dlnActCoeffdlnN_diag_[icat] = dlnActCoeffdlnN_diag_NeutralMolecule_[jNeut];
        }
        break;

    case cIonSolnType_SINGLECATION:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeff_dlnN_diag()", "Unimplemented type");
        break;
    case cIonSolnType_MULTICATIONANION:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeff_dlnN_diag()", "Unimplemented type");
        break;
    default:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeff_dlnN_diag()", "Unimplemented type");
        break;
    }

}

void IonsFromNeutralVPSSTP::s_update_dlnActCoeff_dlnN() const
{
    size_t kcat = 0, kNeut = 0, mcat = 0, mNeut = 0;
    doublereal fmij = 0.0;
    dlnActCoeffdlnN_.zero();
    /*
     * Get the activity coefficients of the neutral molecules
     */
    if (!geThermo) {
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_dlnActCoeff_dlnN()", "dynamic cast failed");
    }
    size_t nsp_ge = geThermo->nSpecies();
    geThermo->getdlnActCoeffdlnN(nsp_ge, &(dlnActCoeffdlnN_NeutralMolecule_(0,0)));

    switch (ionSolnType_) {
    case cIonSolnType_PASSTHROUGH:
        break;
    case cIonSolnType_SINGLEANION:

        // Do the cation list
        for (size_t k = 0; k < cationList_.size(); k++) {
            for (size_t m = 0; m < cationList_.size(); m++) {
                kcat = cationList_[k];

                kNeut = fm_invert_ionForNeutral[kcat];
                fmij =  fm_neutralMolec_ions_[kcat + kNeut * m_kk];
                dlnActCoeffdlnN_diag_[kcat] = dlnActCoeffdlnN_diag_NeutralMolecule_[kNeut]/fmij;

                mcat = cationList_[m];
                mNeut = fm_invert_ionForNeutral[mcat];
                double mfmij = fm_neutralMolec_ions_[mcat + mNeut * m_kk];

                dlnActCoeffdlnN_(kcat,mcat) = dlnActCoeffdlnN_NeutralMolecule_(kNeut,mNeut) * mfmij / fmij;

            }
            for (size_t m = 0; m < passThroughList_.size(); m++) {
                mcat = passThroughList_[m];
                mNeut = fm_invert_ionForNeutral[mcat];
                dlnActCoeffdlnN_(kcat, mcat) =  dlnActCoeffdlnN_NeutralMolecule_(kNeut, mNeut) / fmij;
            }
        }

        // Do the anion list -> anion activity coefficient is one
        kcat = anionList_[0];
        kNeut = fm_invert_ionForNeutral[kcat];
        for (size_t k = 0; k < m_kk; k++) {
            dlnActCoeffdlnN_(kcat, k) = 0.0;
            dlnActCoeffdlnN_(k, kcat) = 0.0;
        }

        // Do the list of neutral molecules
        for (size_t k = 0; k <  passThroughList_.size(); k++) {
            kcat = passThroughList_[k];
            kNeut = fm_invert_ionForNeutral[kcat];
            dlnActCoeffdlnN_diag_[kcat] = dlnActCoeffdlnN_diag_NeutralMolecule_[kNeut];

            for (size_t m = 0; m < m_kk; m++) {
                mcat = passThroughList_[m];
                mNeut = fm_invert_ionForNeutral[mcat];
                dlnActCoeffdlnN_(kcat, mcat) =  dlnActCoeffdlnN_NeutralMolecule_(kNeut, mNeut);
            }


            for (size_t m = 0; m < cationList_.size(); m++) {
                mcat = cationList_[m];
                mNeut = fm_invert_ionForNeutral[mcat];
                dlnActCoeffdlnN_(kcat, mcat) =  dlnActCoeffdlnN_NeutralMolecule_(kNeut,mNeut);
            }

        }
        break;

    case cIonSolnType_SINGLECATION:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeff_dlnN", "Unimplemented type");
        break;
    case cIonSolnType_MULTICATIONANION:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeff_dlnN", "Unimplemented type");
        break;
    default:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_lnActCoeff_dlnN", "Unimplemented type");
        break;
    }
}

}
