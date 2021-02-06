/**
 *  @file IonsFromNeutralVPSSTP.cpp
 *   Definitions for the object which treats ionic liquids as made of ions as species
 *   even though the thermodynamics is obtained from the neutral molecule representation.
 *  (see \ref thermoprops
 *   and class \link Cantera::IonsFromNeutralVPSSTP IonsFromNeutralVPSSTP\endlink).
 *
 * Header file for a derived class of ThermoPhase that handles variable pressure
 * standard state methods for calculating thermodynamic properties that are
 * further based upon expressions for the excess Gibbs free energy expressed as
 * a function of the mole fractions.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

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
    neutralMoleculePhase_(0),
    geThermo(0)
{
}

IonsFromNeutralVPSSTP::IonsFromNeutralVPSSTP(const std::string& inputFile,
        const std::string& id_) :
    ionSolnType_(cIonSolnType_SINGLEANION),
    numNeutralMoleculeSpecies_(0),
    indexSpecialSpecies_(npos)
{
    initThermoFile(inputFile, id_);
}

IonsFromNeutralVPSSTP::IonsFromNeutralVPSSTP(XML_Node& phaseRoot,
        const std::string& id_) :
    ionSolnType_(cIonSolnType_SINGLEANION),
    numNeutralMoleculeSpecies_(0),
    indexSpecialSpecies_(npos)
{
    importPhase(phaseRoot, this);
}

// ------------ Molar Thermodynamic Properties ----------------------

doublereal IonsFromNeutralVPSSTP::enthalpy_mole() const
{
    getPartialMolarEnthalpies(m_work.data());
    return mean_X(m_work);
}

doublereal IonsFromNeutralVPSSTP::entropy_mole() const
{
    getPartialMolarEntropies(m_work.data());
    return mean_X(m_work);
}

doublereal IonsFromNeutralVPSSTP::gibbs_mole() const
{
    getChemPotentials(m_work.data());
    return mean_X(m_work);
}

doublereal IonsFromNeutralVPSSTP::cp_mole() const
{
    getPartialMolarCp(m_work.data());
    return mean_X(m_work);
}

doublereal IonsFromNeutralVPSSTP::cv_mole() const
{
    // Need to revisit this, as it is wrong
    getPartialMolarCp(m_work.data());
    return mean_X(m_work);
}

// -- Activities, Standard States, Activity Concentrations -----------

void IonsFromNeutralVPSSTP::getDissociationCoeffs(vector_fp& coeffs,
        vector_fp& charges, std::vector<size_t>& neutMolIndex) const
{
    coeffs = fm_neutralMolec_ions_;
    charges = m_speciesCharge;
    neutMolIndex = fm_invert_ionForNeutral;
}

void IonsFromNeutralVPSSTP::getActivityCoefficients(doublereal* ac) const
{
    // Update the activity coefficients
    s_update_lnActCoeff();

    // take the exp of the internally stored coefficients.
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = exp(lnActCoeff_Scaled_[k]);
    }
}

// ---------  Partial Molar Properties of the Solution -------------

void IonsFromNeutralVPSSTP::getChemPotentials(doublereal* mu) const
{
    size_t icat, jNeut;
    doublereal xx, fact2;

    // Get the standard chemical potentials of neutral molecules
    neutralMoleculePhase_->getStandardChemPotentials(muNeutralMolecule_.data());

    switch (ionSolnType_) {
    case cIonSolnType_PASSTHROUGH:
        neutralMoleculePhase_->getChemPotentials(mu);
        break;
    case cIonSolnType_SINGLEANION:
        neutralMoleculePhase_->getLnActivityCoefficients(lnActCoeff_NeutralMolecule_.data());
        fact2 = 2.0 * RT() * log(2.0);

        // Do the cation list
        for (size_t k = 0; k < cationList_.size(); k++) {
            // Get the id for the next cation
            icat = cationList_[k];
            jNeut = fm_invert_ionForNeutral[icat];
            xx = std::max(SmallNumber, moleFractions_[icat]);
            mu[icat] = muNeutralMolecule_[jNeut] + fact2 + RT() * (lnActCoeff_NeutralMolecule_[jNeut] + log(xx));
        }

        // Do the anion list
        icat = anionList_[0];
        jNeut = fm_invert_ionForNeutral[icat];
        xx = std::max(SmallNumber, moleFractions_[icat]);
        mu[icat] = RT() * log(xx);

        // Do the list of neutral molecules
        for (size_t k = 0; k < passThroughList_.size(); k++) {
            icat = passThroughList_[k];
            jNeut = fm_invert_ionForNeutral[icat];
            xx = std::max(SmallNumber, moleFractions_[icat]);
            mu[icat] = muNeutralMolecule_[jNeut] + RT() * (lnActCoeff_NeutralMolecule_[jNeut] + log(xx));
        }
        break;

    case cIonSolnType_SINGLECATION:
        throw CanteraError("IonsFromNeutralVPSSTP::getChemPotentials", "Unknown type");
    case cIonSolnType_MULTICATIONANION:
        throw CanteraError("IonsFromNeutralVPSSTP::getChemPotentials", "Unknown type");
    default:
        throw CanteraError("IonsFromNeutralVPSSTP::getChemPotentials", "Unknown type");
    }
}

void IonsFromNeutralVPSSTP::getPartialMolarEnthalpies(doublereal* hbar) const
{
    // Get the nondimensional standard state enthalpies
    getEnthalpy_RT(hbar);

    // dimensionalize it.
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] *= RT();
    }

    // Update the activity coefficients, This also update the internally stored
    // molalities.
    s_update_lnActCoeff();
    s_update_dlnActCoeffdT();
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] -= RT() * temperature() * dlnActCoeffdT_Scaled_[k];
    }
}

void IonsFromNeutralVPSSTP::getPartialMolarEntropies(doublereal* sbar) const
{
    // Get the nondimensional standard state entropies
    getEntropy_R(sbar);

    // Update the activity coefficients, This also update the internally stored
    // molalities.
    s_update_lnActCoeff();
    s_update_dlnActCoeffdT();

    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(moleFractions_[k], SmallNumber);
        sbar[k] += - lnActCoeff_Scaled_[k] -log(xx) - temperature() * dlnActCoeffdT_Scaled_[k];
    }

    // dimensionalize it.
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

void IonsFromNeutralVPSSTP::calcDensity()
{
    // This is a two phase process. First, we calculate the standard states
    // within the neutral molecule phase.
    neutralMoleculePhase_->setState_TP(temperature(), pressure());

    // Calculate the partial molar volumes, and then the density of the fluid
    Phase::assignDensity(neutralMoleculePhase_->density());
}

void IonsFromNeutralVPSSTP::calcIonMoleFractions(doublereal* const mf) const
{
    // Download the neutral mole fraction vector into the vector,
    // NeutralMolecMoleFractions_[]
    neutralMoleculePhase_->getMoleFractions(NeutralMolecMoleFractions_.data());

    // Zero the mole fractions
    for (size_t k = 0; k < m_kk; k++) {
        mf[k] = 0.0;
    }

    // Use the formula matrix to calculate the relative mole numbers.
    for (size_t jNeut = 0; jNeut < numNeutralMoleculeSpecies_; jNeut++) {
        for (size_t k = 0; k < m_kk; k++) {
            double fmij = fm_neutralMolec_ions_[k + jNeut * m_kk];
            mf[k] += fmij * NeutralMolecMoleFractions_[jNeut];
        }
    }

    // Normalize the new mole fractions
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

    // Zero the vector we are trying to find.
    for (size_t k = 0; k < numNeutralMoleculeSpecies_; k++) {
        NeutralMolecMoleFractions_[k] = 0.0;
    }
    sum = -1.0;
    for (size_t k = 0; k < m_kk; k++) {
        sum += moleFractions_[k];
    }
    if (fabs(sum) > 1.0E-11) {
        throw CanteraError("IonsFromNeutralVPSSTP::calcNeutralMoleculeMoleFractions",
                           "molefracts don't sum to one: {}", sum);
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
            // Get the id for the next cation
            icat = cationList_[k];
            jNeut = fm_invert_ionForNeutral[icat];
            if (jNeut != npos) {
                fmij = fm_neutralMolec_ions_[icat + jNeut * m_kk];
                AssertTrace(fmij != 0.0);
                NeutralMolecMoleFractions_[jNeut] += moleFractions_[icat] / fmij;
            }
        }

        for (size_t k = 0; k < passThroughList_.size(); k++) {
            icat = passThroughList_[k];
            jNeut = fm_invert_ionForNeutral[icat];
            fmij = fm_neutralMolec_ions_[ icat + jNeut * m_kk];
            NeutralMolecMoleFractions_[jNeut] += moleFractions_[icat] / fmij;
        }

        for (size_t k = 0; k < m_kk; k++) {
            moleFractionsTmp_[k] = moleFractions_[k];
        }
        for (jNeut = 0; jNeut < numNeutralMoleculeSpecies_; jNeut++) {
            for (size_t k = 0; k < m_kk; k++) {
                fmij = fm_neutralMolec_ions_[k + jNeut * m_kk];
                moleFractionsTmp_[k] -= fmij * NeutralMolecMoleFractions_[jNeut];
            }
        }
        for (size_t k = 0; k < m_kk; k++) {
            if (fabs(moleFractionsTmp_[k]) > 1.0E-13) {
                // Check to see if we have in fact found the inverse.
                if (anionList_[0] != k) {
                    throw CanteraError("IonsFromNeutralVPSSTP::calcNeutralMoleculeMoleFractions",
                                       "neutral molecule calc error");
                } else {
                    // For the single anion case, we will allow some slippage
                    if (fabs(moleFractionsTmp_[k]) > 1.0E-5) {
                        throw CanteraError("IonsFromNeutralVPSSTP::calcNeutralMoleculeMoleFractions",
                                           "neutral molecule calc error - anion");
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

    case cIonSolnType_SINGLECATION:
        throw CanteraError("IonsFromNeutralVPSSTP::calcNeutralMoleculeMoleFractions", "Unknown type");
        break;
    case cIonSolnType_MULTICATIONANION:
        throw CanteraError("IonsFromNeutralVPSSTP::calcNeutralMoleculeMoleFractions", "Unknown type");
        break;
    default:
        throw CanteraError("IonsFromNeutralVPSSTP::calcNeutralMoleculeMoleFractions", "Unknown type");
        break;
    }
}

void IonsFromNeutralVPSSTP::getNeutralMoleculeMoleGrads(const doublereal* const dx, doublereal* const dy) const
{
    doublereal sumy, sumdy;

    // check sum dx = 0
    // Zero the vector we are trying to find.
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
            // Get the id for the next cation
            size_t icat = cationList_[k];
            size_t jNeut = fm_invert_ionForNeutral[icat];
            if (jNeut != npos) {
                double fmij = fm_neutralMolec_ions_[icat + jNeut * m_kk];
                AssertTrace(fmij != 0.0);
                const doublereal temp = 1.0/fmij;
                dy[jNeut] += dx[icat] * temp;
                y_[jNeut] += moleFractions_[icat] * temp;
            }
        }

        for (size_t k = 0; k < passThroughList_.size(); k++) {
            size_t icat = passThroughList_[k];
            size_t jNeut = fm_invert_ionForNeutral[icat];
            double fmij = fm_neutralMolec_ions_[ icat + jNeut * m_kk];
            const doublereal temp = 1.0/fmij;
            dy[jNeut] += dx[icat] * temp;
            y_[jNeut] += moleFractions_[icat] * temp;
        }
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

    case cIonSolnType_SINGLECATION:
        throw CanteraError("IonsFromNeutralVPSSTP::getNeutralMoleculeMoleGrads",
                           "Unknown type");
        break;
    case cIonSolnType_MULTICATIONANION:
        throw CanteraError("IonsFromNeutralVPSSTP::getNeutralMoleculeMoleGrads",
                           "Unknown type");
        break;
    default:
        throw CanteraError("IonsFromNeutralVPSSTP::getNeutralMoleculeMoleGrads",
                           "Unknown type");
        break;
    }
}

void IonsFromNeutralVPSSTP::compositionChanged()
{
    GibbsExcessVPSSTP::compositionChanged();
    calcNeutralMoleculeMoleFractions();
    neutralMoleculePhase_->setMoleFractions(NeutralMolecMoleFractions_.data());
}

// ------------ Partial Molar Properties of the Solution ------------

//!  Return the factor overlap
/*!
 *     @param elnamesVN
 *     @param elemVectorN
 *     @param nElementsN
 *     @param elnamesVI
 *     @param elemVectorI
 *     @param nElementsI
 */
static double factorOverlap(const std::vector<std::string>& elnamesVN ,
                            const vector_fp& elemVectorN,
                            const size_t nElementsN,
                            const std::vector<std::string>& elnamesVI ,
                            const vector_fp& elemVectorI,
                            const size_t nElementsI)
{
    double fMax = 1.0E100;
    for (size_t mi = 0; mi < nElementsI; mi++) {
        if (elnamesVI[mi] != "E" && elemVectorI[mi] > 1.0E-13) {
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
    return fMax;
}

void IonsFromNeutralVPSSTP::setParameters(const AnyMap& phaseNode,
                                          const AnyMap& rootNode)
{
    ThermoPhase::setParameters(phaseNode, rootNode);
    m_rootNode = rootNode;
}

void IonsFromNeutralVPSSTP::initThermo()
{
    if (m_input.hasKey("neutral-phase")) {
        string neutralName = m_input["neutral-phase"].asString();
        const auto& slash = boost::ifind_last(neutralName, "/");
        if (slash) {
            string fileName(neutralName.begin(), slash.begin());
            neutralName = string(slash.end(), neutralName.end());
            AnyMap infile = AnyMap::fromYamlFile(fileName,
                        m_input.getString("__file__", ""));
            AnyMap& phaseNode = infile["phases"].getMapWhere("name", neutralName);
            setNeutralMoleculePhase(newPhase(phaseNode, infile));
        } else {
            AnyMap& phaseNode = m_rootNode["phases"].getMapWhere("name", neutralName);
            setNeutralMoleculePhase(newPhase(phaseNode, m_rootNode));
        }
    }

    if (!neutralMoleculePhase_) {
        throw CanteraError(
            "IonsFromNeutralVPSSTP::initThermo",
            "The neutral phase has not been initialized. Are you missing the "
            "'neutral-phase' key?"
        );
    }

    size_t nElementsN = neutralMoleculePhase_->nElements();
    const std::vector<std::string>& elnamesVN = neutralMoleculePhase_->elementNames();
    vector_fp elemVectorN(nElementsN);

    size_t nElementsI = nElements();
    const std::vector<std::string>& elnamesVI = elementNames();
    vector_fp elemVectorI(nElementsI);

    if (indexSpecialSpecies_ == npos) {
        throw CanteraError(
            "IonsFromNeutralVPSSTP::initThermo",
            "No special-species were specified in the phase."
        );
    }
    for (size_t m = 0; m < nElementsI; m++) {
        elemVectorI[m] = nAtoms(indexSpecialSpecies_, m);
    }

    for (size_t jNeut = 0; jNeut < numNeutralMoleculeSpecies_; jNeut++) {
        for (size_t m = 0; m < nElementsN; m++) {
            elemVectorN[m] = neutralMoleculePhase_->nAtoms(jNeut, m);
        }

        double fac = factorOverlap(elnamesVN, elemVectorN, nElementsN,
                                   elnamesVI ,elemVectorI, nElementsI);
        if (fac > 0.0) {
            for (size_t m = 0; m < nElementsN; m++) {
                for (size_t mi = 0; mi < nElementsI; mi++) {
                    if (elnamesVN[m] == elnamesVI[mi]) {
                        elemVectorN[m] -= fac * elemVectorI[mi];
                    }

                }
            }
        }
        fm_neutralMolec_ions_[indexSpecialSpecies_ + jNeut * m_kk ] += fac;

        for (size_t k = 0; k < m_kk; k++) {
            for (size_t m = 0; m < nElementsI; m++) {
                elemVectorI[m] = nAtoms(k, m);
            }
            fac = factorOverlap(elnamesVN, elemVectorN, nElementsN,
                                elnamesVI ,elemVectorI, nElementsI);
            if (fac > 0.0) {
                for (size_t m = 0; m < nElementsN; m++) {
                    for (size_t mi = 0; mi < nElementsI; mi++) {
                        if (elnamesVN[m] ==  elnamesVI[mi]) {
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
                    throw CanteraError("IonsFromNeutralVPSSTP::initThermo",
                                       "Simple formula matrix generation failed, one cation is shared between two salts");
                }
            }
            fm_neutralMolec_ions_[k + jNeut * m_kk] += fac;
        }

        // Ok check the work
        for (size_t m = 0; m < nElementsN; m++) {
            if (fabs(elemVectorN[m]) > 1.0E-13) {
                throw CanteraError("IonsFromNeutralVPSSTP::initThermo",
                                   "Simple formula matrix generation failed");
            }
        }
    }

    GibbsExcessVPSSTP::initThermo();
}

void IonsFromNeutralVPSSTP::setNeutralMoleculePhase(shared_ptr<ThermoPhase> neutral)
{
    neutralMoleculePhase_ = neutral;
    geThermo = dynamic_cast<GibbsExcessVPSSTP*>(neutralMoleculePhase_.get());
    numNeutralMoleculeSpecies_ = neutralMoleculePhase_->nSpecies();
    fm_neutralMolec_ions_.resize(numNeutralMoleculeSpecies_ * m_kk);
    NeutralMolecMoleFractions_.resize(numNeutralMoleculeSpecies_);
    muNeutralMolecule_.resize(numNeutralMoleculeSpecies_);
    lnActCoeff_NeutralMolecule_.resize(numNeutralMoleculeSpecies_);
    dlnActCoeffdT_NeutralMolecule_.resize(numNeutralMoleculeSpecies_);
    dlnActCoeffdlnX_diag_NeutralMolecule_.resize(numNeutralMoleculeSpecies_);
    dlnActCoeffdlnN_diag_NeutralMolecule_.resize(numNeutralMoleculeSpecies_);
    dlnActCoeffdlnN_NeutralMolecule_.resize(numNeutralMoleculeSpecies_, numNeutralMoleculeSpecies_, 0.0);
    y_.resize(numNeutralMoleculeSpecies_, 0.0);
    dlnActCoeff_NeutralMolecule_.resize(numNeutralMoleculeSpecies_, 0.0);
    dX_NeutralMolecule_.resize(numNeutralMoleculeSpecies_, 0.0);
    for (size_t k = 0; k < nSpecies(); k++) {
        providePDSS(k)->setParent(this, k);
    }
}

shared_ptr<ThermoPhase> IonsFromNeutralVPSSTP::getNeutralMoleculePhase()
{
    return neutralMoleculePhase_;
}

bool IonsFromNeutralVPSSTP::addSpecies(shared_ptr<Species> spec)
{
    bool added = GibbsExcessVPSSTP::addSpecies(spec);
    if (added) {
        moleFractions_.push_back(0.0);
        moleFractionsTmp_.push_back(0.0);
        m_work.push_back(0.0);
        fm_invert_ionForNeutral.push_back(npos);
        fm_neutralMolec_ions_.resize(numNeutralMoleculeSpecies_ * m_kk);

        if (spec->charge > 0) {
            cationList_.push_back(m_kk-1);
        } else if (spec->charge < 0) {
            anionList_.push_back(m_kk-1);
        } else {
            passThroughList_.push_back(m_kk-1);
        }

        if (spec->input.hasKey("equation-of-state")) {
            auto& ss = spec->input["equation-of-state"].getMapWhere(
                "model", "ions-from-neutral-molecule");
            if (ss.getBool("special-species", false)) {
                indexSpecialSpecies_ = m_kk - 1;
            }
        }
    }
    return added;
}

void IonsFromNeutralVPSSTP::setParametersFromXML(const XML_Node& thermoNode)
{
    GibbsExcessVPSSTP::setParametersFromXML(thermoNode);
    // Find the Neutral Molecule Phase
    if (!thermoNode.hasChild("neutralMoleculePhase")) {
        throw CanteraError("IonsFromNeutralVPSSTP::setParametersFromXML",
                           "no neutralMoleculePhase XML node");
    }
    XML_Node& neutralMoleculeNode = thermoNode.child("neutralMoleculePhase");

    XML_Node* neut_ptr = get_XML_Node(neutralMoleculeNode["datasrc"], 0);
    if (!neut_ptr) {
        throw CanteraError("IonsFromNeutralVPSSTP::setParametersFromXML",
                           "neut_ptr = 0");
    }

    setNeutralMoleculePhase(shared_ptr<ThermoPhase>(newPhase(*neut_ptr)));
}

void IonsFromNeutralVPSSTP::s_update_lnActCoeff() const
{
    size_t icat, jNeut;
    // Get the activity coefficients of the neutral molecules
    neutralMoleculePhase_->getLnActivityCoefficients(lnActCoeff_NeutralMolecule_.data());

    switch (ionSolnType_) {
    case cIonSolnType_PASSTHROUGH:
        break;
    case cIonSolnType_SINGLEANION:
        // Do the cation list
        for (size_t k = 0; k < cationList_.size(); k++) {
            // Get the id for the next cation
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
        for (size_t k = 0; k < passThroughList_.size(); k++) {
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
    // Get the activity coefficients of the neutral molecules
    if (!geThermo) {
        for (size_t k = 0; k < m_kk; k++) {
            dlnActCoeffds[k] = dXds[k] / moleFractions_[k];
        }
        return;
    }

    getNeutralMoleculeMoleGrads(dXds, dX_NeutralMolecule_.data());

    // All mole fractions returned to normal
    geThermo->getdlnActCoeffds(dTds, dX_NeutralMolecule_.data(), dlnActCoeff_NeutralMolecule_.data());

    switch (ionSolnType_) {
    case cIonSolnType_PASSTHROUGH:
        break;
    case cIonSolnType_SINGLEANION:
        // Do the cation list
        for (size_t k = 0; k < cationList_.size(); k++) {
            // Get the id for the next cation
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
        throw CanteraError("IonsFromNeutralVPSSTP::getdlnActCoeffds", "Unimplemented type");
        break;
    case cIonSolnType_MULTICATIONANION:
        throw CanteraError("IonsFromNeutralVPSSTP::getdlnActCoeffds", "Unimplemented type");
        break;
    default:
        throw CanteraError("IonsFromNeutralVPSSTP::getdlnActCoeffds", "Unimplemented type");
        break;
    }
}

void IonsFromNeutralVPSSTP::s_update_dlnActCoeffdT() const
{
    size_t icat, jNeut;

    // Get the activity coefficients of the neutral molecules
    if (!geThermo) {
        dlnActCoeffdT_Scaled_.assign(m_kk, 0.0);
        return;
    }

    geThermo->getdlnActCoeffdT(dlnActCoeffdT_NeutralMolecule_.data());

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
        for (size_t k = 0; k < passThroughList_.size(); k++) {
            icat = passThroughList_[k];
            jNeut = fm_invert_ionForNeutral[icat];
            dlnActCoeffdT_Scaled_[icat] = dlnActCoeffdT_NeutralMolecule_[jNeut];
        }
        break;

    case cIonSolnType_SINGLECATION:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_dlnActCoeffdT", "Unimplemented type");
        break;
    case cIonSolnType_MULTICATIONANION:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_dlnActCoeffdT", "Unimplemented type");
        break;
    default:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_dlnActCoeffdT", "Unimplemented type");
        break;
    }
}

void IonsFromNeutralVPSSTP::s_update_dlnActCoeff_dlnX_diag() const
{
    size_t icat, jNeut;

    // Get the activity coefficients of the neutral molecules
    if (!geThermo) {
        dlnActCoeffdlnX_diag_.assign(m_kk, 0.0);
        return;
    }

    geThermo->getdlnActCoeffdlnX_diag(dlnActCoeffdlnX_diag_NeutralMolecule_.data());

    switch (ionSolnType_) {
    case cIonSolnType_PASSTHROUGH:
        break;
    case cIonSolnType_SINGLEANION:
        // Do the cation list
        for (size_t k = 0; k < cationList_.size(); k++) {
            // Get the id for the next cation
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
        for (size_t k = 0; k < passThroughList_.size(); k++) {
            icat = passThroughList_[k];
            jNeut = fm_invert_ionForNeutral[icat];
            dlnActCoeffdlnX_diag_[icat] = dlnActCoeffdlnX_diag_NeutralMolecule_[jNeut];
        }
        break;

    case cIonSolnType_SINGLECATION:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_dlnActCoeff_dlnX_diag", "Unimplemented type");
        break;
    case cIonSolnType_MULTICATIONANION:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_dlnActCoeff_dlnX_diag", "Unimplemented type");
        break;
    default:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_dlnActCoeff_dlnX_diag", "Unimplemented type");
        break;
    }
}

void IonsFromNeutralVPSSTP::s_update_dlnActCoeff_dlnN_diag() const
{
    size_t icat, jNeut;

    // Get the activity coefficients of the neutral molecules
    if (!geThermo) {
        dlnActCoeffdlnN_diag_.assign(m_kk, 0.0);
        return;
    }

    geThermo->getdlnActCoeffdlnN_diag(dlnActCoeffdlnN_diag_NeutralMolecule_.data());

    switch (ionSolnType_) {
    case cIonSolnType_PASSTHROUGH:
        break;
    case cIonSolnType_SINGLEANION:
        // Do the cation list
        for (size_t k = 0; k < cationList_.size(); k++) {
            // Get the id for the next cation
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
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_dlnActCoeff_dlnN_diag", "Unimplemented type");
        break;
    case cIonSolnType_MULTICATIONANION:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_dlnActCoeff_dlnN_diag", "Unimplemented type");
        break;
    default:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_dlnActCoeff_dlnN_diag", "Unimplemented type");
        break;
    }
}

void IonsFromNeutralVPSSTP::s_update_dlnActCoeff_dlnN() const
{
    size_t kcat = 0, kNeut = 0, mcat = 0, mNeut = 0;
    doublereal fmij = 0.0;
    dlnActCoeffdlnN_.zero();
    // Get the activity coefficients of the neutral molecules
    if (!geThermo) {
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_dlnActCoeff_dlnN", "dynamic cast failed");
    }
    size_t nsp_ge = geThermo->nSpecies();
    geThermo->getdlnActCoeffdlnN(nsp_ge, &dlnActCoeffdlnN_NeutralMolecule_(0,0));

    switch (ionSolnType_) {
    case cIonSolnType_PASSTHROUGH:
        break;
    case cIonSolnType_SINGLEANION:
        // Do the cation list
        for (size_t k = 0; k < cationList_.size(); k++) {
            for (size_t m = 0; m < cationList_.size(); m++) {
                kcat = cationList_[k];

                kNeut = fm_invert_ionForNeutral[kcat];
                fmij = fm_neutralMolec_ions_[kcat + kNeut * m_kk];
                dlnActCoeffdlnN_diag_[kcat] = dlnActCoeffdlnN_diag_NeutralMolecule_[kNeut]/fmij;

                mcat = cationList_[m];
                mNeut = fm_invert_ionForNeutral[mcat];
                double mfmij = fm_neutralMolec_ions_[mcat + mNeut * m_kk];

                dlnActCoeffdlnN_(kcat,mcat) = dlnActCoeffdlnN_NeutralMolecule_(kNeut,mNeut) * mfmij / fmij;

            }
            for (size_t m = 0; m < passThroughList_.size(); m++) {
                mcat = passThroughList_[m];
                mNeut = fm_invert_ionForNeutral[mcat];
                dlnActCoeffdlnN_(kcat, mcat) = dlnActCoeffdlnN_NeutralMolecule_(kNeut, mNeut) / fmij;
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
        for (size_t k = 0; k < passThroughList_.size(); k++) {
            kcat = passThroughList_[k];
            kNeut = fm_invert_ionForNeutral[kcat];
            dlnActCoeffdlnN_diag_[kcat] = dlnActCoeffdlnN_diag_NeutralMolecule_[kNeut];

            for (size_t m = 0; m < m_kk; m++) {
                mcat = passThroughList_[m];
                mNeut = fm_invert_ionForNeutral[mcat];
                dlnActCoeffdlnN_(kcat, mcat) = dlnActCoeffdlnN_NeutralMolecule_(kNeut, mNeut);
            }

            for (size_t m = 0; m < cationList_.size(); m++) {
                mcat = cationList_[m];
                mNeut = fm_invert_ionForNeutral[mcat];
                dlnActCoeffdlnN_(kcat, mcat) = dlnActCoeffdlnN_NeutralMolecule_(kNeut,mNeut);
            }
        }
        break;

    case cIonSolnType_SINGLECATION:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_dlnActCoeff_dlnN", "Unimplemented type");
        break;
    case cIonSolnType_MULTICATIONANION:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_dlnActCoeff_dlnN", "Unimplemented type");
        break;
    default:
        throw CanteraError("IonsFromNeutralVPSSTP::s_update_dlnActCoeff_dlnN", "Unimplemented type");
        break;
    }
}

}
