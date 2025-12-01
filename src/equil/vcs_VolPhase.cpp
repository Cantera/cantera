/**
 * @file vcs_VolPhase.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/equil/vcs_VolPhase.h"
#include "cantera/equil/vcs_solve.h"

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/stringUtils.h"

#include <cstdio>

namespace Cantera
{

vcs_VolPhase::vcs_VolPhase(VCS_SOLVE* owningSolverObject) :
    m_owningSolverObject(owningSolverObject)
{
    if (!m_owningSolverObject) {
        throw CanteraError("vcs_VolPhase::vcs_VolPhase",
                           "owningSolverObject must not be null");
    }
}

void vcs_VolPhase::resize(const size_t phaseNum, const size_t nspecies,
                          const size_t numElem, const char* const phaseName)
{
    AssertThrowMsg(nspecies > 0, "vcs_VolPhase::resize", "nspecies Error");
    m_phi = 0.0;
    m_phiVarIndex = npos;

    if (phaseNum == VP_ID_) {
        if (strcmp(PhaseName.c_str(), phaseName)) {
            throw CanteraError("vcs_VolPhase::resize",
                               "Strings are different: " + PhaseName + " " +
                               phaseName + " :unknown situation");
        }
    } else {
        VP_ID_ = phaseNum;
        if (!phaseName) {
            PhaseName = fmt::format("Phase_{}", VP_ID_);
        } else {
            PhaseName = phaseName;
        }
    }
    if (nspecies > 1) {
        m_singleSpecies = false;
    } else {
        m_singleSpecies = true;
    }

    if (m_numSpecies == nspecies && numElem == m_numElemConstraints) {
        return;
    }

    m_numSpecies = nspecies;
    if (nspecies > 1) {
        m_singleSpecies = false;
    }

    IndSpecies.resize(nspecies, npos);

    Xmol_.resize(nspecies, 0.0);
    creationMoleNumbers_.resize(nspecies, 0.0);
    creationGlobalRxnNumbers_.resize(nspecies, npos);
    for (size_t i = 0; i < nspecies; i++) {
        Xmol_[i] = 1.0/nspecies;
        creationMoleNumbers_[i] = 1.0/nspecies;
        if (IndSpecies[i] >= m_numElemConstraints) {
            creationGlobalRxnNumbers_[i] = IndSpecies[i] - m_numElemConstraints;
        } else {
            creationGlobalRxnNumbers_[i] = npos;
        }
    }

    SS0ChemicalPotential.resize(nspecies, -1.0);
    StarChemicalPotential.resize(nspecies, -1.0);
    StarMolarVol.resize(nspecies, -1.0);
    PartialMolarVol.resize(nspecies, -1.0);
    ActCoeff.resize(nspecies, 1.0);
    np_dLnActCoeffdMolNumber.resize(nspecies, nspecies, 0.0);

    m_speciesUnknownType.resize(nspecies, VCS_SPECIES_TYPE_MOLNUM);
    m_UpToDate = false;
    m_vcsStateStatus = VCS_STATECALC_OLD;
    m_UpToDate_AC = false;
    m_UpToDate_VolStar = false;
    m_UpToDate_VolPM = false;
    m_UpToDate_GStar = false;
    m_UpToDate_G0 = false;

    elemResize(numElem);
}

void vcs_VolPhase::elemResize(const size_t numElemConstraints)
{
    m_elementNames.resize(numElemConstraints);
    m_elementActive.resize(numElemConstraints+1, 1);
    m_elementType.resize(numElemConstraints, VCS_ELEM_TYPE_ABSPOS);
    m_formulaMatrix.resize(m_numSpecies, numElemConstraints, 0.0);
    m_elementNames.resize(numElemConstraints, "");
    m_elemGlobalIndex.resize(numElemConstraints, npos);
    m_numElemConstraints = numElemConstraints;
}

void vcs_VolPhase::_updateActCoeff() const
{
    if (m_isIdealSoln) {
        m_UpToDate_AC = true;
        return;
    }
    TP_ptr->getActivityCoefficients(&ActCoeff[0]);
    m_UpToDate_AC = true;
}

void vcs_VolPhase::_updateG0() const
{
    TP_ptr->getGibbs_ref(&SS0ChemicalPotential[0]);
    m_UpToDate_G0 = true;
}

double vcs_VolPhase::G0_calc_one(size_t kspec) const
{
    if (!m_UpToDate_G0) {
        _updateG0();
    }
    return SS0ChemicalPotential[kspec];
}

void vcs_VolPhase::_updateGStar() const
{
    TP_ptr->getStandardChemPotentials(&StarChemicalPotential[0]);
    m_UpToDate_GStar = true;
}

double vcs_VolPhase::GStar_calc_one(size_t kspec) const
{
    if (!m_UpToDate_GStar) {
        _updateGStar();
    }
    return StarChemicalPotential[kspec];
}

void vcs_VolPhase::setMoleFractions(const double* const xmol)
{
    double sum = -1.0;
    for (size_t k = 0; k < m_numSpecies; k++) {
        Xmol_[k] = xmol[k];
        sum+= xmol[k];
    }
    if (std::fabs(sum) > 1.0E-13) {
        for (size_t k = 0; k < m_numSpecies; k++) {
            Xmol_[k] /= sum;
        }
    }
    _updateMoleFractionDependencies();
    m_UpToDate = false;
    m_vcsStateStatus = VCS_STATECALC_TMP;
}

void vcs_VolPhase::_updateMoleFractionDependencies()
{
    if (TP_ptr) {
        TP_ptr->setMoleFractions(&Xmol_[m_MFStartIndex]);
        TP_ptr->setPressure(Pres_);
    }
    if (!m_isIdealSoln) {
        m_UpToDate_AC = false;
        m_UpToDate_VolPM = false;
    }
}

const vector<double> & vcs_VolPhase::moleFractions() const
{
    return Xmol_;
}

double vcs_VolPhase::moleFraction(size_t k) const
{
    return Xmol_[k];
}

void vcs_VolPhase::setMoleFractionsState(const double totalMoles,
        const double* const moleFractions,
        const int vcsStateStatus)
{
    if (totalMoles != 0.0) {
        // There are other ways to set the mole fractions when VCS_STATECALC
        // is set to a normal settting.
        if (vcsStateStatus != VCS_STATECALC_TMP) {
            throw CanteraError("vcs_VolPhase::setMolesFractionsState",
                               "inappropriate usage");
        }
        m_UpToDate = false;
        m_vcsStateStatus = VCS_STATECALC_TMP;
        if (m_existence == VCS_PHASE_EXIST_ZEROEDPHASE) {
            throw CanteraError("vcs_VolPhase::setMolesFractionsState",
                               "inappropriate usage");
        }
        m_existence = VCS_PHASE_EXIST_YES;
    } else {
        m_UpToDate = true;
        m_vcsStateStatus = vcsStateStatus;
        m_existence = std::min(m_existence, VCS_PHASE_EXIST_NO);
    }
    double fractotal = 1.0;
    v_totalMoles = totalMoles;
    double sum = 0.0;
    for (size_t k = 0; k < m_numSpecies; k++) {
        Xmol_[k] = moleFractions[k];
        sum += moleFractions[k];
    }
    if (sum == 0.0) {
        throw CanteraError("vcs_VolPhase::setMolesFractionsState",
                           "inappropriate usage");
    }
    if (sum  != fractotal) {
        for (size_t k = 0; k < m_numSpecies; k++) {
            Xmol_[k] *= (fractotal /sum);
        }
    }
    _updateMoleFractionDependencies();
}

void vcs_VolPhase::setMolesFromVCS(const int stateCalc,
                                   const double* molesSpeciesVCS)
{
    v_totalMoles = 0.0;

    if (molesSpeciesVCS == 0) {
        if (stateCalc == VCS_STATECALC_OLD) {
            molesSpeciesVCS = &m_owningSolverObject->m_molNumSpecies_old[0];
        } else if (stateCalc == VCS_STATECALC_NEW) {
            molesSpeciesVCS = &m_owningSolverObject->m_molNumSpecies_new[0];
        } else {
            throw CanteraError("vcs_VolPhase::setMolesFromVCS", "shouldn't be here");
        }
    }

    for (size_t k = 0; k < m_numSpecies; k++) {
        if (m_speciesUnknownType[k] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            size_t kglob = IndSpecies[k];
            v_totalMoles += std::max(0.0, molesSpeciesVCS[kglob]);
        }
    }
    if (v_totalMoles > 0.0) {
        for (size_t k = 0; k < m_numSpecies; k++) {
            if (m_speciesUnknownType[k] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                size_t kglob = IndSpecies[k];
                double tmp = std::max(0.0, molesSpeciesVCS[kglob]);
                Xmol_[k] = tmp / v_totalMoles;
            }
        }
        m_existence = VCS_PHASE_EXIST_YES;
    } else {
        // This is where we will start to store a better approximation
        // for the mole fractions, when the phase doesn't exist.
        // This is currently unimplemented.
        m_existence = VCS_PHASE_EXIST_NO;
    }

    // Update the electric potential if it is a solution variable in the
    // equation system
    if (m_phiVarIndex != npos) {
        size_t kglob = IndSpecies[m_phiVarIndex];
        if (m_numSpecies == 1) {
            Xmol_[m_phiVarIndex] = 1.0;
        } else {
            Xmol_[m_phiVarIndex] = 0.0;
        }
        double phi = molesSpeciesVCS[kglob];
        setElectricPotential(phi);
        if (m_numSpecies == 1) {
            m_existence = VCS_PHASE_EXIST_YES;
        }
    }
    _updateMoleFractionDependencies();

    // If stateCalc is old and the total moles is positive, then we have a valid
    // state. If the phase went away, it would be a valid starting point for
    // F_k's. So, save the state.
    if (stateCalc == VCS_STATECALC_OLD && v_totalMoles > 0.0) {
        creationMoleNumbers_ = Xmol_;
    }

    // Set flags indicating we are up to date with the VCS state vector.
    m_UpToDate = true;
    m_vcsStateStatus = stateCalc;
}

void vcs_VolPhase::setMolesFromVCSCheck(const int vcsStateStatus,
                                        const double* molesSpeciesVCS,
                                        const double* const TPhMoles)
{
    setMolesFromVCS(vcsStateStatus, molesSpeciesVCS);

    // Check for consistency with TPhMoles[]
    double Tcheck = TPhMoles[VP_ID_];
    if (Tcheck != v_totalMoles) {
        if (vcs_doubleEqual(Tcheck, v_totalMoles)) {
            Tcheck = v_totalMoles;
        } else {
            throw CanteraError("vcs_VolPhase::setMolesFromVCSCheck",
                  "We have a consistency problem: {} {}", Tcheck, v_totalMoles);
        }
    }
}

void vcs_VolPhase::updateFromVCS_MoleNumbers(const int vcsStateStatus)
{
    if ((!m_UpToDate || vcsStateStatus != m_vcsStateStatus) &&
        (vcsStateStatus == VCS_STATECALC_OLD || vcsStateStatus == VCS_STATECALC_NEW)) {
        setMolesFromVCS(vcsStateStatus);
    }
}

void vcs_VolPhase::sendToVCS_ActCoeff(const int vcsStateStatus,
                                      double* const AC)
{
    updateFromVCS_MoleNumbers(vcsStateStatus);
    if (!m_UpToDate_AC) {
        _updateActCoeff();
    }
    for (size_t k = 0; k < m_numSpecies; k++) {
        size_t kglob = IndSpecies[k];
        AC[kglob] = ActCoeff[k];
    }
}

double vcs_VolPhase::sendToVCS_VolPM(double* const VolPM) const
{
    if (!m_UpToDate_VolPM) {
        _updateVolPM();
    }
    for (size_t k = 0; k < m_numSpecies; k++) {
        size_t kglob = IndSpecies[k];
        VolPM[kglob] = PartialMolarVol[k];
    }
    return m_totalVol;
}

void vcs_VolPhase::sendToVCS_GStar(double* const gstar) const
{
    if (!m_UpToDate_GStar) {
        _updateGStar();
    }
    for (size_t k = 0; k < m_numSpecies; k++) {
        size_t kglob = IndSpecies[k];
        gstar[kglob] = StarChemicalPotential[k];
    }
}

void vcs_VolPhase::setElectricPotential(const double phi)
{
    m_phi = phi;
    TP_ptr->setElectricPotential(m_phi);
    // We have changed the state variable. Set uptodate flags to false
    m_UpToDate_AC = false;
    m_UpToDate_VolStar = false;
    m_UpToDate_VolPM = false;
    m_UpToDate_GStar = false;
}

double vcs_VolPhase::electricPotential() const
{
    return m_phi;
}

void vcs_VolPhase::setState_TP(const double temp, const double pres)
{
    if (Temp_ == temp && Pres_ == pres) {
        return;
    }
    TP_ptr->setElectricPotential(m_phi);
    TP_ptr->setState_TP(temp, pres);
    Temp_ = temp;
    Pres_ = pres;
    m_UpToDate_AC = false;
    m_UpToDate_VolStar = false;
    m_UpToDate_VolPM = false;
    m_UpToDate_GStar = false;
    m_UpToDate_G0 = false;
}

void vcs_VolPhase::_updateVolStar() const
{
    TP_ptr->getStandardVolumes(&StarMolarVol[0]);
    m_UpToDate_VolStar = true;
}

double vcs_VolPhase::_updateVolPM() const
{
    TP_ptr->getPartialMolarVolumes(&PartialMolarVol[0]);
    m_totalVol = 0.0;
    for (size_t k = 0; k < m_numSpecies; k++) {
        m_totalVol += PartialMolarVol[k] * Xmol_[k];
    }
    m_totalVol *= v_totalMoles;
    m_UpToDate_VolPM = true;
    return m_totalVol;
}

void vcs_VolPhase::_updateLnActCoeffJac()
{
    double phaseTotalMoles = v_totalMoles;
    if (phaseTotalMoles < 1.0E-14) {
        phaseTotalMoles = 1.0;
    }

    // Evaluate the current base activity coefficients if necessary
    if (!m_UpToDate_AC) {
        _updateActCoeff();
    }
    if (!TP_ptr) {
        return;
    }
    TP_ptr->getdlnActCoeffdlnN(m_numSpecies, &np_dLnActCoeffdMolNumber(0,0));
    for (size_t j = 0; j < m_numSpecies; j++) {
        double moles_j_base = phaseTotalMoles * Xmol_[j];
        double* const np_lnActCoeffCol = np_dLnActCoeffdMolNumber.ptrColumn(j);
        if (moles_j_base < 1.0E-200) {
            moles_j_base = 1.0E-7 * moles_j_base + 1.0E-13 * phaseTotalMoles + 1.0E-150;
        }
        for (size_t k = 0; k < m_numSpecies; k++) {
            np_lnActCoeffCol[k] = np_lnActCoeffCol[k] * phaseTotalMoles / moles_j_base;
        }
    }

    double deltaMoles_j = 0.0;
    // Make copies of ActCoeff and Xmol_ for use in taking differences
    vector<double> ActCoeff_Base(ActCoeff);
    vector<double> Xmol_Base(Xmol_);
    double TMoles_base = phaseTotalMoles;

    // Loop over the columns species to be deltad
    for (size_t j = 0; j < m_numSpecies; j++) {
        // Calculate a value for the delta moles of species j. Note Xmol_[] and
        // Tmoles are always positive or zero quantities.
        double moles_j_base = phaseTotalMoles * Xmol_Base[j];
        deltaMoles_j = 1.0E-7 * moles_j_base + 1.0E-13 * phaseTotalMoles + 1.0E-150;

        // Now, update the total moles in the phase and all of the mole
        // fractions based on this.
        phaseTotalMoles = TMoles_base + deltaMoles_j;
        for (size_t k = 0; k < m_numSpecies; k++) {
            Xmol_[k] = Xmol_Base[k] * TMoles_base / phaseTotalMoles;
        }
        Xmol_[j] = (moles_j_base + deltaMoles_j) / phaseTotalMoles;

        // Go get new values for the activity coefficients. Note this calls
        // setState_PX();
        _updateMoleFractionDependencies();
        _updateActCoeff();

        // Revert to the base case Xmol_, v_totalMoles
        v_totalMoles = TMoles_base;
        Xmol_ = Xmol_Base;
    }

    // Go get base values for the activity coefficients. Note this calls
    // setState_TPX() again; Just wanted to make sure that cantera is in sync
    // with VolPhase after this call.
    setMoleFractions(&Xmol_Base[0]);
    _updateMoleFractionDependencies();
    _updateActCoeff();
}

void vcs_VolPhase::sendToVCS_LnActCoeffJac(Array2D& np_LnACJac_VCS)
{
    // update the Ln Act Coeff Jacobian entries with respect to the mole number
    // of species in the phase -> we always assume that they are out of date.
    _updateLnActCoeffJac();

    // Now copy over the values
    for (size_t j = 0; j < m_numSpecies; j++) {
        size_t jglob = IndSpecies[j];
        for (size_t k = 0; k < m_numSpecies; k++) {
            size_t kglob = IndSpecies[k];
            np_LnACJac_VCS(kglob,jglob) = np_dLnActCoeffdMolNumber(k,j);
        }
    }
}

void vcs_VolPhase::setPtrThermoPhase(ThermoPhase* tp_ptr)
{
    TP_ptr = tp_ptr;
    Temp_ = TP_ptr->temperature();
    Pres_ = TP_ptr->pressure();
    setState_TP(Temp_, Pres_);
    m_phi = TP_ptr->electricPotential();
    size_t nsp = TP_ptr->nSpecies();
    size_t nelem = TP_ptr->nElements();
    if (nsp != m_numSpecies) {
        if (m_numSpecies != 0) {
            warn_user("vcs_VolPhase::setPtrThermoPhase",
                "Nsp != NVolSpeces: {} {}", nsp, m_numSpecies);
        }
        resize(VP_ID_, nsp, nelem, PhaseName.c_str());
    }
    TP_ptr->getMoleFractions(&Xmol_[0]);
    creationMoleNumbers_ = Xmol_;
    _updateMoleFractionDependencies();

    // figure out ideal solution tag
    if (nsp == 1) {
        m_isIdealSoln = true;
    } else {
        m_isIdealSoln = TP_ptr->isIdeal();
    }
}

double vcs_VolPhase::totalMoles() const
{
    return v_totalMoles;
}

double vcs_VolPhase::molefraction(size_t k) const
{
    return Xmol_[k];
}

void vcs_VolPhase::setCreationMoleNumbers(const double* const n_k,
        const vector<size_t> &creationGlobalRxnNumbers)
{
    creationMoleNumbers_.assign(n_k, n_k+m_numSpecies);
    for (size_t k = 0; k < m_numSpecies; k++) {
        creationGlobalRxnNumbers_[k] = creationGlobalRxnNumbers[k];
    }
}

const vector<double>& vcs_VolPhase::creationMoleNumbers(
        vector<size_t> &creationGlobalRxnNumbers) const
{
    creationGlobalRxnNumbers = creationGlobalRxnNumbers_;
    return creationMoleNumbers_;
}

void vcs_VolPhase::setTotalMoles(const double totalMols)
{
    v_totalMoles = totalMols;
    if (m_singleSpecies && (m_phiVarIndex == 0)) {
        m_existence = VCS_PHASE_EXIST_ALWAYS;
    } else {
        if (totalMols > 0.0) {
            m_existence = VCS_PHASE_EXIST_YES;
        } else {
            m_existence = VCS_PHASE_EXIST_NO;
        }
    }
}

void vcs_VolPhase::setMolesOutOfDate(int stateCalc)
{
    m_UpToDate = false;
    if (stateCalc != -1) {
        m_vcsStateStatus = stateCalc;
    }
}

void vcs_VolPhase::setMolesCurrent(int stateCalc)
{
    m_UpToDate = true;
    m_vcsStateStatus = stateCalc;
}

bool vcs_VolPhase::isIdealSoln() const
{
    return m_isIdealSoln;
}

size_t vcs_VolPhase::phiVarIndex() const
{
    return m_phiVarIndex;
}

void vcs_VolPhase::setPhiVarIndex(size_t phiVarIndex)
{
    m_phiVarIndex = phiVarIndex;
    m_speciesUnknownType[m_phiVarIndex] = VCS_SPECIES_TYPE_INTERFACIALVOLTAGE;
    if (m_singleSpecies && m_phiVarIndex == 0) {
        m_existence = VCS_PHASE_EXIST_ALWAYS;
    }
}

int vcs_VolPhase::exists() const
{
    return m_existence;
}

void vcs_VolPhase::setExistence(const int existence)
{
    if (existence == VCS_PHASE_EXIST_NO || existence == VCS_PHASE_EXIST_ZEROEDPHASE) {
        if (v_totalMoles != 0.0) {
            throw CanteraError("vcs_VolPhase::setExistence",
                               "setting false existence for phase with moles");
        }
    } else if (v_totalMoles == 0.0 && (!m_singleSpecies || m_phiVarIndex != 0)) {
        throw CanteraError("vcs_VolPhase::setExistence",
                "setting true existence for phase with no moles");
    }
    if (m_singleSpecies && m_phiVarIndex == 0 && (existence == VCS_PHASE_EXIST_NO || existence == VCS_PHASE_EXIST_ZEROEDPHASE)) {
        throw CanteraError("vcs_VolPhase::setExistence",
                "Trying to set existence of an electron phase to false");
    }
    m_existence = existence;
}

size_t vcs_VolPhase::spGlobalIndexVCS(const size_t spIndex) const
{
    return IndSpecies[spIndex];
}

void vcs_VolPhase::setSpGlobalIndexVCS(const size_t spIndex,
                                       const size_t spGlobalIndex)
{
    IndSpecies[spIndex] = spGlobalIndex;
    if (spGlobalIndex >= m_numElemConstraints) {
        creationGlobalRxnNumbers_[spIndex] = spGlobalIndex - m_numElemConstraints;
    }
}

size_t vcs_VolPhase::elemGlobalIndex(const size_t e) const
{
    AssertThrow(e < m_numElemConstraints, " vcs_VolPhase::elemGlobalIndex");
    return m_elemGlobalIndex[e];
}

void vcs_VolPhase::setElemGlobalIndex(const size_t eLocal, const size_t eGlobal)
{
    AssertThrow(eLocal < m_numElemConstraints,
                "vcs_VolPhase::setElemGlobalIndex");
    m_elemGlobalIndex[eLocal] = eGlobal;
}

size_t vcs_VolPhase::nElemConstraints() const
{
    return m_numElemConstraints;
}

string vcs_VolPhase::elementName(size_t e) const
{

    if (e < m_elementNames.size()) {
        return m_elementNames[e];
    }
    throw IndexError("vcs_VolPhase::elementName", "element", e, m_elementNames.size());
}

//! This function decides whether a phase has charged species or not.
static bool hasChargedSpecies(const ThermoPhase* const tPhase)
{
    for (size_t k = 0; k < tPhase->nSpecies(); k++) {
        if (tPhase->charge(k) != 0.0) {
            return true;
        }
    }
    return false;
}

//! This utility routine decides whether a Cantera ThermoPhase needs
//! a constraint equation representing the charge neutrality of the
//! phase. It does this by searching for charged species. If it
//! finds one, and if the phase needs one, then it returns true.
static bool chargeNeutralityElement(const ThermoPhase* const tPhase)
{
    int hasCharge = hasChargedSpecies(tPhase);
    if (tPhase->chargeNeutralityNecessary() && hasCharge) {
        return true;
    }
    return false;
}

size_t vcs_VolPhase::transferElementsFM(const ThermoPhase* const tPhase)
{
    size_t nebase = tPhase->nElements();
    size_t ne = nebase;
    size_t ns = tPhase->nSpecies();

    // Decide whether we need an extra element constraint for charge
    // neutrality of the phase
    bool cne = chargeNeutralityElement(tPhase);
    if (cne) {
        ChargeNeutralityElement = ne;
        ne++;
    }

    // Assign and resize structures
    elemResize(ne);

    if (ChargeNeutralityElement != npos) {
        m_elementType[ChargeNeutralityElement] = VCS_ELEM_TYPE_CHARGENEUTRALITY;
    }

    size_t eFound = npos;
    if (hasChargedSpecies(tPhase)) {
        if (cne) {
            // We need a charge neutrality constraint. We also have an Electron
            // Element. These are duplicates of each other. To avoid trouble
            // with possible range error conflicts, sometimes we eliminate the
            // Electron condition. Flag that condition for elimination by
            // toggling the ElActive variable. If we find we need it later, we
            // will retoggle ElActive to true.
            for (size_t eT = 0; eT < nebase; eT++) {
                if (tPhase->elementName(eT) == "E") {
                    eFound = eT;
                    m_elementActive[eT] = 0;
                    m_elementType[eT] = VCS_ELEM_TYPE_ELECTRONCHARGE;
                }
            }
        } else {
            for (size_t eT = 0; eT < nebase; eT++) {
                if (tPhase->elementName(eT) == "E") {
                    eFound = eT;
                    m_elementType[eT] = VCS_ELEM_TYPE_ELECTRONCHARGE;
                }
            }
        }
        if (eFound == npos) {
            eFound = ne;
            m_elementType[ne] = VCS_ELEM_TYPE_ELECTRONCHARGE;
            m_elementActive[ne] = 0;
            string ename = "E";
            m_elementNames[ne] = ename;
            ne++;
            elemResize(ne);
        }
    }

    m_formulaMatrix.resize(ns, ne, 0.0);
    m_speciesUnknownType.resize(ns, VCS_SPECIES_TYPE_MOLNUM);
    elemResize(ne);

    size_t e = 0;
    for (size_t eT = 0; eT < nebase; eT++) {
        m_elementNames[e] = tPhase->elementName(eT);
        m_elementType[e] = tPhase->elementType(eT);
        e++;
    }

    if (cne) {
        string pname = tPhase->name();
        if (pname == "") {
            pname = fmt::format("phase{}", VP_ID_);
        }
        e = ChargeNeutralityElement;
        m_elementNames[e] = "cn_" + pname;
    }

    for (size_t k = 0; k < ns; k++) {
        e = 0;
        for (size_t eT = 0; eT < nebase; eT++) {
            m_formulaMatrix(k,e) = tPhase->nAtoms(k, eT);
            e++;
        }
        if (eFound != npos) {
            m_formulaMatrix(k,eFound) = - tPhase->charge(k);
        }
    }

    if (cne) {
        for (size_t k = 0; k < ns; k++) {
            m_formulaMatrix(k,ChargeNeutralityElement) = tPhase->charge(k);
        }
    }

    // Here, we figure out what is the species types are The logic isn't set in
    // stone, and is just for a particular type of problem that I'm solving
    // first.
    if (ns == 1 && tPhase->charge(0) != 0.0) {
        m_speciesUnknownType[0] = VCS_SPECIES_TYPE_INTERFACIALVOLTAGE;
        setPhiVarIndex(0);
    }

    return ne;
}

int vcs_VolPhase::elementType(const size_t e) const
{
    return m_elementType[e];
}

const Array2D& vcs_VolPhase::getFormulaMatrix() const
{
    return m_formulaMatrix;
}

int vcs_VolPhase::speciesUnknownType(const size_t k) const
{
    return m_speciesUnknownType[k];
}

int vcs_VolPhase::elementActive(const size_t e) const
{
    return m_elementActive[e];
}

size_t vcs_VolPhase::nSpecies() const
{
    return m_numSpecies;
}

string vcs_VolPhase::eos_name() const
{
    return TP_ptr->type();
}
}
