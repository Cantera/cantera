/**
 * @file vcs_VolPhase.cpp
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#include "cantera/equil/vcs_VolPhase.h"
#include "cantera/equil/vcs_internal.h"
#include "cantera/equil/vcs_SpeciesProperties.h"
#include "cantera/equil/vcs_species_thermo.h"
#include "cantera/equil/vcs_solve.h"

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/mix_defs.h"

#include <sstream>
#include <cstdio>

namespace VCSnonideal
{

vcs_VolPhase::vcs_VolPhase(VCS_SOLVE* owningSolverObject) :
    m_owningSolverObject(0),
    VP_ID_(npos),
    m_singleSpecies(true),
    m_gasPhase(false),
    m_eqnState(VCS_EOS_CONSTANT),
    ChargeNeutralityElement(npos),
    p_VCS_UnitsFormat(VCS_UNITS_MKS),
    p_activityConvention(0),
    m_numElemConstraints(0),
    m_elemGlobalIndex(0),
    m_numSpecies(0),
    m_totalMolesInert(0.0),
    m_isIdealSoln(false),
    m_existence(VCS_PHASE_EXIST_NO),
    m_MFStartIndex(0),
    IndSpecies(0),
    m_useCanteraCalls(false),
    TP_ptr(0),
    v_totalMoles(0.0),
    creationMoleNumbers_(0),
    creationGlobalRxnNumbers_(0),
    m_phiVarIndex(npos),
    m_totalVol(0.0),
    m_vcsStateStatus(VCS_STATECALC_OLD),
    m_phi(0.0),
    m_UpToDate(false),
    m_UpToDate_AC(false),
    m_UpToDate_VolStar(false),
    m_UpToDate_VolPM(false),
    m_UpToDate_GStar(false),
    m_UpToDate_G0(false),
    Temp_(273.15),
    Pres_(1.01325E5)
{
    m_owningSolverObject = owningSolverObject;
}

vcs_VolPhase::~vcs_VolPhase()
{
    for (size_t k = 0; k < m_numSpecies; k++) {
        vcs_SpeciesProperties* sp = ListSpeciesPtr[k];
        delete sp;
        sp = 0;
    }
}

vcs_VolPhase::vcs_VolPhase(const vcs_VolPhase& b) :
    m_owningSolverObject(b.m_owningSolverObject),
    VP_ID_(b.VP_ID_),
    m_singleSpecies(b.m_singleSpecies),
    m_gasPhase(b.m_gasPhase),
    m_eqnState(b.m_eqnState),
    ChargeNeutralityElement(b.ChargeNeutralityElement),
    p_VCS_UnitsFormat(b.p_VCS_UnitsFormat),
    p_activityConvention(b.p_activityConvention),
    m_numElemConstraints(b.m_numElemConstraints),
    m_numSpecies(b.m_numSpecies),
    m_totalMolesInert(b.m_totalMolesInert),
    m_isIdealSoln(b.m_isIdealSoln),
    m_existence(b.m_existence),
    m_MFStartIndex(b.m_MFStartIndex),
    m_useCanteraCalls(b.m_useCanteraCalls),
    TP_ptr(b.TP_ptr),
    v_totalMoles(b.v_totalMoles),
    creationMoleNumbers_(0),
    creationGlobalRxnNumbers_(0),
    m_phiVarIndex(npos),
    m_totalVol(b.m_totalVol),
    m_vcsStateStatus(VCS_STATECALC_OLD),
    m_phi(b.m_phi),
    m_UpToDate(false),
    m_UpToDate_AC(false),
    m_UpToDate_VolStar(false),
    m_UpToDate_VolPM(false),
    m_UpToDate_GStar(false),
    m_UpToDate_G0(false),
    Temp_(b.Temp_),
    Pres_(b.Pres_)
{
    //! Objects that are owned by this object are deep copied here, except for
    //! the ThermoPhase object. The assignment operator does most of the work.
    *this = b;
}

vcs_VolPhase& vcs_VolPhase::operator=(const vcs_VolPhase& b)
{
    if (&b != this) {
        size_t old_num = m_numSpecies;

        //  Note: we comment this out for the assignment operator
        //        specifically, because it isn't true for the assignment
        //        operator but is true for a copy constructor
        // m_owningSolverObject = b.m_owningSolverObject;

        VP_ID_               = b.VP_ID_;
        m_singleSpecies     = b.m_singleSpecies;
        m_gasPhase            = b.m_gasPhase;
        m_eqnState            = b.m_eqnState;
        ChargeNeutralityElement = b.ChargeNeutralityElement;
        p_VCS_UnitsFormat   = b.p_VCS_UnitsFormat;
        p_activityConvention= b.p_activityConvention;
        m_numSpecies = b.m_numSpecies;
        m_numElemConstraints    = b.m_numElemConstraints;
        m_elementNames.resize(b.m_numElemConstraints);
        for (size_t e = 0; e < b.m_numElemConstraints; e++) {
            m_elementNames[e] = b.m_elementNames[e];
        }
        m_elementActive = b.m_elementActive;
        m_elementType = b.m_elementType;
        m_formulaMatrix.resize(m_numElemConstraints, m_numSpecies, 0.0);
        for (size_t e = 0; e < m_numElemConstraints; e++) {
            for (size_t k = 0; k < m_numSpecies; k++) {
                m_formulaMatrix[e][k] = b.m_formulaMatrix[e][k];
            }
        }
        m_speciesUnknownType = b.m_speciesUnknownType;
        m_elemGlobalIndex    = b.m_elemGlobalIndex;
        PhaseName           = b.PhaseName;
        m_totalMolesInert   = b.m_totalMolesInert;
        m_isIdealSoln       = b.m_isIdealSoln;
        m_existence         = b.m_existence;
        m_MFStartIndex      = b.m_MFStartIndex;
        /*
         * Do a shallow copy because we haven' figured this out.
         */
        IndSpecies = b.IndSpecies;
        //IndSpeciesContig = b.IndSpeciesContig;

        for (size_t k = 0; k < old_num; k++) {
            if (ListSpeciesPtr[k]) {
                delete  ListSpeciesPtr[k];
                ListSpeciesPtr[k] = 0;
            }
        }
        ListSpeciesPtr.resize(m_numSpecies, 0);
        for (size_t k = 0; k < m_numSpecies; k++) {
            ListSpeciesPtr[k] =
                new vcs_SpeciesProperties(*(b.ListSpeciesPtr[k]));
        }
        m_useCanteraCalls   = b.m_useCanteraCalls;
        /*
         * Do a shallow copy of the ThermoPhase object pointer.
         * We don't duplicate the object.
         *  Um, there is no reason we couldn't do a
         *  duplicateMyselfAsThermoPhase() call here. This will
         *  have to be looked into.
         */
        TP_ptr              = b.TP_ptr;
        v_totalMoles              = b.v_totalMoles;
        Xmol_ = b.Xmol_;
        creationMoleNumbers_ = b.creationMoleNumbers_;
        creationGlobalRxnNumbers_ = b.creationGlobalRxnNumbers_;
        m_phiVarIndex       = b.m_phiVarIndex;
        m_totalVol          = b.m_totalVol;
        SS0ChemicalPotential = b.SS0ChemicalPotential;
        StarChemicalPotential = b.StarChemicalPotential;
        StarMolarVol = b.StarMolarVol;
        PartialMolarVol = b.PartialMolarVol;
        ActCoeff = b.ActCoeff;
        np_dLnActCoeffdMolNumber = b.np_dLnActCoeffdMolNumber;
        m_vcsStateStatus      = b.m_vcsStateStatus;
        m_phi               = b.m_phi;
        m_UpToDate            = false;
        m_UpToDate_AC         = false;
        m_UpToDate_VolStar    = false;
        m_UpToDate_VolPM      = false;
        m_UpToDate_GStar      = false;
        m_UpToDate_G0         = false;
        Temp_                = b.Temp_;
        Pres_                = b.Pres_;

        setState_TP(Temp_, Pres_);
        _updateMoleFractionDependencies();
    }
    return *this;
}

void vcs_VolPhase::resize(const size_t phaseNum, const size_t nspecies,
                          const size_t numElem, const char* const phaseName,
                          const double molesInert)
{
#ifdef DEBUG_MODE
    if (nspecies <= 0) {
        plogf("nspecies Error\n");
        exit(EXIT_FAILURE);
    }
    if (phaseNum < 0) {
        plogf("phaseNum should be greater than 0\n");
        exit(EXIT_FAILURE);
    }
#endif
    setTotalMolesInert(molesInert);
    m_phi = 0.0;
    m_phiVarIndex = npos;

    if (phaseNum == VP_ID_) {
        if (strcmp(PhaseName.c_str(), phaseName)) {
            plogf("Strings are different: %s %s :unknown situation\n",
                  PhaseName.c_str(), phaseName);
            exit(EXIT_FAILURE);
        }
    } else {
        VP_ID_ = phaseNum;
        if (!phaseName) {
            std::stringstream sstmp;
            sstmp << "Phase_" << VP_ID_;
            PhaseName = sstmp.str();
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

    if (ListSpeciesPtr.size() >= m_numSpecies) {
        for (size_t i = 0; i < m_numSpecies; i++) {
            if (ListSpeciesPtr[i]) {
                delete ListSpeciesPtr[i];
                ListSpeciesPtr[i] = 0;
            }
        }
    }
    ListSpeciesPtr.resize(nspecies, 0);
    for (size_t i = 0; i < nspecies; i++) {
        ListSpeciesPtr[i] = new vcs_SpeciesProperties(phaseNum, i, this);
    }

    Xmol_.resize(nspecies, 0.0);
    creationMoleNumbers_.resize(nspecies, 0.0);
    creationGlobalRxnNumbers_.resize(nspecies, npos);
    for (size_t i = 0; i < nspecies; i++) {
        Xmol_[i] = 1.0/nspecies;
        creationMoleNumbers_[i] = 1.0/nspecies;
        if (IndSpecies[i] - m_numElemConstraints >= 0) {
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
    m_UpToDate            = false;
    m_vcsStateStatus      = VCS_STATECALC_OLD;
    m_UpToDate_AC         = false;
    m_UpToDate_VolStar    = false;
    m_UpToDate_VolPM      = false;
    m_UpToDate_GStar      = false;
    m_UpToDate_G0         = false;


    elemResize(numElem);

}

void vcs_VolPhase::elemResize(const size_t numElemConstraints)
{

    m_elementNames.resize(numElemConstraints);

    m_elementActive.resize(numElemConstraints+1, 1);
    m_elementType.resize(numElemConstraints, VCS_ELEM_TYPE_ABSPOS);
    m_formulaMatrix.resize(numElemConstraints, m_numSpecies, 0.0);

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
    if (m_useCanteraCalls) {
        TP_ptr->getActivityCoefficients(VCS_DATA_PTR(ActCoeff));
    }
    m_UpToDate_AC = true;
}

double vcs_VolPhase::AC_calc_one(size_t kspec) const
{
    if (! m_UpToDate_AC) {
        _updateActCoeff();
    }
    return ActCoeff[kspec];
}

void vcs_VolPhase::_updateG0() const
{
    if (m_useCanteraCalls) {
        TP_ptr->getGibbs_ref(VCS_DATA_PTR(SS0ChemicalPotential));
    } else {
        double R = vcsUtil_gasConstant(p_VCS_UnitsFormat);
        for (size_t k = 0; k < m_numSpecies; k++) {
            size_t kglob = IndSpecies[k];
            vcs_SpeciesProperties* sProp = ListSpeciesPtr[k];
            VCS_SPECIES_THERMO* sTherm = sProp->SpeciesThermo;
            SS0ChemicalPotential[k] =
                R * (sTherm->G0_R_calc(kglob, Temp_));
        }
    }
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
    if (m_useCanteraCalls) {
        TP_ptr->getStandardChemPotentials(VCS_DATA_PTR(StarChemicalPotential));
    } else {
        double R = vcsUtil_gasConstant(p_VCS_UnitsFormat);
        for (size_t k = 0; k < m_numSpecies; k++) {
            size_t kglob = IndSpecies[k];
            vcs_SpeciesProperties* sProp = ListSpeciesPtr[k];
            VCS_SPECIES_THERMO* sTherm = sProp->SpeciesThermo;
            StarChemicalPotential[k] =
                R * (sTherm->GStar_R_calc(kglob, Temp_, Pres_));
        }
    }
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
    if (m_useCanteraCalls) {
        if (TP_ptr) {
            TP_ptr->setState_PX(Pres_, &(Xmol_[m_MFStartIndex]));
        }
    }
    if (!m_isIdealSoln) {
        m_UpToDate_AC = false;
        m_UpToDate_VolPM = false;
    }
}

const std::vector<double> & vcs_VolPhase::moleFractions() const
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
            printf("vcs_VolPhase::setMolesFractionsState: inappropriate usage\n");
            exit(EXIT_FAILURE);
        }
        m_UpToDate = false;
        m_vcsStateStatus = VCS_STATECALC_TMP;
        if (m_existence == VCS_PHASE_EXIST_ZEROEDPHASE) {
            printf("vcs_VolPhase::setMolesFractionsState: inappropriate usage\n");
            exit(EXIT_FAILURE);
        }
        m_existence = VCS_PHASE_EXIST_YES;
    } else {
        m_UpToDate = true;
        m_vcsStateStatus = vcsStateStatus;
        if (m_existence > VCS_PHASE_EXIST_NO) {
            m_existence = VCS_PHASE_EXIST_NO;
        }
    }
    double fractotal = 1.0;
    v_totalMoles = totalMoles;
    if (m_totalMolesInert > 0.0) {
        if (m_totalMolesInert > v_totalMoles) {
            printf("vcs_VolPhase::setMolesFractionsState: inerts greater than total: %g %g\n",
                   v_totalMoles,  m_totalMolesInert);
            exit(EXIT_FAILURE);
        }
        fractotal = 1.0 - m_totalMolesInert/v_totalMoles;
    }
    double sum = 0.0;
    for (size_t k = 0; k < m_numSpecies; k++) {
        Xmol_[k] = moleFractions[k];
        sum += moleFractions[k];
    }
    if (sum == 0.0) {
        printf("vcs_VolPhase::setMolesFractionsState: inappropriate usage\n");
        exit(EXIT_FAILURE);
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
    size_t kglob;
    double tmp;
    v_totalMoles = m_totalMolesInert;

    if (molesSpeciesVCS == 0) {
#ifdef DEBUG_MODE
        if (m_owningSolverObject == 0) {
            printf("vcs_VolPhase::setMolesFromVCS  shouldn't be here\n");
            exit(EXIT_FAILURE);
        }
#endif
        if (stateCalc == VCS_STATECALC_OLD) {
            molesSpeciesVCS = VCS_DATA_PTR(m_owningSolverObject->m_molNumSpecies_old);
        } else if (stateCalc == VCS_STATECALC_NEW) {
            molesSpeciesVCS = VCS_DATA_PTR(m_owningSolverObject->m_molNumSpecies_new);
        }
#ifdef DEBUG_MODE
        else {
            printf("vcs_VolPhase::setMolesFromVCS shouldn't be here\n");
            exit(EXIT_FAILURE);
        }
#endif
    }
#ifdef DEBUG_MODE
    else {
        if (m_owningSolverObject) {
            if (stateCalc == VCS_STATECALC_OLD) {
                if (molesSpeciesVCS != VCS_DATA_PTR(m_owningSolverObject->m_molNumSpecies_old)) {
                    printf("vcs_VolPhase::setMolesFromVCS shouldn't be here\n");
                    exit(EXIT_FAILURE);
                }
            } else if (stateCalc == VCS_STATECALC_NEW) {
                if (molesSpeciesVCS != VCS_DATA_PTR(m_owningSolverObject->m_molNumSpecies_new)) {
                    printf("vcs_VolPhase::setMolesFromVCS shouldn't be here\n");
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
#endif

    for (size_t k = 0; k < m_numSpecies; k++) {
        if (m_speciesUnknownType[k] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            kglob = IndSpecies[k];
            v_totalMoles += std::max(0.0, molesSpeciesVCS[kglob]);
        }
    }
    if (v_totalMoles > 0.0) {
        for (size_t k = 0; k < m_numSpecies; k++) {
            if (m_speciesUnknownType[k] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                kglob = IndSpecies[k];
                tmp = std::max(0.0, molesSpeciesVCS[kglob]);
                Xmol_[k] = tmp / v_totalMoles;
            }
        }
        m_existence = VCS_PHASE_EXIST_YES;
    } else {
        // This is where we will start to store a better approximation
        // for the mole fractions, when the phase doesn't exist.
        // This is currently unimplemented.
        //for (int k = 0; k < m_numSpecies; k++) {
        //    Xmol_[k] = 1.0 / m_numSpecies;
        //}
        m_existence = VCS_PHASE_EXIST_NO;
    }
    /*
     * Update the electric potential if it is a solution variable
     * in the equation system
     */
    if (m_phiVarIndex != npos) {
        kglob = IndSpecies[m_phiVarIndex];
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
    if (m_totalMolesInert > 0.0) {
        m_existence = VCS_PHASE_EXIST_ALWAYS;
    }

    /*
     * If stateCalc is old and the total moles is positive,
     * then we have a valid state. If the phase went away, it would
     * be a valid starting point for F_k's. So, save the state.
     */
    if (stateCalc == VCS_STATECALC_OLD) {
        if (v_totalMoles > 0.0) {
            vcs_dcopy(VCS_DATA_PTR(creationMoleNumbers_), VCS_DATA_PTR(Xmol_), m_numSpecies);

        }
    }

    /*
     * Set flags indicating we are up to date with the VCS state vector.
     */
    m_UpToDate = true;
    m_vcsStateStatus = stateCalc;
}

void vcs_VolPhase::setMolesFromVCSCheck(const int vcsStateStatus,
                                        const double* molesSpeciesVCS,
                                        const double* const TPhMoles)
{
    setMolesFromVCS(vcsStateStatus, molesSpeciesVCS);
    /*
     * Check for consistency with TPhMoles[]
     */
    double Tcheck = TPhMoles[VP_ID_];
    if (Tcheck != v_totalMoles) {
        if (vcs_doubleEqual(Tcheck, v_totalMoles)) {
            Tcheck = v_totalMoles;
        } else {
            plogf("vcs_VolPhase::setMolesFromVCSCheck: "
                  "We have a consistency problem: %21.16g %21.16g\n",
                  Tcheck, v_totalMoles);
            exit(EXIT_FAILURE);
        }
    }
}

void vcs_VolPhase::updateFromVCS_MoleNumbers(const int vcsStateStatus)
{
    if (!m_UpToDate || (vcsStateStatus != m_vcsStateStatus)) {
        if (vcsStateStatus == VCS_STATECALC_OLD || vcsStateStatus == VCS_STATECALC_NEW) {
            if (m_owningSolverObject) {
                setMolesFromVCS(vcsStateStatus);
            }
        }
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
        (void) _updateVolPM();
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
    if (m_useCanteraCalls) {
        TP_ptr->setElectricPotential(m_phi);
    }
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
    if (Temp_ == temp) {
        if (Pres_ == pres) {
            return;
        }
    }
    if (m_useCanteraCalls) {
        TP_ptr->setElectricPotential(m_phi);
        TP_ptr->setState_TP(temp, pres);
    }
    Temp_ = temp;
    Pres_ = pres;
    m_UpToDate_AC      = false;
    m_UpToDate_VolStar = false;
    m_UpToDate_VolPM   = false;
    m_UpToDate_GStar   = false;
    m_UpToDate_G0      = false;
}

void vcs_VolPhase::setState_T(const double temp)
{
    setState_TP(temp, Pres_);
}

void vcs_VolPhase::_updateVolStar() const
{
    if (m_useCanteraCalls) {
        TP_ptr->getStandardVolumes(VCS_DATA_PTR(StarMolarVol));
    } else {
        for (size_t k = 0; k < m_numSpecies; k++) {
            size_t kglob = IndSpecies[k];
            vcs_SpeciesProperties* sProp = ListSpeciesPtr[k];
            VCS_SPECIES_THERMO* sTherm = sProp->SpeciesThermo;
            StarMolarVol[k] = (sTherm->VolStar_calc(kglob, Temp_, Pres_));
        }
    }
    m_UpToDate_VolStar = true;
}

double vcs_VolPhase::VolStar_calc_one(size_t kspec) const
{
    if (!m_UpToDate_VolStar) {
        _updateVolStar();
    }
    return StarMolarVol[kspec];
}

double vcs_VolPhase::_updateVolPM() const
{
    if (m_useCanteraCalls) {
        TP_ptr->getPartialMolarVolumes(VCS_DATA_PTR(PartialMolarVol));
    } else {
        for (size_t k = 0; k < m_numSpecies; k++) {
            size_t kglob = IndSpecies[k];
            vcs_SpeciesProperties* sProp = ListSpeciesPtr[k];
            VCS_SPECIES_THERMO* sTherm = sProp->SpeciesThermo;
            StarMolarVol[k] = (sTherm->VolStar_calc(kglob, Temp_, Pres_));
        }
        for (size_t k = 0; k < m_numSpecies; k++) {
            PartialMolarVol[k] = StarMolarVol[k];
        }
    }

    m_totalVol = 0.0;
    for (size_t k = 0; k < m_numSpecies; k++) {
        m_totalVol += PartialMolarVol[k] * Xmol_[k];
    }
    m_totalVol *= v_totalMoles;

    if (m_totalMolesInert > 0.0) {
        if (m_gasPhase) {
            double volI = m_totalMolesInert * Cantera::GasConstant * Temp_ / Pres_;
            m_totalVol += volI;
        } else {
            printf("unknown situation\n");
            exit(EXIT_FAILURE);
        }
    }
    m_UpToDate_VolPM = true;
    return m_totalVol;
}

void vcs_VolPhase::_updateLnActCoeffJac()
{
    double phaseTotalMoles = v_totalMoles;
    if (phaseTotalMoles < 1.0E-14) {
        phaseTotalMoles = 1.0;
    }

    /*
     * Evaluate the current base activity coefficients if necessary
     */
    if (!m_UpToDate_AC) {
        _updateActCoeff();
    }
    if (!TP_ptr) {
        return;
    }
    TP_ptr->getdlnActCoeffdlnN(m_numSpecies, &np_dLnActCoeffdMolNumber[0][0]);
    for (size_t j = 0; j < m_numSpecies; j++) {
        double moles_j_base = phaseTotalMoles * Xmol_[j];
        double* const np_lnActCoeffCol = np_dLnActCoeffdMolNumber[j];
        if (moles_j_base < 1.0E-200) {
            moles_j_base = 1.0E-7 * moles_j_base + 1.0E-13 * phaseTotalMoles + 1.0E-150;
        }
        for (size_t k = 0; k < m_numSpecies; k++) {
            np_lnActCoeffCol[k] = np_lnActCoeffCol[k] * phaseTotalMoles / moles_j_base;
        }
    }

    double deltaMoles_j = 0.0;
    // Make copies of ActCoeff and Xmol_ for use in taking differences
    std::vector<double> ActCoeff_Base(ActCoeff);
    std::vector<double> Xmol_Base(Xmol_);
    double TMoles_base = phaseTotalMoles;

    /*
     *  Loop over the columns species to be deltad
     */
    for (size_t j = 0; j < m_numSpecies; j++) {
        /*
         * Calculate a value for the delta moles of species j
         * -> Note Xmol_[] and Tmoles are always positive or zero
         *    quantities.
         */
        double moles_j_base = phaseTotalMoles * Xmol_Base[j];
        deltaMoles_j = 1.0E-7 * moles_j_base + 1.0E-13 * phaseTotalMoles + 1.0E-150;
        /*
         * Now, update the total moles in the phase and all of the
         * mole fractions based on this.
         */
        phaseTotalMoles = TMoles_base + deltaMoles_j;
        for (size_t k = 0; k < m_numSpecies; k++) {
            Xmol_[k] = Xmol_Base[k] * TMoles_base / phaseTotalMoles;
        }
        Xmol_[j] = (moles_j_base + deltaMoles_j) / phaseTotalMoles;

        /*
         * Go get new values for the activity coefficients.
         * -> Note this calls setState_PX();
         */
        _updateMoleFractionDependencies();
        _updateActCoeff();
        /*
         * Calculate the column of the matrix
         */
        double* const np_lnActCoeffCol = np_dLnActCoeffdMolNumber[j];
        for (size_t k = 0; k < m_numSpecies; k++) {
            double tmp;
            tmp = (ActCoeff[k] - ActCoeff_Base[k]) /
                  ((ActCoeff[k] + ActCoeff_Base[k]) * 0.5 * deltaMoles_j);
            if (fabs(tmp - np_lnActCoeffCol[k]) > 1.0E-4 * fabs(tmp) +  fabs(np_lnActCoeffCol[k])) {
                //  printf(" we have an error\n");

            }
            //tmp = lnActCoeffCol[k];

        }
        /*
         * Revert to the base case Xmol_, v_totalMoles
         */
        v_totalMoles = TMoles_base;
        vcs_vdcopy(Xmol_, Xmol_Base, m_numSpecies);
    }
    /*
     * Go get base values for the activity coefficients.
     * -> Note this calls setState_TPX() again;
     * -> Just wanted to make sure that cantera is in sync
     *    with VolPhase after this call.
     */
    setMoleFractions(VCS_DATA_PTR(Xmol_Base));
    _updateMoleFractionDependencies();
    _updateActCoeff();
}

void vcs_VolPhase::sendToVCS_LnActCoeffJac(double* const* const np_LnACJac_VCS)
{
    /*
     * update the Ln Act Coeff jacobian entries with respect to the
     * mole number of species in the phase -> we always assume that
     * they are out of date.
     */
    _updateLnActCoeffJac();

    /*
     *  Now copy over the values
     */
    for (size_t j = 0; j < m_numSpecies; j++) {
        size_t jglob = IndSpecies[j];
        double* const np_lnACJacVCS_col = np_LnACJac_VCS[jglob];
        const double* const np_lnACJac_col = np_dLnActCoeffdMolNumber[j];
        for (size_t k = 0; k < m_numSpecies; k++) {
            size_t kglob = IndSpecies[k];
            np_lnACJacVCS_col[kglob] = np_lnACJac_col[k];
        }
    }
}

void vcs_VolPhase::setPtrThermoPhase(Cantera::ThermoPhase* tp_ptr)
{
    TP_ptr = tp_ptr;
    if (TP_ptr) {
        m_useCanteraCalls = true;
        Temp_ = TP_ptr->temperature();
        Pres_ = TP_ptr->pressure();
        setState_TP(Temp_, Pres_);
        p_VCS_UnitsFormat = VCS_UNITS_MKS;
        m_phi = TP_ptr->electricPotential();
        size_t nsp = TP_ptr->nSpecies();
        size_t nelem = TP_ptr->nElements();
        if (nsp !=  m_numSpecies) {
            if (m_numSpecies != 0) {
                plogf("Warning Nsp != NVolSpeces: %d %d \n", nsp, m_numSpecies);
            }
            resize(VP_ID_, nsp, nelem, PhaseName.c_str());
        }
        TP_ptr->getMoleFractions(VCS_DATA_PTR(Xmol_));
        vcs_dcopy(VCS_DATA_PTR(creationMoleNumbers_), VCS_DATA_PTR(Xmol_), m_numSpecies);
        _updateMoleFractionDependencies();

        /*
         *  figure out ideal solution tag
         */
        if (nsp == 1) {
            m_isIdealSoln = true;
        } else {
            int eos = TP_ptr->eosType();
            switch (eos) {
            case Cantera::cIdealGas:
            case Cantera::cIncompressible:
            case Cantera::cSurf:
            case Cantera::cMetal:
            case Cantera::cStoichSubstance:
            case Cantera::cSemiconductor:
            case Cantera::cLatticeSolid:
            case Cantera::cLattice:
            case Cantera::cEdge:
            case Cantera::cIdealSolidSolnPhase:
                m_isIdealSoln = true;
                break;
            default:
                m_isIdealSoln = false;
            };
        }
    } else {
        m_useCanteraCalls = false;
    }
}

const Cantera::ThermoPhase* vcs_VolPhase::ptrThermoPhase() const
{
    return TP_ptr;
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
        const std::vector<size_t> &creationGlobalRxnNumbers)
{
    vcs_dcopy(VCS_DATA_PTR(creationMoleNumbers_), n_k, m_numSpecies);
    for (size_t k = 0; k < m_numSpecies; k++) { 
        creationGlobalRxnNumbers_[k] = creationGlobalRxnNumbers[k];
    }
}

const std::vector<double> & vcs_VolPhase::creationMoleNumbers(std::vector<size_t> &creationGlobalRxnNumbers) const
{
    creationGlobalRxnNumbers = creationGlobalRxnNumbers_;
    return creationMoleNumbers_;
}

void vcs_VolPhase::setTotalMoles(const double totalMols)
{
    v_totalMoles = totalMols;
    if (m_totalMolesInert > 0.0) {
        m_existence = VCS_PHASE_EXIST_ALWAYS;
#ifdef DEBUG_MODE
        if (totalMols < m_totalMolesInert) {
            printf(" vcs_VolPhase::setTotalMoles:: ERROR totalMoles "
                   "less than inert moles: %g %g\n",
                   totalMols, m_totalMolesInert);
            exit(EXIT_FAILURE);
        }
#endif
    } else {
        if (m_singleSpecies && (m_phiVarIndex == 0)) {
            m_existence =  VCS_PHASE_EXIST_ALWAYS;
        } else {
            if (totalMols > 0.0) {
                m_existence = VCS_PHASE_EXIST_YES;
            } else {
                m_existence = VCS_PHASE_EXIST_NO;
            }
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

std::string string16_EOSType(int EOSType)
{
    char st[32];
    st[16] = '\0';
    switch (EOSType) {
    case VCS_EOS_CONSTANT:
        sprintf(st,"Constant        ");
        break;
    case VCS_EOS_IDEAL_GAS:
        sprintf(st,"Ideal Gas       ");
        break;
    case  VCS_EOS_STOICH_SUB:
        sprintf(st,"Stoich Sub      ");
        break;
    case VCS_EOS_IDEAL_SOLN:
        sprintf(st,"Ideal Soln      ");
        break;
    case VCS_EOS_DEBEYE_HUCKEL:
        sprintf(st,"Debeye Huckel   ");
        break;
    case VCS_EOS_REDLICK_KWONG:
        sprintf(st,"Redlick_Kwong   ");
        break;
    case VCS_EOS_REGULAR_SOLN:
        sprintf(st,"Regular Soln    ");
        break;
    default:
        sprintf(st,"UnkType: %-7d", EOSType);
        break;
    }
    st[16] = '\0';
    return st;
}

bool vcs_VolPhase::isIdealSoln() const
{
    return m_isIdealSoln;
}

bool vcs_VolPhase::usingCanteraCalls() const
{
    return m_useCanteraCalls;
}

size_t vcs_VolPhase::phiVarIndex() const
{
    return m_phiVarIndex;
}

void vcs_VolPhase::setPhiVarIndex(size_t phiVarIndex)
{
    m_phiVarIndex = phiVarIndex;
    m_speciesUnknownType[m_phiVarIndex] = VCS_SPECIES_TYPE_INTERFACIALVOLTAGE;
    if (m_singleSpecies) {
        if (m_phiVarIndex == 0) {
            m_existence = VCS_PHASE_EXIST_ALWAYS;
        }
    }
}

vcs_SpeciesProperties* vcs_VolPhase::speciesProperty(const size_t kindex)
{
    return  ListSpeciesPtr[kindex];
}

int vcs_VolPhase::exists() const
{
    return m_existence;
}

void vcs_VolPhase::setExistence(const int existence)
{
    if (existence == VCS_PHASE_EXIST_NO || existence == VCS_PHASE_EXIST_ZEROEDPHASE) {
        if (v_totalMoles != 0.0) {
#ifdef DEBUG_MODE
            plogf("vcs_VolPhase::setExistence setting false existence for phase with moles");
            plogendl();
            exit(EXIT_FAILURE);
#else
            v_totalMoles = 0.0;
#endif
        }
    }
#ifdef DEBUG_MODE
    else {
        if (m_totalMolesInert == 0.0) {
            if (v_totalMoles == 0.0) {
                if (!m_singleSpecies  || m_phiVarIndex != 0) {
                    plogf("vcs_VolPhase::setExistence setting true existence for phase with no moles");
                    plogendl();
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
#endif
#ifdef DEBUG_MODE
    if (m_singleSpecies) {
        if (m_phiVarIndex == 0) {
            if (existence == VCS_PHASE_EXIST_NO || existence == VCS_PHASE_EXIST_ZEROEDPHASE) {
                plogf("vcs_VolPhase::Trying to set existence of an electron phase to false");
                plogendl();
                exit(EXIT_FAILURE);
            }
        }
    }
#endif
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

void vcs_VolPhase::setTotalMolesInert(const double tMolesInert)
{
    if (m_totalMolesInert != tMolesInert) {
        m_UpToDate = false;
        m_UpToDate_AC = false;
        m_UpToDate_VolStar = false;
        m_UpToDate_VolPM = false;
        m_UpToDate_GStar = false;
        m_UpToDate_G0 = false;
        v_totalMoles += (tMolesInert - m_totalMolesInert);
        m_totalMolesInert = tMolesInert;
    }
    if (m_totalMolesInert > 0.0) {
        m_existence = VCS_PHASE_EXIST_ALWAYS;
    } else if (m_singleSpecies && (m_phiVarIndex == 0)) {
        m_existence = VCS_PHASE_EXIST_ALWAYS;
    } else {
        if (v_totalMoles > 0.0) {
            m_existence = VCS_PHASE_EXIST_YES;
        } else {
            m_existence = VCS_PHASE_EXIST_NO;
        }
    }
}

double vcs_VolPhase::totalMolesInert() const
{
    return m_totalMolesInert;
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

std::string vcs_VolPhase::elementName(const size_t e) const
{
    return m_elementNames[e];
}

//! This function decides whether a phase has charged species or not.
static bool hasChargedSpecies(const Cantera::ThermoPhase* const tPhase)
{
    for (size_t k = 0; k < tPhase->nSpecies(); k++) {
        if (tPhase->charge(k) != 0.0) {
            return true;
        }
    }
    return false;
}

/*!
 *  This utility routine decides whether a Cantera ThermoPhase needs
 *  a constraint equation representing the charge neutrality of the
 *  phase. It does this by searching for charged species. If it
 *  finds one, and if the phase needs one, then it returns true.
 */
static bool chargeNeutralityElement(const Cantera::ThermoPhase* const tPhase)
{
    int hasCharge = hasChargedSpecies(tPhase);
    if (tPhase->chargeNeutralityNecessary()) {
        if (hasCharge) {
            return true;
        }
    }
    return false;
}

size_t vcs_VolPhase::transferElementsFM(const Cantera::ThermoPhase* const tPhase)
{
    size_t e, k, eT;
    std::string ename;
    size_t eFound = npos;
    size_t nebase = tPhase->nElements();
    size_t ne = nebase;
    size_t ns = tPhase->nSpecies();

    /*
     * Decide whether we need an extra element constraint for charge
     * neutrality of the phase
     */
    bool cne = chargeNeutralityElement(tPhase);
    if (cne) {
        ChargeNeutralityElement = ne;
        ne++;
    }

    /*
     * Assign and malloc structures
     */
    elemResize(ne);


    if (ChargeNeutralityElement != npos) {
        m_elementType[ChargeNeutralityElement] = VCS_ELEM_TYPE_CHARGENEUTRALITY;
    }

    if (hasChargedSpecies(tPhase)) {
        if (cne) {
            /*
             * We need a charge neutrality constraint.
             * We also have an Electron Element. These are
             * duplicates of each other. To avoid trouble with
             * possible range error conflicts, sometimes we eliminate
             * the Electron condition. Flag that condition for elimination
             * by toggling the ElActive variable. If we find we need it
             * later, we will retoggle ElActive to true.
             */
            for (eT = 0; eT < nebase; eT++) {
                ename = tPhase->elementName(eT);
                if (ename == "E") {
                    eFound = eT;
                    m_elementActive[eT] = 0;
                    m_elementType[eT] = VCS_ELEM_TYPE_ELECTRONCHARGE;
                }
            }
        } else {
            for (eT = 0; eT < nebase; eT++) {
                ename = tPhase->elementName(eT);
                if (ename == "E") {
                    eFound = eT;
                    m_elementType[eT] = VCS_ELEM_TYPE_ELECTRONCHARGE;
                }
            }
        }
        if (eFound == npos) {
            eFound = ne;
            m_elementType[ne] = VCS_ELEM_TYPE_ELECTRONCHARGE;
            m_elementActive[ne] = 0;
            std::string ename = "E";
            m_elementNames[ne] = ename;
            ne++;
            elemResize(ne);
        }

    }

    m_formulaMatrix.resize(ne, ns, 0.0);

    m_speciesUnknownType.resize(ns, VCS_SPECIES_TYPE_MOLNUM);

    elemResize(ne);

    e = 0;
    for (eT = 0; eT < nebase; eT++) {
        ename = tPhase->elementName(eT);
        m_elementNames[e] = ename;
        m_elementType[e] = tPhase->elementType(eT);
        e++;
    }

    if (cne) {
        std::string pname = tPhase->id();
        if (pname == "") {
            std::stringstream sss;
            sss << "phase" << VP_ID_;
            pname = sss.str();
        }
        ename = "cn_" + pname;
        e = ChargeNeutralityElement;
        m_elementNames[e] = ename;
    }

    double* const* const fm = m_formulaMatrix.baseDataAddr();
    for (k = 0; k < ns; k++) {
        e = 0;
        for (eT = 0; eT < nebase; eT++) {
            fm[e][k] = tPhase->nAtoms(k, eT);
            e++;
        }
        if (eFound != npos) {
            fm[eFound][k] = - tPhase->charge(k);
        }
    }

    if (cne) {
        for (k = 0; k < ns; k++) {
            fm[ChargeNeutralityElement][k] = tPhase->charge(k);
        }
    }


    /*
     * Here, we figure out what is the species types are
     * The logic isn't set in stone, and is just for a particular type
     * of problem that I'm solving first.
     */
    if (ns == 1) {
        if (tPhase->charge(0) != 0.0) {
            m_speciesUnknownType[0] = VCS_SPECIES_TYPE_INTERFACIALVOLTAGE;
            setPhiVarIndex(0);
        }
    }

    return ne;
}

int vcs_VolPhase::elementType(const size_t e) const
{
    return m_elementType[e];
}

void vcs_VolPhase::setElementType(const size_t e, const int eType)
{
    m_elementType[e] = eType;
}

double const* const* vcs_VolPhase::getFormulaMatrix() const
{
    return m_formulaMatrix.constBaseDataAddr();
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

}
