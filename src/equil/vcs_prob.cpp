/**
 * @file vcs_prob.cpp
 *  Implementation for the Interface class for the vcs thermo
 *  equilibrium solver package,
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/equil/vcs_prob.h"
#include "cantera/equil/vcs_VolPhase.h"
#include "cantera/equil/vcs_species_thermo.h"
#include "cantera/equil/vcs_internal.h"

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/MolalityVPSSTP.h"

#include <cstdio>

using namespace std;

namespace VCSnonideal
{

VCS_PROB::VCS_PROB(size_t nsp, size_t nel, size_t nph) :
    prob_type(VCS_PROBTYPE_TP),
    nspecies(nsp),
    NSPECIES0(0),
    ne(nel),
    NE0(0),
    NPhase(nph),
    NPHASE0(0),
    T(298.15),
    PresPA(1.0),
    Vol(0.0),
    m_VCS_UnitsFormat(VCS_UNITS_UNITLESS),
/* Set the units for the chemical potential data to be
 * unitless */
    iest(-1),    /* The default is to not expect an initial estimate
                  * of the species concentrations */
    tolmaj(1.0E-8),
    tolmin(1.0E-6),
    m_Iterations(0),
    m_NumBasisOptimizations(0),
    m_printLvl(0),
    vcs_debug_print_lvl(0)
{
    NSPECIES0 = nspecies;
    if (nspecies <= 0) {
        plogf("number of species is zero or neg\n");
        exit(EXIT_FAILURE);
    }
    NE0       = ne;
    if (ne <= 0) {
        plogf("number of elements is zero or neg\n");
        exit(EXIT_FAILURE);
    }
    NPHASE0   = NPhase;
    if (NPhase <= 0) {
        plogf("number of phases is zero or neg\n");
        exit(EXIT_FAILURE);
    }
    if (nspecies < NPhase) {
        plogf("number of species is less than number of phases\n");
        exit(EXIT_FAILURE);
    }

    m_gibbsSpecies.resize(nspecies, 0.0);
    w.resize(nspecies, 0.0);
    mf.resize(nspecies, 0.0);
    gai.resize(ne, 0.0);
    FormulaMatrix.resize(ne, nspecies, 0.0);
    SpeciesUnknownType.resize(nspecies, VCS_SPECIES_TYPE_MOLNUM);
    VolPM.resize(nspecies, 0.0);
    PhaseID.resize(nspecies, npos);
    SpName.resize(nspecies, "");
    ElName.resize(ne, "");
    m_elType.resize(ne, VCS_ELEM_TYPE_ABSPOS);
    ElActive.resize(ne, 1);
    WtSpecies.resize(nspecies, 0.0);
    Charge.resize(nspecies, 0.0);
    SpeciesThermo.resize(nspecies,0);
    for (size_t kspec = 0; kspec < nspecies; kspec++) {
        VCS_SPECIES_THERMO* ts_tmp = new VCS_SPECIES_THERMO(0, 0);
        if (ts_tmp == 0) {
            plogf("Failed to init a ts struct\n");
            exit(EXIT_FAILURE);
        }
        SpeciesThermo[kspec] = ts_tmp;
    }
    VPhaseList.resize(nph, 0);
    for (size_t iphase = 0; iphase < NPhase; iphase++) {
        VPhaseList[iphase] = new vcs_VolPhase();
    }
}

VCS_PROB::~VCS_PROB()
{
    for (size_t i = 0; i < nspecies; i++) {
        delete SpeciesThermo[i];
        SpeciesThermo[i] = 0;
    }
    for (size_t iph = 0; iph < NPhase; iph++) {
        delete VPhaseList[iph];
        VPhaseList[iph] = 0;
    }
}

void VCS_PROB::resizePhase(size_t nPhase, int force)
{
    if (force || nPhase > NPHASE0) {
        NPHASE0 = nPhase;
    }
}

void VCS_PROB::resizeSpecies(size_t nsp, int force)
{
    if (force || nsp > NSPECIES0) {
        m_gibbsSpecies.resize(nsp, 0.0);
        w.resize(nsp, 0.0);
        mf.resize(nsp, 0.0);
        FormulaMatrix.resize(NE0, nsp, 0.0);
        SpeciesUnknownType.resize(nsp, VCS_SPECIES_TYPE_MOLNUM);
        VolPM.resize(nsp, 0.0);
        PhaseID.resize(nsp, 0);
        SpName.resize(nsp, "");
        WtSpecies.resize(nsp, 0.0);
        Charge.resize(nsp, 0.0);
        NSPECIES0 = nsp;
        if (nspecies > NSPECIES0) {
            nspecies = NSPECIES0;
            plogf("shouldn't be here\n");
            exit(EXIT_FAILURE);
        }
    }
}

void VCS_PROB::resizeElements(size_t nel, int force)
{
    if (force || nel > NE0) {
        gai.resize(nel, 0.0);
        FormulaMatrix.resize(nel, NSPECIES0, 0.0);
        ElName.resize(nel, "");
        m_elType.resize(nel, VCS_ELEM_TYPE_ABSPOS);
        ElActive.resize(nel, 1);
        NE0 = nel;
        if (ne > NE0) {
            ne = NE0;
        }
    }
}

void VCS_PROB::set_gai()
{
    double* ElemAbund = VCS_DATA_PTR(gai);
    double* const* const fm = FormulaMatrix.baseDataAddr();
    vcs_dzero(ElemAbund, ne);

    for (size_t j = 0; j < ne; j++) {
        for (size_t kspec = 0; kspec < nspecies; kspec++) {
            if (SpeciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                ElemAbund[j] += fm[j][kspec] * w[kspec];
            }
        }
    }
}

static void print_space(int num)
{
    for (int j = 0; j < num; j++) {
        (void) plogf(" ");
    }
}

static void print_char(const char letter, const int num)
{
    for (int i = 0; i < num; i++) {
        plogf("%c", letter);
    }
}

void VCS_PROB::prob_report(int print_lvl)
{
    m_printLvl = print_lvl;
    vcs_VolPhase* Vphase = 0;
    /*
     *          Printout the species information: PhaseID's and mole nums
     */
    if (m_printLvl > 0) {
        plogf("\n");
        print_char('=', 80);
        plogf("\n");
        print_char('=', 20);
        plogf(" VCS_PROB: PROBLEM STATEMENT ");
        print_char('=', 31);
        plogf("\n");
        print_char('=', 80);
        plogf("\n");

        plogf("\n");
        if (prob_type == 0) {
            plogf("\tSolve a constant T, P problem:\n");
            plogf("\t\tT    = %g K\n", T);
            double pres_atm = PresPA / 1.01325E5;

            plogf("\t\tPres = %g atm\n", pres_atm);
        } else {
            plogf("\tUnknown problem type\n");
            exit(EXIT_FAILURE);
        }
        plogf("\n");
        plogf("             Phase IDs of species\n");
        plogf("            species     phaseID        phaseName   ");
        plogf(" Initial_Estimated_Moles   Species_Type\n");
        for (size_t i = 0; i < nspecies; i++) {
            Vphase = VPhaseList[PhaseID[i]];
            plogf("%16s      %5d   %16s", SpName[i].c_str(), PhaseID[i],
                  Vphase->PhaseName.c_str());
            if (iest >= 0) {
                plogf("             %-10.5g",  w[i]);
            } else {
                plogf("                N/A");
            }
            if (SpeciesUnknownType[i] == VCS_SPECIES_TYPE_MOLNUM) {
                plogf("                 Mol_Num");
            } else if (SpeciesUnknownType[i] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                plogf("                 Voltage");
            } else {
                plogf("                        ");
            }
            plogf("\n");
        }

        /*
         *   Printout of the Phase structure information
         */
        plogf("\n");
        print_char('-', 80);
        plogf("\n");
        plogf("             Information about phases\n");
        plogf("  PhaseName    PhaseNum SingSpec  GasPhase   "
              " EqnState    NumSpec");
        plogf("  TMolesInert      TKmoles\n");

        for (size_t iphase = 0; iphase < NPhase; iphase++) {
            Vphase = VPhaseList[iphase];
            std::string EOS_cstr = string16_EOSType(Vphase->m_eqnState);
            plogf("%16s %5d %5d %8d ", Vphase->PhaseName.c_str(),
                  Vphase->VP_ID_, Vphase->m_singleSpecies, Vphase->m_gasPhase);
            plogf("%16s %8d %16e ", EOS_cstr.c_str(),
                  Vphase->nSpecies(), Vphase->totalMolesInert());
            if (iest >= 0) {
                plogf("%16e\n",  Vphase->totalMoles());
            } else {
                plogf("   N/A\n");
            }
        }

        plogf("\nElemental Abundances:    ");
        plogf("         Target_kmol    ElemType ElActive\n");
        double fac = 1.0;
        if (m_VCS_UnitsFormat == VCS_UNITS_MKS) {
            //fac = 1.0E3;
            fac = 1.0;
        }
        for (size_t i = 0; i < ne; ++i) {
            print_space(26);
            plogf("%-2.2s", ElName[i].c_str());
            plogf("%20.12E  ", fac * gai[i]);
            plogf("%3d       %3d\n", m_elType[i], ElActive[i]);
        }

        plogf("\nChemical Potentials:  ");
        if (m_VCS_UnitsFormat == VCS_UNITS_UNITLESS) {
            plogf("(unitless)");
        } else if (m_VCS_UnitsFormat == VCS_UNITS_KCALMOL) {
            plogf("(kcal/gmol)");
        } else if (m_VCS_UnitsFormat == VCS_UNITS_KJMOL) {
            plogf("(kJ/gmol)");
        } else if (m_VCS_UnitsFormat == VCS_UNITS_KELVIN) {
            plogf("(Kelvin)");
        } else if (m_VCS_UnitsFormat == VCS_UNITS_MKS) {
            plogf("(J/kmol)");
        }
        plogf("\n");
        plogf("             Species       (phase)    "
              "    SS0ChemPot       StarChemPot\n");
        for (size_t iphase = 0; iphase < NPhase; iphase++) {
            Vphase = VPhaseList[iphase];
            Vphase->setState_TP(T, PresPA);
            for (size_t kindex = 0; kindex < Vphase->nSpecies(); kindex++) {
                size_t kglob = Vphase->spGlobalIndexVCS(kindex);
                plogf("%16s ", SpName[kglob].c_str());
                if (kindex == 0) {
                    plogf("%16s", Vphase->PhaseName.c_str());
                } else {
                    plogf("                ");
                }

                plogf("%16g   %16g\n", Vphase->G0_calc_one(kindex),
                      Vphase->GStar_calc_one(kindex));
            }
        }
        plogf("\n");
        print_char('=', 80);
        plogf("\n");
        print_char('=', 20);
        plogf(" VCS_PROB: END OF PROBLEM STATEMENT ");
        print_char('=', 24);
        plogf("\n");
        print_char('=', 80);
        plogf("\n\n");
    }
}

void VCS_PROB::addPhaseElements(vcs_VolPhase* volPhase)
{
    size_t e, eVP;
    size_t foundPos = npos;
    size_t neVP = volPhase->nElemConstraints();
    std::string en;
    std::string enVP;
    /*
     * Loop through the elements in the vol phase object
     */
    for (eVP = 0; eVP < neVP; eVP++) {
        foundPos = npos;
        enVP = volPhase->elementName(eVP);
        /*
         * Search for matches with the existing elements.
         * If found, then fill in the entry in the global
         * mapping array.
         */
        for (e = 0; e < ne; e++) {
            en = ElName[e];
            if (!strcmp(enVP.c_str(), en.c_str())) {
                volPhase->setElemGlobalIndex(eVP, e);
                foundPos = e;
            }
        }
        if (foundPos == npos) {
            int elType = volPhase->elementType(eVP);
            int elactive = volPhase->elementActive(eVP);
            e = addElement(enVP.c_str(), elType, elactive);
            volPhase->setElemGlobalIndex(eVP, e);
        }
    }
}

size_t VCS_PROB::addElement(const char* elNameNew, int elType, int elactive)
{
    if (!elNameNew) {
        plogf("error: element must have a name\n");
        exit(EXIT_FAILURE);
    }
    size_t nel = ne + 1;
    resizeElements(nel, 1);
    ne = nel;
    ElName[ne-1] = elNameNew;
    m_elType[ne-1] = elType;
    ElActive[ne-1] = elactive;
    return ne - 1;
}

size_t VCS_PROB::addOnePhaseSpecies(vcs_VolPhase* volPhase, size_t k, size_t kT)
{
    size_t e, eVP;
    if (kT > nspecies) {
        /*
         * Need to expand the number of species here
         */
        plogf("Shouldn't be here\n");
        exit(EXIT_FAILURE);
    }
    double const* const* const fm = volPhase->getFormulaMatrix();
    for (eVP = 0; eVP < volPhase->nElemConstraints(); eVP++) {
        e = volPhase->elemGlobalIndex(eVP);
#ifdef DEBUG_MODE
        if (e == npos) {
            exit(EXIT_FAILURE);
        }
#endif
        FormulaMatrix[e][kT] = fm[eVP][k];
    }
    /*
     * Tell the phase object about the current position of the
     * species within the global species vector
     */
    volPhase->setSpGlobalIndexVCS(k, kT);
    return kT;
}

void VCS_PROB::reportCSV(const std::string& reportFile)
{
    size_t k;
    size_t istart;

    double vol = 0.0;
    string sName;

    FILE* FP = fopen(reportFile.c_str(), "w");
    if (!FP) {
        plogf("Failure to open file\n");
        exit(EXIT_FAILURE);
    }
    double Temp = T;

    std::vector<double> volPM(nspecies, 0.0);
    std::vector<double> activity(nspecies, 0.0);
    std::vector<double> ac(nspecies, 0.0);
    std::vector<double> mu(nspecies, 0.0);
    std::vector<double> mu0(nspecies, 0.0);
    std::vector<double> molalities(nspecies, 0.0);

    vol = 0.0;
    size_t iK = 0;
    for (size_t iphase = 0; iphase < NPhase; iphase++) {
        istart = iK;
        vcs_VolPhase* volP = VPhaseList[iphase];
        //const Cantera::ThermoPhase *tptr = volP->ptrThermoPhase();
        size_t nSpeciesPhase = volP->nSpecies();
        volPM.resize(nSpeciesPhase, 0.0);
        volP->sendToVCS_VolPM(VCS_DATA_PTR(volPM));

        double TMolesPhase = volP->totalMoles();
        double VolPhaseVolumes = 0.0;
        for (k = 0; k < nSpeciesPhase; k++) {
            iK++;
            VolPhaseVolumes += volPM[istart + k] * mf[istart + k];
        }
        VolPhaseVolumes *= TMolesPhase;
        vol += VolPhaseVolumes;
    }

    fprintf(FP,"--------------------- VCS_MULTIPHASE_EQUIL FINAL REPORT"
            " -----------------------------\n");
    fprintf(FP,"Temperature  = %11.5g kelvin\n", Temp);
    fprintf(FP,"Pressure     = %11.5g Pascal\n", PresPA);
    fprintf(FP,"Total Volume = %11.5g m**3\n", vol);
    fprintf(FP,"Number Basis optimizations = %d\n", m_NumBasisOptimizations);
    fprintf(FP,"Number VCS iterations = %d\n", m_Iterations);

    iK = 0;
    for (size_t iphase = 0; iphase < NPhase; iphase++) {
        istart = iK;

        vcs_VolPhase* volP = VPhaseList[iphase];
        const Cantera::ThermoPhase* tp = volP->ptrThermoPhase();
        string phaseName = volP->PhaseName;
        size_t nSpeciesPhase = volP->nSpecies();
        volP->sendToVCS_VolPM(VCS_DATA_PTR(volPM));
        double TMolesPhase = volP->totalMoles();
        //AssertTrace(TMolesPhase == m_mix->phaseMoles(iphase));
        activity.resize(nSpeciesPhase, 0.0);
        ac.resize(nSpeciesPhase, 0.0);

        mu0.resize(nSpeciesPhase, 0.0);
        mu.resize(nSpeciesPhase, 0.0);
        volPM.resize(nSpeciesPhase, 0.0);
        molalities.resize(nSpeciesPhase, 0.0);

        int actConvention = tp->activityConvention();
        tp->getActivities(VCS_DATA_PTR(activity));
        tp->getActivityCoefficients(VCS_DATA_PTR(ac));
        tp->getStandardChemPotentials(VCS_DATA_PTR(mu0));

        tp->getPartialMolarVolumes(VCS_DATA_PTR(volPM));
        tp->getChemPotentials(VCS_DATA_PTR(mu));
        double VolPhaseVolumes = 0.0;
        for (k = 0; k < nSpeciesPhase; k++) {
            VolPhaseVolumes += volPM[k] * mf[istart + k];
        }
        VolPhaseVolumes *= TMolesPhase;
        vol += VolPhaseVolumes;


        if (actConvention == 1) {
            const Cantera::MolalityVPSSTP* mTP = static_cast<const Cantera::MolalityVPSSTP*>(tp);
            tp->getChemPotentials(VCS_DATA_PTR(mu));
            mTP->getMolalities(VCS_DATA_PTR(molalities));
            tp->getChemPotentials(VCS_DATA_PTR(mu));

            if (iphase == 0) {
                fprintf(FP,"        Name,      Phase,  PhaseMoles,  Mole_Fract, "
                        "Molalities,  ActCoeff,   Activity,"
                        "ChemPot_SS0,   ChemPot,   mole_num,       PMVol, Phase_Volume\n");

                fprintf(FP,"            ,           ,      (kmol),            , "
                        "          ,          ,           ,"
                        "   (J/kmol),  (J/kmol),     (kmol), (m**3/kmol),     (m**3)\n");
            }
            for (k = 0; k < nSpeciesPhase; k++) {
                sName = tp->speciesName(k);
                fprintf(FP,"%12s, %11s, %11.3e, %11.3e, %11.3e, %11.3e, %11.3e,"
                        "%11.3e, %11.3e, %11.3e, %11.3e, %11.3e\n",
                        sName.c_str(),
                        phaseName.c_str(), TMolesPhase,
                        mf[istart + k], molalities[k], ac[k], activity[k],
                        mu0[k]*1.0E-6, mu[k]*1.0E-6,
                        mf[istart + k] * TMolesPhase,
                        volPM[k],  VolPhaseVolumes);
            }

        } else {
            if (iphase == 0) {
                fprintf(FP,"        Name,       Phase,  PhaseMoles,  Mole_Fract,  "
                        "Molalities,   ActCoeff,    Activity,"
                        "  ChemPotSS0,     ChemPot,   mole_num,       PMVol, Phase_Volume\n");

                fprintf(FP,"            ,            ,      (kmol),            ,  "
                        "          ,           ,            ,"
                        "    (J/kmol),    (J/kmol),     (kmol), (m**3/kmol),       (m**3)\n");
            }
            for (k = 0; k < nSpeciesPhase; k++) {
                molalities[k] = 0.0;
            }
            for (k = 0; k < nSpeciesPhase; k++) {
                sName = tp->speciesName(k);
                fprintf(FP,"%12s, %11s, %11.3e, %11.3e, %11.3e, %11.3e, %11.3e, "
                        "%11.3e, %11.3e,% 11.3e, %11.3e, %11.3e\n",
                        sName.c_str(),
                        phaseName.c_str(), TMolesPhase,
                        mf[istart + k],  molalities[k], ac[k],
                        activity[k], mu0[k]*1.0E-6, mu[k]*1.0E-6,
                        mf[istart + k] * TMolesPhase,
                        volPM[k],  VolPhaseVolumes);
            }
        }

#ifdef DEBUG_MODE
        /*
         * Check consistency: These should be equal
         */
        tp->getChemPotentials(VCS_DATA_PTR(m_gibbsSpecies)+istart);
        for (k = 0; k < nSpeciesPhase; k++) {
            if (!vcs_doubleEqual(m_gibbsSpecies[istart+k], mu[k])) {
                fprintf(FP,"ERROR: incompatibility!\n");
                fclose(FP);
                plogf("ERROR: incompatibility!\n");
                exit(EXIT_FAILURE);
            }
        }
#endif
        iK += nSpeciesPhase;
    }
    fclose(FP);
}

void VCS_PROB::setDebugPrintLvl(int lvl)
{
    vcs_debug_print_lvl = lvl;
}

}
