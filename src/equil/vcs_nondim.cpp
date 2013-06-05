/**
 *  @file vcs_nondim.cpp
 *     Nondimensionalization routines within VCSnonideal
 */
/*
 * Copyright (2007) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/equil/vcs_solve.h"
#include "cantera/equil/vcs_internal.h"
#include "cantera/equil/vcs_VolPhase.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"

namespace VCSnonideal
{
double VCS_SOLVE::vcs_nondim_Farad(int mu_units, double TKelvin) const
{
    double Farad;
    if (TKelvin <= 0.0) {
        TKelvin = 293.15;
    }
    switch (mu_units) {
    case VCS_UNITS_MKS:
    case VCS_UNITS_KJMOL:
    case VCS_UNITS_KCALMOL:
        Farad = Cantera::ElectronCharge * Cantera::Avogadro /
                (TKelvin * Cantera::GasConstant);
        break;
    case VCS_UNITS_UNITLESS:
        Farad = Cantera::ElectronCharge * Cantera::Avogadro;
        break;
    case VCS_UNITS_KELVIN:
        Farad = Cantera::ElectronCharge * Cantera::Avogadro/ TKelvin;
        break;
    default:
        plogf("vcs_nondim_Farad error: unknown units: %d\n", mu_units);
        plogendl();
        exit(EXIT_FAILURE);
    }
    return Farad;
}

double VCS_SOLVE::vcs_nondimMult_TP(int mu_units, double TKelvin) const
{
    double rt;
    if (TKelvin <= 0.0) {
        TKelvin = 293.15;
    }
    switch (mu_units) {
    case VCS_UNITS_KCALMOL:
        rt = TKelvin * Cantera::GasConst_cal_mol_K * 1e-3;
        break;
    case VCS_UNITS_UNITLESS:
        rt = 1.0;
        break;
    case VCS_UNITS_KJMOL:
        rt = TKelvin * Cantera::GasConstant * 1e-6;
        break;
    case VCS_UNITS_KELVIN:
        rt = TKelvin;
        break;
    case VCS_UNITS_MKS:
        rt = TKelvin * Cantera::GasConstant;
        break;
    default:
        plogf("vcs_nondimMult_TP error: unknown units: %d\n", mu_units);
        plogendl();
        exit(EXIT_FAILURE);
    }
    return rt;
}

void VCS_SOLVE::vcs_nondim_TP()
{
    double tf;
    if (m_unitsState == VCS_DIMENSIONAL_G) {
        m_unitsState = VCS_NONDIMENSIONAL_G;
        tf = 1.0 / vcs_nondimMult_TP(m_VCS_UnitsFormat, m_temperature);
        for (size_t i = 0; i < m_numSpeciesTot; ++i) {
            /*
             *        Modify the standard state and total chemical potential data,
             *        FF(I),  to make it dimensionless, i.e.,  mu / RT.
             *        Thus, we may divide it by the temperature.
             */
            m_SSfeSpecies[i] *= tf;
            m_deltaGRxn_new[i] *= tf;
            m_deltaGRxn_old[i] *= tf;
            m_feSpecies_old[i] *= tf;
        }

        m_Faraday_dim =  vcs_nondim_Farad(m_VCS_UnitsFormat, m_temperature);

        /*
         * Scale the total moles if necessary:
         *  First find out the total moles
         */
        double tmole_orig = vcs_tmoles();

        /*
         * Then add in the total moles of elements that are goals. Either one
         * or the other is specified here.
         */
        double esum = 0.0;
        for (size_t i = 0; i < m_numElemConstraints; ++i) {
            if (m_elType[i] == VCS_ELEM_TYPE_ABSPOS) {
                esum += fabs(m_elemAbundancesGoal[i]);
            }
        }
        tmole_orig += esum;

        /*
         * Ok now test out the bounds on the total moles that this program can
         * handle. These are a bit arbitrary. However, it would seem that any
         * reasonable input would be between these two numbers below.
         */
        if (tmole_orig < 1.0E-200 || tmole_orig > 1.0E200) {
            plogf(" VCS_SOLVE::vcs_nondim_TP ERROR: Total input moles , %g,  is outside the range handled by vcs. exit",
                  tmole_orig);
            plogendl();
            throw Cantera::CanteraError("VCS_SOLVE::vcs_nondim_TP",
                                        " Total input moles ," + Cantera::fp2str(tmole_orig) +
                                        "is outside the range handled by vcs.\n");
        }

        // Determine the scale of the problem
        if (tmole_orig > 1.0E4) {
            m_totalMoleScale =  tmole_orig / 1.0E4;
        } else if (tmole_orig < 1.0E-4) {
            m_totalMoleScale =  tmole_orig / 1.0E-4;
        } else {
            m_totalMoleScale = 1.0;
        }

        if (m_totalMoleScale != 1.0) {
            if (m_VCS_UnitsFormat == VCS_UNITS_MKS) {
#ifdef DEBUG_MODE
                if (m_debug_print_lvl >= 2) {
                    plogf("  --- vcs_nondim_TP() called: USING A MOLE SCALE OF %g until further notice", m_totalMoleScale);
                    plogendl();
                }
#endif
                for (size_t i = 0; i < m_numSpeciesTot; ++i) {
                    if (m_speciesUnknownType[i] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                        m_molNumSpecies_old[i] *= (1.0 / m_totalMoleScale);
                    }
                }
                for (size_t i = 0; i < m_numElemConstraints; ++i) {
                    m_elemAbundancesGoal[i] *= (1.0 / m_totalMoleScale);
                }

                for (size_t iph = 0; iph < m_numPhases; iph++) {
                    TPhInertMoles[iph] *= (1.0 / m_totalMoleScale);
                    if (TPhInertMoles[iph] != 0.0) {
                        vcs_VolPhase* vphase = m_VolPhaseList[iph];
                        vphase->setTotalMolesInert(TPhInertMoles[iph]);
                    }
                }
            }
            vcs_tmoles();
        }
    }
}

void VCS_SOLVE::vcs_redim_TP(void)
{
    double tf;
    if (m_unitsState != VCS_DIMENSIONAL_G) {
        m_unitsState = VCS_DIMENSIONAL_G;
        tf = vcs_nondimMult_TP(m_VCS_UnitsFormat, m_temperature);
        for (size_t i = 0; i < m_numSpeciesTot; ++i) {
            /*
             *        Modify the standard state and total chemical potential data,
             *        FF(I),  to make it have units, i.e. mu = RT * mu_star
             */
            m_SSfeSpecies[i] *= tf;
            m_deltaGRxn_new[i] *= tf;
            m_deltaGRxn_old[i] *= tf;
            m_feSpecies_old[i] *= tf;
        }
        m_Faraday_dim *= tf;
    }
    if (m_totalMoleScale != 1.0) {
        if (m_VCS_UnitsFormat == VCS_UNITS_MKS) {
#ifdef DEBUG_MODE
            if (m_debug_print_lvl >= 2) {
                plogf("  --- vcs_redim_TP() called: getting rid of mole scale of %g", m_totalMoleScale);
                plogendl();
            }
#endif
            for (size_t i = 0; i < m_numSpeciesTot; ++i) {
                if (m_speciesUnknownType[i] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                    m_molNumSpecies_old[i] *= m_totalMoleScale;
                }
            }
            for (size_t i = 0; i < m_numElemConstraints; ++i) {
                m_elemAbundancesGoal[i] *= m_totalMoleScale;
            }

            for (size_t iph = 0; iph < m_numPhases; iph++) {
                TPhInertMoles[iph] *= m_totalMoleScale;
                if (TPhInertMoles[iph] != 0.0) {
                    vcs_VolPhase* vphase = m_VolPhaseList[iph];
                    vphase->setTotalMolesInert(TPhInertMoles[iph]);
                }
            }
            vcs_tmoles();
        }
    }
}

void VCS_SOLVE::vcs_printChemPotUnits(int unitsFormat) const
{
    switch (unitsFormat) {
    case VCS_UNITS_KCALMOL:
        plogf("kcal/gmol");
        break;
    case VCS_UNITS_UNITLESS:
        plogf("dimensionless");
        break;
    case VCS_UNITS_KJMOL:
        plogf("kJ/gmol");
        break;
    case VCS_UNITS_KELVIN:
        plogf("Kelvin");
        break;
    case VCS_UNITS_MKS:
        plogf("J/kmol");
        break;
    default:
        plogf("unknown units!");
        exit(EXIT_FAILURE);
    }
}

}
