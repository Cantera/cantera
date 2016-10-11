/**
 * @file vcs_species_thermo.cpp Implementation for the VCS_SPECIES_THERMO
 *   object.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/equil/vcs_species_thermo.h"
#include "cantera/equil/vcs_defs.h"
#include "cantera/equil/vcs_VolPhase.h"

#include "cantera/base/ctexceptions.h"
#include "cantera/equil/vcs_internal.h"

using namespace std;
namespace Cantera
{
VCS_SPECIES_THERMO::VCS_SPECIES_THERMO(size_t indexPhase,
                                       size_t indexSpeciesPhase) :
    IndexPhase(indexPhase),
    IndexSpeciesPhase(indexSpeciesPhase),
    OwningPhase(0),
    SS0_Model(VCS_SS0_CONSTANT),
    SS0_feSave(0.0),
    SS0_TSave(-90.0),
    SS0_T0(273.15),
    SS0_H0(0.0),
    SS0_S0(0.0),
    SS0_Cp0(0.0),
    SS0_Pref(1.01325E5),
    SSStar_Model(VCS_SSSTAR_CONSTANT),
    SSStar_Vol_Model(VCS_SSVOL_IDEALGAS),
    SSStar_Vol0(-1.0)
{
    SS0_Pref = 1.01325E5;
}

VCS_SPECIES_THERMO* VCS_SPECIES_THERMO::duplMyselfAsVCS_SPECIES_THERMO()
{
    return new VCS_SPECIES_THERMO(*this);
}

double VCS_SPECIES_THERMO::GStar_R_calc(size_t kglob, double TKelvin,
                                        double pres)
{
    warn_deprecated("VCS_SPECIES_THERMO::GStar_R_calc",
        "Unused. To be removed after Cantera 2.3.");
    double fe = G0_R_calc(kglob, TKelvin);
    OwningPhase->setState_TP(TKelvin, pres);
    fe = OwningPhase->GStar_calc_one(IndexSpeciesPhase);
    return fe / GasConstant;
}

double VCS_SPECIES_THERMO::VolStar_calc(size_t kglob, double TKelvin,
                                        double presPA)
{
    warn_deprecated("VCS_SPECIES_THERMO::VolStar_calc",
        "Unused. To be removed after Cantera 2.3.");
    OwningPhase->setState_TP(TKelvin, presPA);
    return OwningPhase->VolStar_calc_one(IndexSpeciesPhase);
}

double VCS_SPECIES_THERMO::G0_R_calc(size_t kglob, double TKelvin)
{
    warn_deprecated("VCS_SPECIES_THERMO::G0_R_calc",
        "Unused. To be removed after Cantera 2.3.");
    if (SS0_Model == VCS_SS0_CONSTANT) {
        return SS0_feSave;
    }
    if (TKelvin == SS0_TSave) {
        return SS0_feSave;
    }
    OwningPhase->setState_T(TKelvin);
    double fe = OwningPhase->G0_calc_one(IndexSpeciesPhase);
    fe /= GasConstant;
    SS0_feSave = fe;
    SS0_TSave = TKelvin;
    return fe;
}

double VCS_SPECIES_THERMO::eval_ac(size_t kglob)
{
    // Activity coefficients are frequently evaluated on a per phase basis. If
    // they are, then the currPhAC[] boolean may be used to reduce repeated
    // work. Just set currPhAC[iph], when the activity coefficients for all
    // species in the phase are reevaluated.
    warn_deprecated("VCS_SPECIES_THERMO::eval_ac",
                    "Unused. To be removed after Cantera 2.3.");
    return OwningPhase->AC_calc_one(IndexSpeciesPhase);
}

}
