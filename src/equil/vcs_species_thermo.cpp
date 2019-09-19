/**
 * @file vcs_species_thermo.cpp Implementation for the VCS_SPECIES_THERMO
 *   object.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/equil/vcs_species_thermo.h"
#include "cantera/equil/vcs_defs.h"
#include "cantera/equil/vcs_VolPhase.h"

#include "cantera/base/ctexceptions.h"
#include "cantera/equil/vcs_internal.h"

using namespace std;
namespace Cantera
{
VCS_SPECIES_THERMO::VCS_SPECIES_THERMO() :
    IndexPhase(0),
    IndexSpeciesPhase(0),
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
}

}
