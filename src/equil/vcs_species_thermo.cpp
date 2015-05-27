/**
 * @file vcs_species_thermo.cpp Implementation for the VCS_SPECIES_THERMO
 *   object.
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

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
    SSStar_Vol0(-1.0),
    UseCanteraCalls(false),
    m_VCS_UnitsFormat(VCS_UNITS_UNITLESS)
{
    SS0_Pref = 1.01325E5;
}

VCS_SPECIES_THERMO::VCS_SPECIES_THERMO(const VCS_SPECIES_THERMO& b) :
    IndexPhase(b.IndexPhase),
    IndexSpeciesPhase(b.IndexSpeciesPhase),
    OwningPhase(b.OwningPhase),
    SS0_Model(b.SS0_Model),
    SS0_feSave(b.SS0_feSave),
    SS0_TSave(b.SS0_TSave),
    SS0_T0(b.SS0_T0),
    SS0_H0(b.SS0_H0),
    SS0_S0(b.SS0_S0),
    SS0_Cp0(b.SS0_Cp0),
    SS0_Pref(b.SS0_Pref),
    SSStar_Model(b.SSStar_Model),
    SSStar_Vol_Model(b.SSStar_Vol_Model),
    SSStar_Vol0(b.SSStar_Vol0),
    UseCanteraCalls(b.UseCanteraCalls),
    m_VCS_UnitsFormat(b.m_VCS_UnitsFormat)
{
}

VCS_SPECIES_THERMO&
VCS_SPECIES_THERMO::operator=(const VCS_SPECIES_THERMO& b)
{
    if (&b != this) {
        IndexPhase            = b.IndexPhase;
        IndexSpeciesPhase     = b.IndexSpeciesPhase;
        OwningPhase           = b.OwningPhase;
        SS0_Model             = b.SS0_Model;
        SS0_feSave            = b.SS0_feSave;
        SS0_TSave             = b.SS0_TSave;
        SS0_T0                = b.SS0_T0;
        SS0_H0                = b.SS0_H0;
        SS0_S0                = b.SS0_S0;
        SS0_Cp0               = b.SS0_Cp0;
        SS0_Pref              = b.SS0_Pref;
        SSStar_Model          = b.SSStar_Model;
        SSStar_Vol_Model      = b.SSStar_Vol_Model;
        SSStar_Vol0           = b.SSStar_Vol0;
        UseCanteraCalls       = b.UseCanteraCalls;
        m_VCS_UnitsFormat     = b.m_VCS_UnitsFormat;
    }
    return *this;
}

VCS_SPECIES_THERMO* VCS_SPECIES_THERMO::duplMyselfAsVCS_SPECIES_THERMO()
{
    return new VCS_SPECIES_THERMO(*this);
}

double VCS_SPECIES_THERMO::GStar_R_calc(size_t kglob, double TKelvin,
                                        double pres)
{
    double fe = G0_R_calc(kglob, TKelvin);
    double T = TKelvin;
    if (UseCanteraCalls) {
        if (m_VCS_UnitsFormat != VCS_UNITS_MKS) {
            throw CanteraError("VCS_SPECIES_THERMO::GStar_R_calc",
                               "Possible inconsistency");
        }
        size_t kspec = IndexSpeciesPhase;
        OwningPhase->setState_TP(TKelvin, pres);
        fe = OwningPhase->GStar_calc_one(kspec);
        double R = vcsUtil_gasConstant(m_VCS_UnitsFormat);
        fe /= R;
    } else {
        double pref = SS0_Pref;
        switch (SSStar_Model) {
        case VCS_SSSTAR_CONSTANT:
            break;
        case VCS_SSSTAR_IDEAL_GAS:
            fe += T * log(pres/ pref);
            break;
        default:
            throw CanteraError("VCS_SPECIES_THERMO::GStar_R_calc",
                               "unknown SSStar model");
        }
    }
    return fe;
}

double VCS_SPECIES_THERMO::VolStar_calc(size_t kglob, double TKelvin,
                                        double presPA)
{
    double vol;

    double T = TKelvin;
    if (UseCanteraCalls) {
        if (m_VCS_UnitsFormat != VCS_UNITS_MKS) {
            throw CanteraError("VCS_SPECIES_THERMO::VolStar_calc",
                               "Possible inconsistency");
        }
        size_t kspec = IndexSpeciesPhase;
        OwningPhase->setState_TP(TKelvin, presPA);
        vol = OwningPhase->VolStar_calc_one(kspec);
    } else {
        switch (SSStar_Vol_Model) {
        case VCS_SSVOL_CONSTANT:
            vol = SSStar_Vol0;
            break;
        case VCS_SSVOL_IDEALGAS:
            vol= GasConstant * T / presPA;
            break;
        default:
            throw CanteraError("VCS_SPECIES_THERMO::VolStar_calc",
                               "unknown SSVol model");
        }
    }
    return vol;
}

double VCS_SPECIES_THERMO::G0_R_calc(size_t kglob, double TKelvin)
{
    double fe, H, S;
    if (SS0_Model == VCS_SS0_CONSTANT) {
        return SS0_feSave;
    }
    if (TKelvin == SS0_TSave) {
        return SS0_feSave;
    }
    if (UseCanteraCalls) {
        if (m_VCS_UnitsFormat != VCS_UNITS_MKS) {
            throw CanteraError("VCS_SPECIES_THERMO::G0_R_calc",
                               "Possible inconsistency");
        }
        size_t kspec = IndexSpeciesPhase;
        OwningPhase->setState_T(TKelvin);
        fe = OwningPhase->G0_calc_one(kspec);
        double R = vcsUtil_gasConstant(m_VCS_UnitsFormat);
        fe /= R;
    } else {
        switch (SS0_Model) {
        case VCS_SS0_CONSTANT:
            fe = SS0_feSave;
            break;
        case VCS_SS0_CONSTANT_CP:
            H  = SS0_H0 + (TKelvin - SS0_T0) * SS0_Cp0;
            S  = SS0_Cp0 + SS0_Cp0 * log((TKelvin / SS0_T0));
            fe = H - TKelvin * S;
            break;
        default:
            throw CanteraError("VCS_SPECIES_THERMO::G0_R_calc",
                               "unknown model");
        }
    }
    SS0_feSave = fe;
    SS0_TSave = TKelvin;
    return fe;
}

double VCS_SPECIES_THERMO::eval_ac(size_t kglob)
{
    double ac;
    /*
     *  Activity coefficients are frequently evaluated on a per phase
     *  basis. If they are, then the currPhAC[] boolean may be used
     *  to reduce repeated work. Just set currPhAC[iph], when the
     *  activity coefficients for all species in the phase are reevaluated.
     */
    if (UseCanteraCalls) {
        size_t kspec = IndexSpeciesPhase;
        ac = OwningPhase->AC_calc_one(kspec);
    } else {
        ac = 1.0;
    }
    return ac;
}

}
