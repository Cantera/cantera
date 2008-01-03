/*
 * $Id$
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vcs_solve.h"
#include "vcs_species_thermo.h"
#include "vcs_defs.h"
#include "vcs_VolPhase.h"
#include "vcs_nasa_poly.h"
#include "vcs_Exception.h"
#include "vcs_internal.h"

using namespace std;

namespace VCSnonideal {

/*****************************************************************************
 *
 * constructor():
 */

VCS_SPECIES_THERMO::VCS_SPECIES_THERMO(int indexPhase, 
				       int indexSpeciesPhase) :
    
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
  SS0_Pref(1.0),
  SS0_Params(0),
  SSStar_Model(VCS_SSSTAR_CONSTANT),
  SSStar_Params(0),
  Activity_Coeff_Model(VCS_AC_CONSTANT),
  Activity_Coeff_Params(0),
  SSStar_Vol_Model(VCS_SSVOL_IDEALGAS),
  SSStar_Vol_Params(0),
  SSStar_Vol0(-1.0),
  UseCanteraCalls(false),
  m_VCS_UnitsFormat(VCS_UNITS_UNITLESS)
{
  /*
   * Set up the numerical value for P_reference, based on the current
   * global units choice.
   */
  if (m_VCS_UnitsFormat == VCS_UNITS_MKS) {
    SS0_Pref = 1.0133E5;
  } else {
    SS0_Pref = 1.0;
  }
}


/******************************************************************************
 *
 * destructor
 */
VCS_SPECIES_THERMO::~VCS_SPECIES_THERMO() 
{
  if (SS0_Model  == VCS_SS0_NASA_POLY) {
    vcs_nasa_poly_destroy((VCS_NASA_POLY **) &(this->SS0_Params));
    SS0_Params = 0;
  }
}

/*****************************************************************************
 *
 * Copy Constructor VCS_SPECIES_THERMO
 */
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
  SS0_Params(0),
  SSStar_Model(b.SSStar_Model),
  SSStar_Params(0),
  Activity_Coeff_Model(b.Activity_Coeff_Model),
  Activity_Coeff_Params(0),
  SSStar_Vol_Model(b.SSStar_Vol_Model),
  SSStar_Vol_Params(0),
  SSStar_Vol0(b.SSStar_Vol0),
  UseCanteraCalls(b.UseCanteraCalls),
  m_VCS_UnitsFormat(b.m_VCS_UnitsFormat)
{
	VCS_NASA_POLY *ppp = 0;	
  switch (SS0_Model) {
  case VCS_SS0_NASA_POLY:
    ppp = (VCS_NASA_POLY *) b.SS0_Params;
    SS0_Params = (void *) new VCS_NASA_POLY(*ppp);
    break;
  default:
	  ppp = 0;
	  SS0_Params = 0;
    break;
  }
}

/*****************************************************************************
 *
 * Assignment operator for VCS_SPECIES_THERMO
 */
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

    VCS_NASA_POLY *ppp= 0;
    switch (SS0_Model) {
    case VCS_SS0_NASA_POLY:
      ppp = (VCS_NASA_POLY *) b.SS0_Params;
      SS0_Params = (void *) new VCS_NASA_POLY(*ppp);
      break;
    default:
      break;
    }   

    SSStar_Model          = b.SSStar_Model;
    /*
     * shallow copy because function is undeveloped.
     */
    SSStar_Params         = b.SSStar_Params;
    Activity_Coeff_Model  = b.Activity_Coeff_Model;
    /*
     * shallow copy because function is undeveloped.
     */
    Activity_Coeff_Params = b.Activity_Coeff_Params;
    SSStar_Vol_Model      = b.SSStar_Vol_Model;
    /*
     * shallow copy because function is undeveloped.
     */
    SSStar_Vol_Params     = b.SSStar_Vol_Params; 
    SSStar_Vol0           = b.SSStar_Vol0;
    UseCanteraCalls       = b.UseCanteraCalls;
    m_VCS_UnitsFormat     = b.m_VCS_UnitsFormat;
  }
  return *this;
}

/******************************************************************************
 *
 * duplMyselfAsVCS_SPECIES_THERMO():                (virtual)
 *
 *    This routine can duplicate inherited objects given a base class
 *    pointer. It relies on valid copy constructors.
 */

VCS_SPECIES_THERMO* VCS_SPECIES_THERMO::duplMyselfAsVCS_SPECIES_THERMO() {
  VCS_SPECIES_THERMO* ptr = new VCS_SPECIES_THERMO(*this);
  return  ptr;
}


/**************************************************************************
 *
 * GStar_R_calc();
 *
 *  This function calculates the standard state Gibbs free energy
 *  for species, kspec, at the solution temperature TKelvin and
 *  solution pressure, Pres.
 *  
 *
 *  Input
 *   kglob = species global index.
 *   TKelvin = Temperature in Kelvin
 *   pres = pressure is given in units specified by if__ variable.
 *
 *
 * Output
 *    return value = standard state free energy in units of Kelvin.
 */
double VCS_SPECIES_THERMO::GStar_R_calc(int kglob, double TKelvin, 
					double pres)
{
  char yo[] = "VCS_SPECIES_THERMO::GStar_R_calc ";
  double fe, T;
  fe = G0_R_calc(kglob, TKelvin);
  T = TKelvin;
  if (UseCanteraCalls) {
    AssertThrowVCS(m_VCS_UnitsFormat == VCS_UNITS_MKS, "Possible inconsistency");
    int kspec = IndexSpeciesPhase;
    fe = OwningPhase->GStar_calc_one(kspec, TKelvin, pres);
    double R = vcsUtil_gasConstant(m_VCS_UnitsFormat);
    fe /= R;
  } else {
    double pref = SS0_Pref;
    switch(SSStar_Model) {
    case VCS_SSSTAR_CONSTANT:
      break;
    case VCS_SSSTAR_IDEAL_GAS:
      fe += T * log( pres/ pref );	 
      break;
    default:
      plogf("%sERROR: unknown SSStar model\n", yo);
      exit(-1);
    }
  }
  return fe;
}
   
/**************************************************************************
 *
 * VolStar_calc:
 *
 *  This function calculates the standard state molar volume
 *  for species, kspec, at the temperature TKelvin and pressure, Pres,
 * 
 *  Input
 *
 * Output
 *    return value = standard state volume in    cm**3 per mol.
 *                   (VCS_UNITS_MKS)              m**3 / kmol
 */
double VCS_SPECIES_THERMO::
VolStar_calc(int kglob, double TKelvin, double pres)
{
  char yo[] = "VCS_SPECIES_THERMO::VStar_calc ";
  double vol, T;
   
  T = TKelvin;
  if (UseCanteraCalls) {
    AssertThrowVCS(m_VCS_UnitsFormat == VCS_UNITS_MKS, "Possible inconsistency");
    int kspec = IndexSpeciesPhase;
    vol = OwningPhase->VolStar_calc_one(kspec, TKelvin, pres);
  } else {
    switch(SSStar_Vol_Model) {
    case VCS_SSVOL_CONSTANT:
      vol = SSStar_Vol0;
      break;
    case VCS_SSVOL_IDEALGAS:
      if (m_VCS_UnitsFormat == VCS_UNITS_MKS) {
	vol = 8.31451E3 * T / pres;
      } else {
	vol= 83.14510 / 1.01325 * T / pres;	 
      }
      break;
    default:     
      plogf("%sERROR: unknown SSVol model\n", yo);
      exit(-1);
    } 
  }
  return vol;
} 

/**************************************************************************
 *
 * G0_R_calc:
 *
 *  This function calculates the naught state Gibbs free energy
 *  for species, kspec, at the temperature TKelvin
 *
 *  Input
 *   kglob = species global index.
 *   TKelvin = Temperature in Kelvin
 *
 * Output
 *    return value = naught state free energy in Kelvin.
 */
double VCS_SPECIES_THERMO::G0_R_calc(int kglob, double TKelvin)
{
#ifdef DEBUG
  char yo[] = "VS_SPECIES_THERMO::G0_R_calc ";
#endif
  double fe, H, S;
  if (SS0_Model == VCS_SS0_CONSTANT) {
    fe = SS0_feSave;  
    return fe;
  }
  if (TKelvin == SS0_TSave) {
    fe = SS0_feSave;
    return fe;
  }
  if (UseCanteraCalls) {
    AssertThrowVCS(m_VCS_UnitsFormat == VCS_UNITS_MKS, "Possible inconsistency");
    int kspec = IndexSpeciesPhase;
    fe = OwningPhase->G0_calc_one(kspec, TKelvin);
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
    case VCS_SS0_NASA_POLY: 
      fe = vcs_G0_NASA(TKelvin, (VCS_NASA_POLY *) SS0_Params);
      break;
    default:
#ifdef DEBUG
      plogf("%sERROR: unknown model\n", yo);
#endif
      exit(-1);
    }
  }
  SS0_feSave = fe;
  SS0_TSave = TKelvin;
  return fe;
} 

/**************************************************************************
 *
 * eval_ac:
 *
 *  This function evaluates the activity coefficient
 *  for species, kspec
 *
 *  Input
 *      kglob -> integer value of the species in the global 
 *            species list within VCS_GLOB. Phase and local species id
 *             can be looked up within object.
 * 
 *   Note, T, P and mole fractions are obtained from the
 *   single private instance of VCS_GLOB
 *   
 *
 * Output
 *    return value = activity coefficient for species kspec
 */
double VCS_SPECIES_THERMO::eval_ac(int kglob)
{
#ifdef DEBUG
  char yo[] = "VCS_SPECIES_THERMO::eval_ac ";
#endif
  double ac;
  /*
   *  Activity coefficients are frequently evaluated on a per phase
   *  basis. If they are, then the currPhAC[] boolean may be used
   *  to reduce repeated work. Just set currPhAC[iph], when the 
   *  activity coefficients for all species in the phase are reevaluated.
   */
  if (UseCanteraCalls) {
    int kspec = IndexSpeciesPhase;
    ac = OwningPhase->AC_calc_one(kspec);
  } else {
    switch (Activity_Coeff_Model) {
    case VCS_AC_CONSTANT:
      ac = 1.0;
      break;
    case VCS_AC_DEBYE_HUCKEL:
	 
      plogf("Not implemented Yet\n");
      exit(-1);
      break;
	 
    case VCS_AC_REGULAR_SOLN:
	 
      plogf("Not implemented Yet\n");
      exit(-1);
      break;
	 
    case VCS_AC_MARGULES:
	 
      plogf("Not implemented Yet\n");
      exit(-1);
      break;
    default:
#ifdef DEBUG
      plogf("%sERROR: unknown model\n", yo);
#endif
      exit(-1);
    }
  }
  return ac;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double VCS_SOLVE::vcs_Gxs_phase_calc(vcs_VolPhase *Vphase, double *mf_PO)
   
  /**************************************************************************
   *
   *  vcs_Gxs_calc:
   *
   *  This function evaluates the Gibbs Excess free energy function for
   *  the phase pointed to by Vphase.
   *
   *  There are two ways. They may be evaluated from the 
   *  activity coefficients themselves
   *
   *    Gxs/RT = sum_i_inphase( X_i * ln (ActCoeff_i)) 
   *
   *  Or, the actual formulas for the excess Gibbs free energy may
   *  be used (which the activity coefficients probably came from anyway.
   *
   *  Input
   *   phase_ptr => Pointer to the phase that we want to calculate
   *                the 
   *   mf_PO     => Vector of mole fractions in the phase
   *                in "Phase Order" order. This must sum to one. However
   *                this condition is not checked.
   *
   * Output
   *    return value = activity coefficient for species kspec
   ***************************************************************************/
{
  int kspec, kglob;
  double Gxs = 0.0, ac;
  VCS_SPECIES_THERMO *ts_ptr;
  if (Vphase->Activity_Coeff_Model != VCS_AC_CONSTANT) {
    for (kspec = 0; kspec < Vphase->NVolSpecies; kspec++) {
      kglob = Vphase->IndSpecies[kspec];
      ts_ptr = SpeciesThermo[kglob];
      ac = ts_ptr->eval_ac(kspec);
      Gxs += mf_PO[kspec] * log(ac);
    }
  }
  return Gxs;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double VCS_SOLVE::vcs_Gxs_calc(int iphase)

  /**************************************************************************
   *
   *  vcs_Gxs_calc:
   *
   *    This function evaluates the Gibbs Excess free energy function.
   *
   *  There are two ways. They may be evaluated from the
   *  activity coefficients themselves
   *
   *    Gxs/RT = sum_i_inphase( X_i * ln (ActCoeff_i))
   *
   *  Or, the actual formulas for the excess Gibbs free energy may
   *  be used (which the activity coefficients probably came from anyway.
   *
   *  Input
   *
   *
   * Output
   *    return value = activity coefficient for species kspec
   ***************************************************************************/
{
  int kspec;
  double Gxs = 0.0, ac;
  double totmol = TPhMoles[iphase];
  vcs_VolPhase *Vphase = VPhaseList[iphase];
  VCS_SPECIES_THERMO *ts_ptr;

  if (totmol != 0.0 && Vphase->Activity_Coeff_Model != VCS_AC_CONSTANT) {
    for (kspec = 0; kspec < m_numSpeciesRdc; kspec++) {
      if (PhaseID[kspec] == iphase) {
	if (SpeciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	  ts_ptr = SpeciesThermo[kspec];
	  ac = ts_ptr->eval_ac(kspec);
	  Gxs += soln[kspec]/totmol * log(ac);
	} else {
	  plogf("FILL IN\n");
	  exit(-1);
	}
      }
    }
  }
  return Gxs;
}
/*****************************************************************************/
}
