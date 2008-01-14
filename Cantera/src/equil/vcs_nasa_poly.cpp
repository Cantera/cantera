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
#include <string.h>

#include "vcs_defs.h"
#include "vcs_nasa_poly.h"
#include "vcs_species_thermo.h" 
#include "vcs_internal.h"

#ifdef WIN32
#pragma warning(disable:4996)
#endif

namespace VCSnonideal {

/******************************************************************************
 *
 *  Constructor
 */
VCS_NASA_POLY::VCS_NASA_POLY(int numTempRegions, int numEl) :
  NumTempRegions(numTempRegions),
  NumEl(numEl),
  PhType(' ')
{
  Date[0]   = '\0';
  SpName[0] = '\0';   
  PhName[0] = '\0';
  ElName.resize(numEl, "");
  ElComp.resize(numEl, 0.0);
  
  if (NumTempRegions < 1) NumTempRegions = 1;
  Tlimits.resize(NumTempRegions+1, 0.0);
  Acoeff.resize(NumTempRegions, 7, 0.0);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

VCS_NASA_POLY *vcs_nasa_poly_create(int numTempRegions, int numEl)
   
  /**************************************************************************
   *
   * vcs_nasa_poly_create:
   *
   *    Constructor routine for the nasa polynomial structure.
   *    It initializes all data to zero. The number of temperature regions
   *    malloced is storred within the structure itself.
   *
   * Input
   *     numTempRegions: Number of temperature regions
   *
   * Return
   *     Pointer to the newly malloced structure.
   *     If NULL, then an out of memory condition occurred
   ***************************************************************************/
   
{
  VCS_NASA_POLY *poly_ptr;   
  poly_ptr = new VCS_NASA_POLY(numTempRegions, numEl);
  return poly_ptr;
}

/***************************************************************************
 * Copy Constructor
 */
VCS_NASA_POLY::VCS_NASA_POLY(const VCS_NASA_POLY &b) :
  NumTempRegions(0),
  NumEl(0)
{
  *this = b;
}
/******************************************************************************
 *
 * operator=()
 *
 */
VCS_NASA_POLY& VCS_NASA_POLY::operator=(const VCS_NASA_POLY &b) {
  if (&b != this) {
    NumTempRegions = b.NumTempRegions;
    Tlimits = b.Tlimits;
    Acoeff = b.Acoeff;
    NumEl = b.NumEl;
    ElComp = b.ElComp;
    ElName = b.ElName;
    strcpy(Date, b.Date);
    PhType = b.PhType;
    strcpy(SpName, b.SpName);
    strcpy(PhName, b.PhName);
  }
  return *this;
}


/*****************************************************************************
 *
 * ~VCS_NASA_POLY():
 *
 *  Destructor for class
 */
VCS_NASA_POLY::~VCS_NASA_POLY() {
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void vcs_nasa_poly_destroy(VCS_NASA_POLY **poly_hdl)
   
  /**************************************************************************
   *
   * vcs_nasa_poly_destroy:
   *
   *    Destructor routine for the nasa polynomial structure.
   *************************************************************************/ 
{
  VCS_NASA_POLY *poly_ptr = *poly_hdl;
  if (poly_ptr) {
    delete poly_ptr;
    poly_ptr = 0;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double vcs_G0_NASA(double TKelvin, VCS_NASA_POLY *poly_ptr)
   
  /**************************************************************************
   *
   * vcs_GibbsFE_NASA:
   *
   *   Calculate the Gibbs free energy (in Kelvin) for a single species 
   *   using the Nasa polynomial format.
   *
   * Input
   *    TKelvin  = Temperature in Kelvin.
   *    poly_ptr = Pointer to structure containing the NASA Polynomial
   *               coefficients   
   *
   * Return
   *    gibbsFE = Gibbs free energy / R -> units of kelvin
   *
   * Error Conditions
   *    VCS_THERMO_OUTOFRANGE:
   *         If the input temperature, is out of range of the polynomials,
   *    an error Flag is set, and the temperature is storred in the 
   *    error structure.
   ***************************************************************************/
{
  int iRegion;
  double *a, gibbsFE;
  double *Tlim = VCS_DATA_PTR(poly_ptr->Tlimits);
  static double Tsave = -10., C0, C1, C2, C3, C4, C5;
  //extern CPC_ERR_STRUCT cpcE;
  /*
   *  Find the temperature region
   */
  if (TKelvin <= *Tlim) {
#ifdef DEBUG_MODE
    plogf("vcs_G0_NASA error: TKelvin below lowest bounds %g\n", *Tlim);
#endif
    iRegion = 0;
  
    if (TKelvin <= 0.0) {
      gibbsFE = poly_ptr->Acoeff[0][5];
      return gibbsFE;
    }
    goto L_FOUNDREGION;
  }
  for (iRegion = 0; iRegion < poly_ptr->NumTempRegions; iRegion++) {
    Tlim++;
    if (TKelvin <= *Tlim) goto L_FOUNDREGION;
  }
#ifdef DEBUG_MODE
  plogf("vcs_G0_NASA error: TKelvin above highest bounds %g\n", *(Tlim));
#endif
 
  iRegion--;
 L_FOUNDREGION:;
  a = poly_ptr->Acoeff[iRegion];
  if (Tsave != TKelvin) {
    Tsave = TKelvin;
    C0 = 1.0 - log(TKelvin);
    C1 = TKelvin * 0.5;
    C2 = TKelvin * TKelvin;
    C3 = C2 * TKelvin;
    C4 = C3 * TKelvin;
    C2 /= 6.0;
    C3 /= 12.0;
    C4 /= 20.0;
    C5 = 1.0 / TKelvin;   
  }
  gibbsFE = a[0]*C0 - a[1]*C1 - a[2]*C2 - a[3]*C3 - a[4]*C4 + a[5]*C5 - a[6];
  gibbsFE *= TKelvin;
  return gibbsFE;
} /***************************************************************************/

double vcs_H0_NASA(double TKelvin, VCS_NASA_POLY *poly_ptr)
   
  /**************************************************************************
   *
   * vcs_H0_NASA:
   *
   *   Calculate the standard state Enthalpy (in Kelvin) for a single species 
   *   using the Nasa polynomial format.
   *
   * Input
   *    TKelvin  = Temperature in Kelvin.
   *    poly_ptr = Pointer to structure containing the NASA Polynomial
   *               coefficients
   *
   * Return
   *     H0 = Standard State Enthalpy / R -> units of kelvin
   *
   * Error Conditions
   *    VCS_THERMO_OUTOFRANGE:
   *         If the input temperature, is out of range of the polynomials,
   *    an error Flag is set, and the temperature is storred in the 
   *    error structure.
   ***************************************************************************/
{
  int iRegion;
  double *a, H0;
  double *Tlim = VCS_DATA_PTR(poly_ptr->Tlimits);
  static double Tsave = -10., C1, C2, C3, C4, C5;
  /*
   *  Find the temperature region
   */
  if (TKelvin <= *Tlim) {
#ifdef DEBUG_MODE
    plogf("vcs_H0_NASA error: TKelvin below lowest bounds\n");
#endif
    iRegion = 0;
    if (TKelvin <= 0.0) {
      H0 = poly_ptr->Acoeff[0][6];
      return H0;
    }
    goto L_FOUNDREGION;
  }
  for (iRegion = 0; iRegion < poly_ptr->NumTempRegions; iRegion++) {
    Tlim++;
    if (TKelvin <= *Tlim) goto L_FOUNDREGION;
  }
#ifdef DEBUG_MODE
  plogf("vcs_H0_NASA error: TKelvin above highest bounds\n");
#endif
  iRegion--;
 L_FOUNDREGION:;
  a = poly_ptr->Acoeff[iRegion];
  if (Tsave != TKelvin) {
    Tsave = TKelvin;
    C1 = TKelvin * 0.5;
    C2 = TKelvin * TKelvin;
    C3 = C2 * TKelvin;
    C4 = C3 * TKelvin;
    C2 /= 3.0;
    C3 /= 4.0;
    C4 /= 5.0;
    C5 = 1.0 / TKelvin;   
  }
  H0 = a[0] + a[1]*C1 + a[2]*C2 + a[3]*C3 + a[4]*C4 + a[5]*C5;
  return H0;
} /***************************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double vcs_Cp0_NASA(double TKelvin, VCS_NASA_POLY *poly_ptr)
   
  /**************************************************************************
   *
   * vcs_Cp0_NASA:
   *
   *   Calculate the standard state Heat Capacity at constant pressure
   *   (in Kelvin) for a single species using the Nasa polynomial format.
   *
   * Input
   *    TKelvin  = Temperature in Kelvin.
   *    poly_ptr = Pointer to structure containing the NASA Polynomial
   *               coefficients 
   *
   * Return
   *     Cp0 = Heat Capacity at constant Pressure / R -> dimensionless
   *
   * Error Conditions
   *    VCS_THERMO_OUTOFRANGE:
   *         If the input temperature, is out of range of the polynomials,
   *    an error Flag is set, and the temperature is storred in the 
   *    error structure.
   ***************************************************************************/
{
  int iRegion;
  double *a, Cp0;
  double *Tlim = VCS_DATA_PTR(poly_ptr->Tlimits);
  static double Tsave = -10., C2, C3, C4;
  /*
   *  Find the temperature region
   */
  if (TKelvin <= *Tlim) {
#ifdef DEBUG_MODE
    plogf("vcs_Cp0_NASA error: TKelvin below lowest bounds\n");
#endif
    iRegion = 0;
    if (TKelvin <= 0.0) {
      Cp0 = poly_ptr->Acoeff[0][0];
      return Cp0;
    }
    goto L_FOUNDREGION;
  }
  for (iRegion = 0; iRegion < poly_ptr->NumTempRegions; iRegion++) {
    Tlim++;
    if (TKelvin <= *Tlim) goto L_FOUNDREGION;
  }
#ifdef DEBUG_MODE
  plogf("vcs_Cp0_NASA error: TKelvin above highest bounds\n");
#endif

  iRegion--;
 L_FOUNDREGION:;
  a = poly_ptr->Acoeff[iRegion];
  if (Tsave != TKelvin) {
    Tsave = TKelvin;
    C2 = TKelvin * TKelvin;
    C3 = C2 * TKelvin;
    C4 = C3 * TKelvin;
  }
  Cp0 = a[0] + a[1]*TKelvin + a[2]*C2 + a[3]*C3 + a[4]*C4;
  return Cp0;
} /***************************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double vcs_S0_NASA(double TKelvin, VCS_NASA_POLY *poly_ptr)
   
  /**************************************************************************
   *
   * vcs_S0_NASA:
   *
   *   Calculates the standard state Entropy (in Kelvin) for a single species 
   *   using the Nasa polynomial format.
   *
   * Input
   *    TKelvin  = Temperature in Kelvin.
   *    poly_ptr = Pointer to structure containing the NASA Polynomial
   *               coefficients   
   *
   * Return
   *     S0 = Standard State Entropy / R -> unitless
   *
   * Error Conditions
   *    VCS_THERMO_OUTOFRANGE:
   *         If the input temperature, is out of range of the polynomials,
   *    an error Flag is set, and the temperature is storred in the 
   *    error structure.
   ***************************************************************************/
{
  int iRegion;
  double *a, S0;
  double *Tlim = VCS_DATA_PTR(poly_ptr->Tlimits);
  static double Tsave = -10., C0, C2, C3, C4;
  /*
   *  Find the temperature region
   */
  if (TKelvin <= *Tlim) {
#ifdef DEBUG_MODE
    plogf("vcs_S0_NASA error: TKelvin below lowest bounds\n");
#endif
    iRegion = 0;
    if (TKelvin <= 0.0) {
      S0 = 0.0;
      return S0;
    }
    goto L_FOUNDREGION;
  }
  for (iRegion = 0; iRegion < poly_ptr->NumTempRegions; iRegion++) {
    Tlim++;
    if (TKelvin <= *Tlim) goto L_FOUNDREGION;
  }
#ifdef DEBUG_MODE
  plogf("vcs_S0_NASA error: TKelvin above highest bounds\n");
#endif
  iRegion--;

   
 L_FOUNDREGION:;
  a = poly_ptr->Acoeff[iRegion];
  if (Tsave != TKelvin) {
    Tsave = TKelvin;
    C0 = log(TKelvin);
    C2 = TKelvin * TKelvin;
    C3 = C2 * TKelvin;
    C4 = C3 * TKelvin;
    C2 /= 2.0;
    C3 /= 3.0;
    C4 /= 4.0;  
  }
  S0 = a[0]*C0 + a[1]*TKelvin + a[2]*C2 + a[3]*C3 + a[4]*C4 + a[7];
  return S0;
} /***************************************************************************/

}

