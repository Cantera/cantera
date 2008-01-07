/*
 * $Id$
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef VCS_NASA_POLY_H
#define VCS_NASA_POLY_H

#include <string>
#include <vector>
#include "vcs_DoubleStarStar.h"

namespace VCSnonideal {

/*
*                  NASA Polynomial Form for Standard state Thermo Functions.
*
*
*  NumberTempRegions
*              Number of temperature regions in the fits: 
*              Must be greater or equal to one.
*
*  Tlimits[NumberTempRegions+1]:
*              Temperature limits of the regions. At the intersection of
*              the regions, the polynomial formulas are suppose to be
*              C1 continuous.
*              To Locate Region i for current temperature, TKelvin:
*                   Tlimits[i] <= TKelvin < Tlimits[i+1]
*
*  Acoeff[NumberTempRegions][7]
*        Coefficients for calculation of the standard state thermodynamic
*        functions.
*
*        double *a;
*        for i such that Tlimits[i] <= T < Tlimits[i+1]:
*        a = Acoeff[i];
*
*        C_p/R = a[0] + a[1]*T + a[2] * T^2 + a[3] * T^3 + a[4] * T^4 
*
*        H/RT  = a[0] + a[1]/2*T + a[2]/3 * T^2 + a[3]/4 * T^3 + a[4]/5 * T^4 
*                     + a[5]/T
*
*        S/R   = a[0] * log(T) + a[1] * T + a[2]/2 * T^2 + a[3]/3 * T^3
*                              + a[4]/4 * T^4  + a[6]
*
*/
class VCS_NASA_POLY {
public:
  VCS_NASA_POLY(int, int);
  VCS_NASA_POLY(const VCS_NASA_POLY &b);
  VCS_NASA_POLY& operator=(const VCS_NASA_POLY &);

  ~VCS_NASA_POLY();
  int NumTempRegions;
/* Vector Of Temperature Limits -> One More Than
		       The Number Of Regions */
  std::vector<double> Tlimits;

  DoubleStarStar Acoeff;

  int     NumEl;
  std::vector<double> ElComp;
  std::vector<std::string> ElName;
  char Date[12];
  char PhType;
  char SpName[24];
  char PhName[24];
};

/* Externals for vcs_nasa_poly.c */

extern VCS_NASA_POLY *vcs_nasa_poly_create(int, int);
extern void vcs_nasa_poly_free(VCS_NASA_POLY *);
extern void vcs_nasa_poly_destroy(VCS_NASA_POLY **);

extern double vcs_G0_NASA(double, VCS_NASA_POLY *);
extern double vcs_H0_NASA(double, VCS_NASA_POLY *);
extern double vcs_Cp0_NASA(double, VCS_NASA_POLY *);
extern double vcs_S0_NASA(double, VCS_NASA_POLY *);

}

#endif
