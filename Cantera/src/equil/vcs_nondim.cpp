/**
 *  @file vcs_nondim.cpp
 *     Nondimensionalization routines with VCSnonideal
 */
/*
 * $Id$
 */
/*
 * Copywrite (2007) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vcs_solve.h"
#include "vcs_internal.h" 

namespace VCSnonideal {

/**************************************************************************
 *
 * vcs_nondimMult:
 *
 *   Returns the multiplier for the nondimensionalization of the equations
 *    (this is basically equal to RT)
 **************************************************************************/
double VCS_SOLVE::vcs_nondim_Farad(int mu_units, double TKelvin)
{
  double Farad;
   if (TKelvin <= 0.0) TKelvin = 293.15;
   switch (mu_units) {
   case VCS_UNITS_MKS:
   case VCS_UNITS_KJMOL: 
   case VCS_UNITS_KCALMOL:  
     Farad = 1.602E-19 * 6.022136736e26/ (TKelvin * 8.314472E3);
     break;
   case VCS_UNITS_UNITLESS: 
     Farad = 1.602E-19 * 6.022136736e26;
     break;   
   case VCS_UNITS_KELVIN: 
     Farad = 1.602E-19 * 6.022136736e26/ (TKelvin);  
     break;
   default:
       plogf("vcs_nondim_Farad error: unknown units: %d\n", mu_units);
       exit(-1);
   }
   return Farad;
}

double VCS_SOLVE::vcs_nondimMult_TP(int mu_units, double TKelvin)
{
   double rt;
   if (TKelvin <= 0.0) TKelvin = 293.15;
   switch (mu_units) {
   case VCS_UNITS_KCALMOL:  
       rt = TKelvin * 8.314472E-3 / 4.184;
       break;
   case VCS_UNITS_UNITLESS: 
       rt = 1.0;
       break;
   case VCS_UNITS_KJMOL: 
       rt = TKelvin * 0.008314472;
       break;
   case VCS_UNITS_KELVIN: 
       rt = TKelvin;
       break;
   case VCS_UNITS_MKS:
       rt = TKelvin * 8.314472E3;
       break;
   default:
       plogf("vcs_nondimMult_TP error: unknown units: %d\n", mu_units);
       exit(-1);
   }
   return rt;
}

/**************************************************************************
 *
 *  vcs_nondim_TP:
 *        Nondimensionalize the problem data:
 *          ->nondimensionalize the free energies using
 *            the divisor, R * T
 *
 *
 * HKM -> I don't think we need to modify the mole nubmers by 1E3 for the
 * case of MKS units. However, what we need to do is to add a scale
 * factor so that the number of moles or kmoles is ~ 1.0. Many of the
 * algorithms rely on this I think in a subtle way. This is the perfect
 * place to add this in.
 **************************************************************************/
void VCS_SOLVE::vcs_nondim_TP(void) {
    int i;
    double tf;
    if (UnitsState == VCS_DIMENSIONAL_G) {
      UnitsState = VCS_NONDIMENSIONAL_G;
      tf = 1.0 / vcs_nondimMult_TP(m_VCS_UnitsFormat, T);
      for (i = 0; i < m_numSpeciesTot; ++i) {
	/* 
	 *        Modify the standard state and total chemical potential data,
	 *        FF(I),  to make it dimensionless, i.e.,  mu / RT.
	 *        Thus, we may divide it by the temperature.
	 */  
	m_SSfeSpecies[i] *= tf;
	m_gibbsSpecies[i] *= tf;
	dg[i] *= tf;
	dgl[i] *= tf;
	fel[i] *= tf;
      }

      Faraday_dim =  vcs_nondim_Farad(m_VCS_UnitsFormat, T);
      if (m_VCS_UnitsFormat == VCS_UNITS_MKS) {
	for (i = 0; i < m_numSpeciesTot; ++i) {
	  if (SpeciesUnknownType[i] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	    soln[i] *= 1.0E3;
	  }
	}
	for (i = 0; i < m_numElemConstraints; ++i) {
	  gai[i] *= 1.0E3;
	}
      }
    }
} /* vcs_nondim_TP() *********************************************************/

/**************************************************************************
 *
 *  vcs_nondim_TP:
 *        Redimensionalize the problem data:
 *          ->redimensionalize the free energies using the reverse
 *            of vcs_nondim_TP
 **************************************************************************/
void VCS_SOLVE::vcs_redim_TP(void)
{
    int i;
    double tf;
    if (UnitsState != VCS_DIMENSIONAL_G) {
      UnitsState = VCS_DIMENSIONAL_G;
      tf = vcs_nondimMult_TP(m_VCS_UnitsFormat, T);
      for (i = 0; i < m_numSpeciesTot; ++i) {
	/* 
	 *        Modify the standard state and total chemical potential data,
	 *        FF(I),  to make it have units, i.e. mu = RT * mu_star
	 */  
	m_SSfeSpecies[i] *= tf;
	m_gibbsSpecies[i] *= tf;
	dg[i] *= tf;
	dgl[i] *= tf;
	fel[i] *= tf;
      } 
      Faraday_dim *= tf;
    }
    if (m_VCS_UnitsFormat == VCS_UNITS_MKS) {
	for (i = 0; i < m_numSpeciesTot; ++i) {
	  if (SpeciesUnknownType[i] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	    soln[i] /= 1.0E3;
	  }
	}
	for (i = 0; i < m_numElemConstraints; ++i) {
	  gai[i] /= 1.0E3;
	}
    }

} /* vcs_redim_TP() **********************************************************/

void VCS_SOLVE::vcs_printChemPotUnits(int unitsFormat) {
    switch(unitsFormat) {
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
	exit(-1);
    }
}

}

