/**
 *  @file vcs_nondim.cpp
 *     Nondimensionalization routines within VCSnonideal
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

  //  Returns the multiplier for electric charge terms
  /*
   *   This is basically equal to F/RT
   *
   * @param mu_units integer representing the dimensional units system
   * @param TKelvin  double  Temperature in Kelvin
   *
   * @return Returns the value of F/RT
   */
  double VCS_SOLVE::vcs_nondim_Farad(int mu_units, double TKelvin) const {
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
      plogendl();
      exit(-1);
    }
    return Farad;
  }

  //  Returns the multiplier for the nondimensionalization of the equations
  /*
   *   This is basically equal to RT
   *
   * @param mu_units integer representing the dimensional units system
   * @param TKelvin  double  Temperature in Kelvin
   *
   * @return Returns the value of RT
   */
  double VCS_SOLVE::vcs_nondimMult_TP(int mu_units, double TKelvin) const {
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
      plogendl();
      exit(-1);
    }
    return rt;
  }

  // Nondimensionalize the problem data
  /*
   *   Nondimensionalize the free energies using the divisor, R * T
   *
   *  Essentially the internal data can either be in dimensional form
   *  or in nondimensional form. This routine switches the data from 
   *  dimensional form into nondimensional form.
   *
   *  What we do is to divide by RT.
   *
   *  @todo Add a scale factor based on the total mole numbers.
   *        The algorithm contains hard coded numbers based on the
   *        total mole number. If we ever were faced with a problem
   *        with significantly different total kmol numbers than one
   *        the algorithm would have problems.
   */
  void VCS_SOLVE::vcs_nondim_TP() {
    int i;
    double tf;
    if (m_unitsState == VCS_DIMENSIONAL_G) {
      m_unitsState = VCS_NONDIMENSIONAL_G;
      tf = 1.0 / vcs_nondimMult_TP(m_VCS_UnitsFormat, m_temperature);
      for (i = 0; i < m_numSpeciesTot; ++i) {
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
      if (m_VCS_UnitsFormat == VCS_UNITS_MKS) {
	for (i = 0; i < m_numSpeciesTot; ++i) {
	  if (m_speciesUnknownType[i] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	    //m_molNumSpecies_old[i] *= 1.0E3;
	    m_molNumSpecies_old[i] *= 1.0;
	  }
	}
	for (i = 0; i < m_numElemConstraints; ++i) {
	  //m_elemAbundancesGoal[i] *= 1.0E3;
	  m_elemAbundancesGoal[i] *= 1.0;
	}
      }
    }
  }

  // Redimensionalize the problem data
  /*
   *  Redimensionalize the free energies using the multiplier R * T
   *
   *  Essentially the internal data can either be in dimensional form
   *  or in nondimensional form. This routine switches the data from 
   *  nondimensional form into dimensional form.
   *
   *  What we do is to multiply by RT.
   */
  void VCS_SOLVE::vcs_redim_TP(void)
  {
    int i;
    double tf;
    if (m_unitsState != VCS_DIMENSIONAL_G) {
      m_unitsState = VCS_DIMENSIONAL_G;
      tf = vcs_nondimMult_TP(m_VCS_UnitsFormat, m_temperature);
      for (i = 0; i < m_numSpeciesTot; ++i) {
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
    if (m_VCS_UnitsFormat == VCS_UNITS_MKS) {
      for (i = 0; i < m_numSpeciesTot; ++i) {
	if (m_speciesUnknownType[i] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	  //m_molNumSpecies_old[i] /= 1.0E3;
	  m_molNumSpecies_old[i] /= 1.0;
	}
      }
      for (i = 0; i < m_numElemConstraints; ++i) {
	//m_elemAbundancesGoal[i] /= 1.0E3;
	m_elemAbundancesGoal[i] /= 1.0;
      }
    }

  } 

  // Computes the current elemental abundances vector
  /*
   *   Computes the elemental abundances vector, m_elemAbundances[], and stores it
   *   back into the global structure
   */
  void VCS_SOLVE::vcs_printChemPotUnits(int unitsFormat) const {
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

