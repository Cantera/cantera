/* ======================================================================= */
/* -------------------------------------------------- */
/* | RCS Head Information on zuzax.pchem.sandia.gov | */
/* -------------------------------------------------- */
/* $RCSfile$ */
/* $Author$ */
/* $Date$ */
/* $Revision$ */
/* ======================================================================= */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vcs_solve.h"
#include "vcs_internal.h" 
#include "vcs_VolPhase.h"
#include "vcs_species_thermo.h"

namespace VCSnonideal {

/**************************************************************************
 *
 * vcs_TV:
 *
 *      Solve an equilibrium problem at a particular fixed temperature 
 *      and volume.
 *      This is done as a root finder problem, solving repetative
 *      calls to solve_TP.
 *
 *     ipr = 1 -> Print results to standard output 
 *           0 -> don't report on anything 
 *     ip1 = 1 -> Print intermediate results. 
 *     maxit -> Maximum number of iterations for the algorithm 
 *      T = Temperature (Kelvin)
 *     Pres = Pressure (units specififed by if__ variable)
 */
int VCS_SOLVE::vcs_TV(int ipr, int ip1, int maxit, double T_arg, double VolRequest)
{
   int iconv, varID;
   double Pmin, Pmax, Preturn;
   VCS_FUNC_PTR func;

   /*
   *        Store the temperature in the private global variables
   */
   T    = T_arg;
   
   /*
   *         Set the unknown variable to the pressure
   */
   varID = 1;
   
   /*
   *         Set the function to the volume function
   */
   func =  vcs_funcVtot;
   
   /*
   * Set max and min Pressures
   */
   Pmin = 0.0;
   Pmax = 1.0E30;
   Preturn = 1.0;
   
   iconv =  vcsUtil_root1d(Pmin, Pmax, maxit, func, (void *) this,
			   VolRequest, varID, &Preturn);
 
   /*
   *        Return the convergence success flag.
   */
   return iconv;
}

 /**************************************************************************
  *
  * vcs_VolTotal
  *
  *  This function calculates the partial molar volume
  *  for all species, kspec, in the thermo problem
  *  at the temperature TKelvin and pressure, Pres, pres is in atm.
  *  And, it calculates the total volume of the combined system.
  *
  * Input
  *    iphase
  *    TKelvin
  *    pres
  *    w[] => vector containing the current mole numbers.
  *
  * Output
  *    VolPM[] => For species in all phase, the entries are the
  *                partial molar volumes
  *    return value = Total volume of the phase in L**3 / MOL_UNITS
  *
  *  (L and MOL_UNITS determined from global units value if__)
  */
double VCS_SOLVE::vcs_VolTotal(double tkelvin, double pres, double w[],
			       double VolPM[])
{
   double volTot = 0.0;
   for (int iphase = 0; iphase < NPhase; iphase++) {
     vcs_VolPhase *Vphase = VPhaseList[iphase];
     Vphase->setState_TP(tkelvin, pres);
     Vphase->setMolesFromVCS(w);
     double volp = Vphase->VolPM_calc();
     (void) Vphase->sendToVCSVolPM(VolPM);
     volTot += volp;
   }
   return volTot;
} 
/**************************************************************************/
}

