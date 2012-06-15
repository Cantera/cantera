/*
 * $Id: Electrode_PhaseStability.cpp,v 1.5 2011/01/25 03:21:33 hkmoffa Exp $
 */

//#ifndef OLDWAY
//#define OLDWAY
//#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <iomanip>

//#include "tok_input_util.h"
//#include "cantera/kernel/mdp_allo.h"

//#include "Electrode_MultiPlateau_NoDiff.h"
//#include "cantera/kernel/NonlinearSolver.h"
//#include "cantera/kernel/FixedChemPotSSTP.h"

#include "vcs_phaseStabilitySolve.h"

#include "vcs_solve.h"
#include "vcs_internal.h"
#include "vcs_VolPhase.h"
#include "vcs_species_thermo.h"
#include "vcs_prob.h"

#include "mdp_allo.h"
#include "clockWC.h"


using namespace Cantera;
using namespace std;


#ifndef SAFE_DELETE
#define SAFE_DELETE(x)  if (x) { delete x;  x = 0;}
#endif

namespace VCSnonideal {

  //====================================================================================================================
  vcs_PhaseStabilitySolve::vcs_PhaseStabilitySolve() :
    neq_(0),
    fValue_(0.0),
    nPhasesToPop_(0),
    phasePopNames_(0),
    jac_(0),
    pSolve_(0),
    printLvl_(0),
    enableExtraPrinting_(0),
    detailedResidPrintFlag_(0)
  {
    // m_resid = new calcPhaseStabFunc_ResidJacEval(this);
  }
 
  //====================================================================================================================
  vcs_PhaseStabilitySolve::vcs_PhaseStabilitySolve(const vcs_PhaseStabilitySolve &right) :
    neq_(0),
    fValue_(0.0),
    nPhasesToPop_(0), 
    phasePopNames_(0),
    jac_(0),
    pSolve_(0),
    printLvl_(0),
    enableExtraPrinting_(0),
    detailedResidPrintFlag_(0)
  { 
    /*
     * Call the assignment operator.
     */
    *this = operator=(right);
  } 
  //====================================================================================================================
  // Destructor
  vcs_PhaseStabilitySolve::~vcs_PhaseStabilitySolve()
  {
    SAFE_DELETE(jac_);
    SAFE_DELETE(pSolve_);
  }
  //======================================================================================================================
  // Assignment operator
  /*
   *  @param right object to be copied
   */
  vcs_PhaseStabilitySolve & vcs_PhaseStabilitySolve::operator=(const vcs_PhaseStabilitySolve &right)
  { 
    /*
     * Check for self assignment.
     */
    if (this == &right) return *this;

 
    neq_ = right.neq_;
    fValue_ = right.fValue_;
    nPhasesToPop_ = right.nPhasesToPop_;
    phasePopIndexList_ = right.phasePopIndexList_;

    phasePopNames_ = right.phasePopNames_;


    phaseMFBig_ = right.phaseMFBig_;
    residVector_ = right.residVector_;
    solnVector_ = right.solnVector_;
 
    atol_= right.atol_;
    ylow_ = right.ylow_;
    yhigh_ = right.yhigh_;
    yval_ = right.yval_;

    SAFE_DELETE(jac_);
    jac_ = new SquareMatrix(*right.jac_);
    SAFE_DELETE(pSolve_);
    pSolve_ = new NonlinearSolver(*right.pSolve_);


    printLvl_  = right.printLvl_;
    enableExtraPrinting_ = right.enableExtraPrinting_;
    detailedResidPrintFlag_ = right.detailedResidPrintFlag_;
   
    /*
     * Return the reference to the current object
     */
    return *this; 
  }

  //====================================================================================================================
  void vcs_PhaseStabilitySolve::setup(const std::vector<int> & phasePopIndexList) 
  {
    /*
     *  The number of phases to pop is equal to the length of the input vector
     */
    nPhasesToPop_ = phasePopIndexList.size();


  }
  //====================================================================================================================

  int vcs_PhaseStabilitySolve::determinePhaseStability(double &retnFunc) 
  {
    printLvl_ = 9;
    /*
     *  In order to get correct rates of progress for phases, it's necessary to put some mass into them
     *  in order to get them to potentially to go backwards. We then have to fix this later.
     */
    for (int iii = 0; iii <  nPhasesToPop_; iii++) {
      int iph = phasePopIndexList_[iii];
      //      m_tPhaseMoles_old_[iii] = fVS_->m_tPhaseMoles_old[iph];
      if (fVS_->m_tPhaseMoles_old[iph] <= 0.0) {
	fVS_->m_tPhaseMoles_old[iph] = 1.0;
      }
      std::vector<double> &mfVector = moleFraction_pl_[iii];
      /*
       *  The default value didn't converge. However, these initial conditions converged robustly. 
       *  Therefore, we have a problem with this method. It seems to be very dependent on the initial conditions
       *  of the solution. I tried a bunch of methods but nothing made it robust. These initial conditions
       *  came from the equilibrium solve of the system. 
       */
      if (iii == 0) {
	mfVector[0] = 0.34;
	mfVector[1] = 0.66;
      }
      if (iii == 1) {
	mfVector[1] = 6.55E-2;
	mfVector[0] = 1.0 - mfVector[1];
      }
    }


    packNonlinSolnVector(DATA_PTR(yval_));
    double time_curr = 0.0;
    int num_newt_its = 0;
    int num_linear_solves = 0;
    int numBacktracks = 0;
    double *ydot = 0;
    pSolve_->setRtol(1.0E-8);
    int nonlinearFlag = pSolve_->solve_nonlinear_problem(NSOLN_TYPE_STEADY_STATE, &yval_[0], ydot, 0.0,
							 time_curr, *jac_,  num_newt_its, num_linear_solves, 
							 numBacktracks, printLvl_);
    if (nonlinearFlag < 0) {
      printf(" vcs_PhaseStabilitySolve::determinePhaseStability():  Unsuccessful Nonlinear Solve\n");
      exit(-1);
    }

    for (int iii = 0; iii <  nPhasesToPop_; iii++) {
      int iph = phasePopIndexList_[iii];  
      std::vector<double> &mfVector = moleFraction_pl_[iii];
      fVS_->m_tPhaseMoles_old[iph] = mfVector[iii];
    }
 

    retnFunc = fValue_;

    if (fValue_ < 1.0) {
      return 1;
    }
    return 0;
  }




  //====================================================================================================================
  std::vector<doublereal> & vcs_PhaseStabilitySolve::moleFractions(int iphase) 
  {
    for (int iii = 0; iii <  nPhasesToPop_; iii++) {
      int ii = phasePopIndexList_[iii];
      if (ii == iphase) {
	return moleFraction_pl_[ii];
      }
    }
    throw CanteraError("", "");
  }
  //====================================================================================================================
  int vcs_PhaseStabilitySolve::nResidEquations() const
  {
    // We have already determined the number of equations
    return neq_;
  }

  //====================================================================================================================
  // Unpack the soln vector
  /*
   *  This function unpacks the solution vector into  m_tPhaseMoles_old,  spMoles_final_, and spMf_final_[]
   */
  int vcs_PhaseStabilitySolve::unpackNonlinSolnVector(const double * const y) {
    //   int index = 0;
    for (int iii = 0; iii < nPhasesToPop_; iii++) {


    }

   
    return 0;
  }

  //====================================================================================================================
  // Pack the soln vector
  /*
   *  This function packs the solution vector
   */
  int vcs_PhaseStabilitySolve::packNonlinSolnVector(double * const y) {
    //int index = 0;
    for (int iii = 0; iii < nPhasesToPop_; iii++) {
 
    }

    return 0;
  } 
  //====================================================================================================================
  // Seed the soln vector 
  void vcs_PhaseStabilitySolve::seedMfVector() {

  }

  //====================================================================================================================
  void vcs_PhaseStabilitySolve::updatePhaseMoleFractions(int iii, int iph) {


  }
  //====================================================================================================================
  /*
   *  Calculate CDot, DDot, DDotLin
   */
  void vcs_PhaseStabilitySolve::extractInfo()
  {
    //   int iph, jph, kph;



  }
 //====================================================================================================================
  double vcs_PhaseStabilitySolve::vcs_phaseStabilitySubSolve(const int iph) {

    /*
     * We will use the _new state calc here
     */
    int kspec, irxn, k, i, kc, kc_spec;
    vcs_VolPhase *Vphase = fVS_->m_VolPhaseList[iph];
    doublereal deltaGRxn;

    // We will do a full newton calculation later, but for now, ...
    bool doSuccessiveSubstitution = true;
    double funcPhaseStability;
    vector<doublereal> X_est(Vphase->nSpecies(), 0.0);
    vector<doublereal> delFrac(Vphase->nSpecies(), 0.0);
    vector<doublereal> E_phi(Vphase->nSpecies(), 0.0);
    vector<doublereal> fracDelta_new(Vphase->nSpecies(), 0.0);
    vector<doublereal> fracDelta_old(Vphase->nSpecies(), 0.0);
    vector<doublereal> fracDelta_raw(Vphase->nSpecies(), 0.0);
    vector<int>        creationGlobalRxnNumbers(Vphase->nSpecies(), -1);
    vcs_dcopy(VCS_DATA_PTR(fVS_->m_deltaGRxn_Deficient), VCS_DATA_PTR(fVS_->m_deltaGRxn_old), fVS_->m_numRxnRdc);

    std::vector<doublereal> m_feSpecies_Deficient(fVS_->m_numComponents, 0.0);
    doublereal damp = 1.0;
    doublereal dampOld = 1.0;
    doublereal normUpdate = 1.0;
    doublereal normUpdateOld = 1.0;
    doublereal sum = 0.0;
    doublereal dirProd = 0.0;
    doublereal dirProdOld = 0.0;

    // get the activity coefficients
    Vphase->sendToVCS_ActCoeff(VCS_STATECALC_OLD, VCS_DATA_PTR(fVS_->m_actCoeffSpecies_new));

    // Get the storred estimate for the composition of the phase if 
    // it gets created
    fracDelta_new = Vphase->creationMoleNumbers(creationGlobalRxnNumbers);


    bool oneIsComponent = false;
    std::vector<int> componentList;

    for (k = 0; k < Vphase->nSpecies(); k++) {
      kspec = Vphase->spGlobalIndexVCS(k);
      if (kspec < fVS_->m_numComponents) {
        oneIsComponent = true;
        componentList.push_back(k);
      }
    }

    for (k = 0; k < fVS_->m_numComponents; k++) {
      m_feSpecies_Deficient[k]  = fVS_->m_feSpecies_old[k];
    }
    normUpdate = 0.1 * vcs_l2norm(fracDelta_new);
    damp = 1.0E-2;

    if (doSuccessiveSubstitution) {

#ifdef DEBUG_MODE
      int KP = 0;
      if (fVS_->m_debug_print_lvl >= 2) {
        plogf("   --- vcs_phaseStabilityTest() called\n");
        plogf("   ---  Its   X_old[%2d]  FracDel_old[%2d]  deltaF[%2d] FracDel_new[%2d]"
              "  normUpdate     damp     FuncPhaseStability\n", KP, KP, KP, KP);
        plogf("   --------------------------------------------------------------"
              "--------------------------------------------------------\n");
      } else if (fVS_->m_debug_print_lvl == 1) {
        plogf("   --- vcs_phaseStabilityTest() called for phase %d\n", iph);
      }
#endif

      for (k = 0; k < Vphase->nSpecies(); k++) {
        if (fracDelta_new[k] < 1.0E-13) {
          fracDelta_new[k] = 1.0E-13;
        }
      }
      bool converged = false;
      for (int its = 0; its < 200  && (!converged); its++) {

        dampOld = damp;
        normUpdateOld = normUpdate;
        fracDelta_old = fracDelta_new;
        dirProdOld = dirProd;



        // Given a set of fracDelta's, we calculate the fracDelta's
        // for the component species, if any
        for (i = 0; i < (int) componentList.size(); i++) {
          kc = componentList[i];
          kc_spec = Vphase->spGlobalIndexVCS(kc);
          fracDelta_old[kc] = 0.0;
          for (k = 0; k <  Vphase->nSpecies(); k++) {
            kspec = Vphase->spGlobalIndexVCS(k);
            irxn = kspec - fVS_->m_numComponents;
            if (irxn >= 0) {
              fracDelta_old[kc] += fVS_->m_stoichCoeffRxnMatrix[irxn][kc_spec] *  fracDelta_old[k];
            }
          }
        }
        // Now, calculate the predicted mole fractions, X_est[k]
        double sumFrac = 0.0;
        for (k = 0; k < Vphase->nSpecies(); k++) {
          sumFrac += fracDelta_old[k];
        }
        // Necessary because this can be identically zero. -> we need to fix this algorithm!
        if (sumFrac <= 0.0) {
          sumFrac = 1.0;
        }
        double sum_Xcomp = 0.0;
        for (k = 0; k < Vphase->nSpecies(); k++) {
          X_est[k] = fracDelta_old[k] / sumFrac;
          kc_spec = Vphase->spGlobalIndexVCS(k);
          if (kc_spec < fVS_->m_numComponents) {
            sum_Xcomp += X_est[k];
          }
        }


        /*
         * Feed the newly formed estimate of the mole fractions back into the
         * ThermoPhase object
         */
        Vphase->setMoleFractionsState(0.0, VCS_DATA_PTR(X_est), VCS_STATECALC_PHASESTABILITY);

        /*
         *   get the activity coefficients
         */
        Vphase->sendToVCS_ActCoeff(VCS_STATECALC_OLD, VCS_DATA_PTR(fVS_->m_actCoeffSpecies_new));

        /*
         * First calculate altered chemical potentials for component species
         * belonging to this phase.
         */
        for (i = 0; i < (int) componentList.size(); i++) {
          kc = componentList[i];
          kc_spec = Vphase->spGlobalIndexVCS(kc);
          if ( X_est[kc] > VCS_DELETE_MINORSPECIES_CUTOFF) {
            m_feSpecies_Deficient[kc_spec] = fVS_->m_feSpecies_old[kc_spec] + log(fVS_->m_actCoeffSpecies_new[kc_spec] * X_est[kc]);
          } else {
            m_feSpecies_Deficient[kc_spec] = fVS_->m_feSpecies_old[kc_spec] + log(fVS_->m_actCoeffSpecies_new[kc_spec] * VCS_DELETE_MINORSPECIES_CUTOFF);
          }
        }

        for (i = 0; i < (int) componentList.size(); i++) {
          kc = componentList[i];
          kc_spec = Vphase->spGlobalIndexVCS(kc);

          for (k = 0; k <  Vphase->nSpecies(); k++) {
            kspec = Vphase->spGlobalIndexVCS(k);
            irxn = kspec - fVS_->m_numComponents;
            if (irxn >= 0) {
              if (i == 0) {
                fVS_->m_deltaGRxn_Deficient[irxn] = fVS_->m_deltaGRxn_old[irxn];
              }
              double *dtmp_ptr = fVS_->m_stoichCoeffRxnMatrix[irxn];
              if (dtmp_ptr[kc_spec] != 0.0) {
                fVS_->m_deltaGRxn_Deficient[irxn] += dtmp_ptr[kc_spec] * (m_feSpecies_Deficient[kc_spec]- fVS_->m_feSpecies_old[kc_spec]);
              }
            }

          }
        }
        /*
         *  Calculate the E_phi's
         */
        sum = 0.0;
        funcPhaseStability = sum_Xcomp - 1.0;
        for (k = 0; k <  Vphase->nSpecies(); k++) {
          kspec = Vphase->spGlobalIndexVCS(k);
          irxn = kspec - fVS_->m_numComponents;
          if (irxn >= 0) {
            deltaGRxn = fVS_->m_deltaGRxn_Deficient[irxn];
            if (deltaGRxn >  50.0) deltaGRxn =  50.0;
            if (deltaGRxn < -50.0) deltaGRxn = -50.0;
            E_phi[k] = std::exp(-deltaGRxn) / fVS_->m_actCoeffSpecies_new[kspec];
            sum +=  E_phi[k];
            funcPhaseStability += E_phi[k];
          } else {
            E_phi[k] = 0.0;
          }
        }

        /*
         * Calculate the raw estimate of the new fracs
         */
        for (k = 0; k <  Vphase->nSpecies(); k++) {
          kspec = Vphase->spGlobalIndexVCS(k);
          irxn = kspec - fVS_->m_numComponents;
          double b =  E_phi[k] / sum * (1.0 - sum_Xcomp);
          if (irxn >= 0) {
            fracDelta_raw[k] = b;
          }
        }


        // Given a set of fracDelta's, we calculate the fracDelta's
        // for the component species, if any
        for (i = 0; i < (int) componentList.size(); i++) {
          kc = componentList[i];
          kc_spec = Vphase->spGlobalIndexVCS(kc);
          fracDelta_raw[kc] = 0.0;
          for (k = 0; k <  Vphase->nSpecies(); k++) {
            kspec = Vphase->spGlobalIndexVCS(k);
            irxn = kspec - fVS_->m_numComponents;
            if (irxn >= 0) {
              fracDelta_raw[kc] += fVS_->m_stoichCoeffRxnMatrix[irxn][kc_spec] * fracDelta_raw[k];
            }
          }
        }

        /*
         * Now possibly dampen the estimate.
         */
        doublereal sumADel = 0.0;
        for (k = 0; k <  Vphase->nSpecies(); k++) {
          delFrac[k] = fracDelta_raw[k] - fracDelta_old[k];
          sumADel += fabs(delFrac[k]);
        }
        normUpdate = vcs_l2norm(delFrac);

        dirProd = 0.0;
        for (k = 0; k <  Vphase->nSpecies(); k++) {
          dirProd += fracDelta_old[k] * delFrac[k];
        }
        bool crossedSign = false;
        if (dirProd * dirProdOld < 0.0) {
          crossedSign = true;
        }


        damp = 0.5;
        if (dampOld < 0.25) {
          damp = 2.0 * dampOld;
        }
        if (crossedSign) {
          if (normUpdate *1.5 > normUpdateOld) {
            damp = 0.5 * dampOld;
          } else if (normUpdate *2.0 > normUpdateOld) {
            damp = 0.8 * dampOld;
          }
        } else {
          if (normUpdate > normUpdateOld * 2.0) {
            damp = 0.6 * dampOld;
          } else if (normUpdate > normUpdateOld * 1.2) {
            damp = 0.9 * dampOld;
          }
        }

        for (k = 0; k < Vphase->nSpecies(); k++) {
          if (fabs(damp * delFrac[k]) > 0.3*fabs(fracDelta_old[k])) {
            damp = MAX(0.3*fabs(fracDelta_old[k]) / fabs( delFrac[k]),
                       1.0E-8/fabs( delFrac[k]));
          }
          if (delFrac[k] < 0.0) {
            if (2.0 * damp * (-delFrac[k]) > fracDelta_old[k]) {
              damp = fracDelta_old[k] / (2.0 * (-delFrac[k]));
            }
          }
          if (delFrac[k] > 0.0) {
            if (2.0 * damp * delFrac[k] > fracDelta_old[k]) {
              damp = fracDelta_old[k] / (2.0 * delFrac[k]);
            }
          }
        }
        if (damp < 0.000001) {
          damp = 0.000001;
        }

        for (k = 0; k <  Vphase->nSpecies(); k++) {
          fracDelta_new[k] = fracDelta_old[k] + damp * (delFrac[k]);
        }

#ifdef DEBUG_MODE
        if (fVS_->m_debug_print_lvl >= 2) {
          plogf("  --- %3d %12g %12g %12g %12g %12g %12g %12g\n", its, X_est[KP], fracDelta_old[KP],
                delFrac[KP], fracDelta_new[KP], normUpdate, damp, funcPhaseStability);
        }
#endif

        if (normUpdate < 1.0E-5) {
          converged = true;
        }

      }

      if (converged) {
        Vphase->setMoleFractionsState(0.0, VCS_DATA_PTR(X_est),
                                      VCS_STATECALC_PHASESTABILITY);
        Vphase->setCreationMoleNumbers(VCS_DATA_PTR(fracDelta_new), creationGlobalRxnNumbers);
      }


    } else {
      printf("not done yet\n");
      exit(-1);
    }
#ifdef DEBUG_MODE
    if (fVS_->m_debug_print_lvl >= 2) {
      plogf("  ------------------------------------------------------------"
            "-------------------------------------------------------------\n");
    } else if (fVS_->m_debug_print_lvl == 1) {
      if (funcPhaseStability > 0.0) {
        plogf("  --- phase %d with func = %g is to be born\n", iph, funcPhaseStability);
      } else {
        plogf("  --- phase %d with func = %g stays dead\n", iph, funcPhaseStability);
      }
    }
#endif
    return funcPhaseStability;
  }
  //====================================================================================================================





  //====================================================================================================================
  // Evaluate the residual function
  /*
   * @param t             Time                    (input) 
   * @param delta_t       The current value of the time step (input)
   * @param y             Solution vector (input, do not modify)
   * @param ydot          Rate of change of solution vector. (input, do not modify)
   * @param resid         Value of the residual that is computed (output)
   * @param evalType      Type of the residual being computed (defaults to Base_ResidEval)
   * @param id_x          Index of the variable that is being numerically differenced to find
   *                      the jacobian (defaults to -1, which indicates that no variable is being
   *                      differenced or that the residual doesn't take this issue into account)
   * @param delta_x       Value of the delta used in the numerical differencing
   */
  int vcs_PhaseStabilitySolve::optResid(const doublereal tdummy, const doublereal delta_t_dummy,
					 const doublereal * const y,
					 const doublereal * const ySolnDot,
					 doublereal * const resid,
					 const ResidEval_Type_Enum evalType,
					 const int id_x, 
					 const doublereal delta_x)
  {

  
    if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
      printf("\t\t================================================================================================="
	     "==============================\n");
      printf("\t\t  EXTRA PRINTING FROM NONLINEAR RESIDUAL: ");
      if (evalType == Base_ResidEval || evalType == Base_LaggedSolutionComponents) {
	printf(" BASE RESIDUAL");
      } else if (evalType == JacBase_ResidEval) {
	printf(" BASE JAC RESIDUAL");
      } else  if  (evalType == JacDelta_ResidEval) {
	printf(" DELTA JAC RESIDUAL");
	printf(" var = %d delta_x = %12.4e Y_del = %12.4e Y_base = %12.4e", id_x, delta_x, y[id_x], y[id_x] - delta_x);
      } else  if  (evalType == Base_ShowSolution) {
	printf(" BASE RESIDUAL - SHOW SOLUTION");
      }
      printf("\n");
    } 
    /*
     *  Current the solution vector are the nsp -1 mole fractions in multiple phases
     */
    unpackNonlinSolnVector(y);


    extractInfo();


    if (enableExtraPrinting_ && detailedResidPrintFlag_ > 1) {
      printf("\t\t=============================================================="
	     "=================================================================\n");
    }

    return 1;
  }


//====================================================================================================================
// Evaluate the residual function
/*
 * @param t             Time                    (input) 
 * @param delta_t       The current value of the time step (input)
 * @param y             Solution vector (input, do not modify)
 * @param ydot          Rate of change of solution vector. (input, do not modify)
 * @param resid         Value of the residual that is computed (output)
 * @param evalType      Type of the residual being computed (defaults to Base_ResidEval)
 * @param id_x          Index of the variable that is being numerically differenced to find
 *                      the jacobian (defaults to -1, which indicates that no variable is being
 *                      differenced or that the residual doesn't take this issue into account)
 * @param delta_x       Value of the delta used in the numerical differencing
 */
  int vcs_PhaseStabilitySolve::evalResidNJ(const doublereal tdummy, const doublereal delta_t_dummy,
	    const doublereal * const y,
	    const doublereal * const ySolnDot,
	    doublereal * const resid,
	    const ResidEval_Type_Enum evalType,
	    const int id_x, 
	    const doublereal delta_x)
{
    /*
     *  Return control for calculating the residual to the controlling vcs_PhaseStabilitySolve object.
     *  This object will calculate the residual at the current conditions
     */
  int retn =  optResid(tdummy, delta_t_dummy, y, ySolnDot, resid, evalType, id_x, delta_x);
  return retn;
}
//  -----------------------------------------------------------------------------------------------------------------
  int vcs_PhaseStabilitySolve::getInitialConditions(const doublereal t0, doublereal * const y, doublereal * const ydot)
{
  for (int k = 0; k < neq_; k++) {
    y[k] = 0.0;
  }
  return 1;
}
//  -----------------------------------------------------------------------------------------------------------------
int vcs_PhaseStabilitySolve::nEquations() const
{
  return neq_;
}
//====================================================================================================================
  //  Apply a filtering process to the step
  /*
   *  @param timeCurrent    Current value of the time
   *  @param ybase          current value of the solution
   *  @param step0          Value of the step in the solution vector that will be filtered.
   *                        The filter is applied to the step values.
   *
   *  @return Returns the norm of the value of the amount filtered
   */
doublereal vcs_PhaseStabilitySolve::filterNewStep(const doublereal timeCurrent, const doublereal * const ybase,
						  doublereal * const step0)
{
  return 0.0;
}
//====================================================================================================================
  //   Determine the big mole fraction in the phase
  void  vcs_PhaseStabilitySolve::determineBigMoleFractions() {

  }
  //====================================================================================================================
} // End of namespace Cantera
//======================================================================================================================

