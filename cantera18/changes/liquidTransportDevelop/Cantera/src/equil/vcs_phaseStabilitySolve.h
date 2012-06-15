/**
 * Predict the stability of a single phase that is currently zereod.
 *
 *  This routine fills in an estimate for the solution
 *  Return 1 if the phases are stable and 0 if they are not
 *
 */
/*
 * $Id: Electrode_PhaseStability.h,v 1.5 2011/01/25 03:21:33 hkmoffa Exp $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _VCS_PHASESTABILITYSOLVE_H
#define _VCS_PHASESTABILITYSOLVE_H

#include <vector>
#include <string>

#include "ct_defs.h"
#include "vcs_defs.h"
#include "vcs_DoubleStarStar.h"
#include "vcs_IntStarStar.h"
#include "vcs_solve.h"

#include "ResidJacEval.h"
#include "NonlinearSolver.h"

namespace VCSnonideal {


  //! Class which Solves the stability problem for a phase
  /*!
   *  The class doubles for the residual equation for the nonlinear solver.
   *
   *  In one part of the class, we set up the nonlinear problem. Then, we hand the calculation
   *  off to the nonlienar solver, which calls back to this class to calculate the residual.
   */
  class vcs_PhaseStabilitySolve : public Cantera::ResidJacEval { 

  public:

    //! Constructor
    /*!
     * 
     * @param elect  Electrode object pertaining to the stability problem
     *               The electrode object assigns this object as a "friend"
     */
    vcs_PhaseStabilitySolve();

    //! Destructor
    virtual ~vcs_PhaseStabilitySolve();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    vcs_PhaseStabilitySolve(const vcs_PhaseStabilitySolve &right);

    //! Assignment operator
    /*!
     *  @param right object to be copied
     */
    vcs_PhaseStabilitySolve & operator=(const vcs_PhaseStabilitySolve &right);

    double vcs_phaseStabilitySubSolve(const int iph);
    //! Setup the stability problem
    void setup(const std::vector<int> & phasePopIndexList);

    int determinePhaseStability(doublereal &retnFunc);

    std::vector<doublereal> & moleFractions(int iphase);

    int  nResidEquations() const;

    //! Unpack the soln vector
    /*!
     *  This function unpacks the solution vector into  phaseMoles_final_,  spMoles_final_, and spMf_final_[]
     */
    int unpackNonlinSolnVector(const double * const y);

    int packNonlinSolnVector(double * const y);

    void  seedMfVector();
    void updatePhaseMoleFractions(int iii, int iph);

    void extractInfo();
   
    //! Residual Evaluation program
    /*!
     *
     */
    int optResid(const doublereal tdummy, const doublereal delta_t_dummy,
		 const doublereal * const y,
		 const doublereal * const ySolnDot,
		 doublereal * const resid,
		 const Cantera::ResidEval_Type_Enum evalType,
		 const int id_x, 
		 const doublereal delta_x);

    void determineBigMoleFractions();

    // ----------------------------------------------------------------------------------------------
    //         Functions which are inherited from the Residual Interface
    // ----------------------------------------------------------------------------------------------

    //! Evaluate the residual function
    /*!
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
     *
     * @return 
     */
    virtual int evalResidNJ(const doublereal t, const doublereal delta_t,
			    const doublereal * const y,
			    const doublereal * const ydot,
			    doublereal * const resid,
			    const Cantera::ResidEval_Type_Enum evalType = Cantera::Base_ResidEval,
			    const int id_x = -1, 
			    const doublereal delta_x = 0.0);

    //! Fill in the initial conditions
    /*!
     * Values for both the solution and the value of ydot may be provided.
     *
     * @param t0            Time                    (input) 
     * @param y             Solution vector (output)
     * @param ydot          Rate of change of solution vector. (output)
     */
    virtual int getInitialConditions(const doublereal t0, doublereal * const y,
				     doublereal * const ydot);
      
      
    //! Return the number of equations in the equation system
    virtual  int nEquations() const;


    //! Apply a filtering process to the step
    /*!
     *  @param timeCurrent    Current value of the time
     *  @param ybase          current value of the solution
     *  @param step0          Value of the step in the solution vector that will be filtered.
     *                        The filter is applied to the step values.
     *
     *  @return Returns the norm of the value of the amount filtered
     */
    virtual doublereal filterNewStep(const doublereal timeCurrent, const doublereal * const ybase,
				     doublereal * const step0);

  protected:

    //! This is a reference to the friend object, where we will pull most of the data from.
    VCS_SOLVE *fVS_;

    //! Number of equations in the nonlinear problem
    int neq_;


    //! Resulting Final value expressing whether the phase set is ready to pop.
    double fValue_;


    //! Number of phases to check
    int nPhasesToPop_;

    //! Vector of indexes of phases to check in the phase pop problem
    std::vector<int> phasePopIndexList_;

  
    std::vector<std::string> phasePopNames_;
    


    std::vector<int> phaseMFBig_;

    //! Residual vector for the nonlinear residual
    /*!
     *  The ordering of the residual vector is always by phase id first. Then
     *  by local species number within the phase
     */
    std::vector<doublereal> residVector_;

    //! This is the solution vector, the vector of independent unknowns
    /*!
     *  The solution vector is the list of reaction rates of progress for the reactions
     *  that pop the phase into existence
     *
     *  The ordering of the residual vector is always by phase id first. Then
     *  by local species number within the phase
     */
    std::vector<doublereal> solnVector_;

    //! Vector of Mole Number estimates for the predicted solution that maximizes
    //! the potential for the phase(s) to be stable
    /*
     *  Outer loop is over the number of phases
     */
    std::vector< std::vector<doublereal> > moleNumber_pl_;

    //! Vector of Mole fraction estimates for the predicted solution that maximizes
    //! the potential for the phase(s) to be stable
    /*
     *  Outer loop is over the number of phases
     */
    std::vector< std::vector<doublereal> > moleFraction_pl_;

    //! Vector of integers which identify the reaction id in the VCS_SOLVE problem
    //! that the solution vector is identified with
    /*!
     *  The mapping will be from this object's solution index to the VCS_SOLVE reaction index
     *  This mapping must be unique.
     *  
     */
    std::vector<int> IrxnIndex_;

    //! Number of reaction indeces associated with each phase
    /*!
     *   Note, this does not have to equal the number of species in the phase.
     *   There can be more Reaction indeces that species, if one of the species 
     *   in the phase is a component, and the component reactions are linearly independent.
     *
     *   This is a vector of length number of phases.
     */
    std::vector<int> numIrxn_pl_;


 
    std::vector<doublereal> atol_;
    std::vector<doublereal> ylow_;
    std::vector<doublereal> yhigh_;
    std::vector<doublereal> yval_;

    Cantera::SquareMatrix *jac_;
    Cantera::NonlinearSolver *pSolve_;

    
    int printLvl_;
  public:
    int enableExtraPrinting_;
    int detailedResidPrintFlag_;
  };

}
#endif
/*****************************************************************************/
                                                                                            
