/**
 * @file vcs_solve.h
 *    Header file for the internal object that holds the problem
 */
/*
 * $Id$
 */
/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */


#ifndef _VCS_SOLVE_H
#define _VCS_SOLVE_H

/*     
* Index of Symbols
* -------------------
*     irxn -> refers to the species or rxn between the species and
*             the components in the problem
*     k    -> refers to the species
*     j    -> refers to the element or component 
*
*     ###  -> to be eliminated
*/
#include <vector>
#include <string>

#include "ct_defs.h"
#include "vcs_defs.h"
#include "vcs_DoubleStarStar.h"
#include "vcs_IntStarStar.h"

namespace VCSnonideal {
/*
 * Forward references
 */
class vcs_VolPhase;
class VCS_SPECIES_THERMO;
class VCS_PROB;
class VCS_COUNTERS;


//!    This is the main structure used to hold the internal data
//!    used in vcs_solve_TP(), and to solve TP systems.
/*!
 *      The indecises of information in this
 *      structure may change when the species basis changes or when
 *      phases pop in and out of existence. Both of these operations
 *      change the species ordering.
 *        
 */
class VCS_SOLVE {
public:
  //! Constructor for the VCS_SOLVE class
  VCS_SOLVE();

  //! Destructor
  ~VCS_SOLVE();

  void InitSizes(int nspecies0, int nelements, int nphase0);

  //! Solve an equilibrium problem
  /*!
   *  This is the main interface routine to the equilibrium solver
   *
   * Input:
   *   @param vprob Object containing the equilibrium Problem statement
   *   
   *   @param ifunc Determines the operation to be done: Valid values:
   *            0 -> Solve a new problem by initializing structures
   *                 first. An initial estimate may or may not have
   *                 been already determined. This is indicated in the
   *                 VCS_PROB structure.
   *            1 -> The problem has already been initialized and 
   *                 set up. We call this routine to resolve it 
   *                 using the problem statement and
   *                 solution estimate contained in 
   *                 the VCS_PROB structure.
   *            2 -> Don't solve a problem. Destroy all the private 
   *                 structures.
   *
   *  @param ipr Printing of results
   *     ipr = 1 -> Print problem statement and final results to 
   *                standard output 
   *           0 -> don't report on anything 
   *  @param ip1 Printing of intermediate results
   *     IP1 = 1 -> Print intermediate results. 
   *
   *  @param maxit  Maximum number of iterations for the algorithm 
   *
   * Output:
   *
   *    @return
   *       nonzero value: failure to solve the problem at hand.
   *       zero : success
   */
  int vcs(VCS_PROB *vprob, int ifunc, int ipr, int ip1, int maxit);

  int vcs_solve_TP(int, int, int);

  void vcs_reinsert_deleted(int kspec);
  int vcs_basopt(int ifirst, double aw[], double sa[], double sm[], 
		 double ss[], double test, int *usedZeroedSpecies);

  //!  Choose a species for the next component
  /*!
   *   
   */
  int vcs_basisOptMax(const double *const x, const int j, const int n);

  int vcs_species_type(int kspec);
  void vcs_chemPotPhase(int iph, const double *const molNum, 
			double * const ac, double * const mu_i,
			bool do_deleted = false);
  void vcs_dfe(double *z, int kk, int ll, int lbot, int ltop);
  void vcs_updateVP(int place); 
  int vcs_RxnStepSizes(void);
  void vcs_tmoles(void);
  void vcs_deltag(int l, bool doDeleted);
  void vcs_switch_pos(int ifunc, int k1, int k2);
  void vcs_deltag_Phase(int iphase, bool doDeleted);

  //! birthGuess returns the number of moles of a species
  //! that is coming back to life or whose concentration has
  //! been forced to zero by a constraint for some reason, and needs
  //! to be reinitialized.
  /*!
   *  Do a minor alt calculation. But, cap the mole numbers at
   *  1.0E-15.
   *  For SS phases use VCS_DELETE_SPECIES_CUTOFF * 100.
   *
   *  The routine makes sure the guess doesn't reduce the concentration
   *  of a component by more than 1/3. Note this may mean that
   *  the vlaue coming back from this routine is zero or a 
   *  very small number.
   *
   *  @param kspec Species number that is coming back to life
   *  @return number of moles of the species
   */
  double vcs_birthGuess(int kspec);

  //! Solve an equilibrium problem at a particular fixed temperature 
  //! and pressure
  /*!
   *  The actual problem statement is assumed to be in the structure
   *  already.  This is a wrapper around the solve_TP() function.
   *  In this wrapper, we nondimensionalize the system 
   *  we calculate the standard state gibbs free energies of the
   *  species, and we decide whether to we need to use the
   *  initial guess algorithm.
   *
   * @param ipr = 1 -> Print results to standard output 
   *              0 -> don't report on anything 
   * @param ip1 = 1 -> Print intermediate results. 
   *              0 -> Dont print any intermediate results
   * @param maxit  Maximum number of iterations for the algorithm 
   * @param T    Value of the Temperature (Kelvin)
   * Param pres Value of the Pressure (units given by m_VCS_UnitsFormat variable
   *
   * @return Returns an integer representing the success of the algorithm
   *   0 = Equilibrium Achieved
   *   1 = Range space error encountered. The element abundance criteria are
   *       only partially satisfied. Specifically, the first NC= (number of
   *       components) conditions are satisfied. However, the full NE 
   *       (number of elements) conditions are not satisfied. The equilibrirum
   *       condition is returned.
   * -1 = Maximum number of iterations is exceeded. Convergence was not
   *      found.
   */
  int vcs_TP(int ipr, int ip1, int maxit, double T, double pres);

  int vcs_evalSS_TP(int ipr, int ip1, double Temp, double pres);
  void  vcs_fePrep_TP(void);
  double vcs_VolTotal(double, double, double [], double []);

  int vcs_prep_oneTime(int printLvl);

  //! Prepare the object for resolution
  /*!
   *  This routine is mostly concerned with changing the private data  
   *  to be consistent with that needed for solution. It is called for
   *  every invocation of the vcs_solve() except for the cleanup invocation.
   *
   * Tasks:
   *  1)  Initialization of arrays to zero.
   *
   *  return code
   *     VCS_SUCCESS = everything went OK
   *     VCS_PUB_BAD = There is an irreconcilable difference in the 
   *                   public data structure from when the problem was
   *                   initially set up.
   */
  int vcs_prep(void);
 
  bool vcs_wellPosed(VCS_PROB *vprob);

  int vcs_elem_rearrange(double *aw, double *sa, double *sm, double *ss);
  void vcs_switch_elem_pos(int ipos, int jpos);

  int    vcs_rxn_adj_cg(void);
  double vcs_Hessian_diag_adj(int, double);
  double vcs_Hessian_actCoeff_diag(int irxn);
  void vcs_CalcLnActCoeffJac(const double * const moleSpeciesVCS);
#ifdef DEBUG_MODE
  double vcs_line_search(int irxn, double dx_orig, char *ANOTE);
#else
  double vcs_line_search(int irxn, double dx_orig);
#endif

  int vcs_report(int);

  int vcs_rearrange(void);


  double vcs_nondim_Farad(int mu_units, double TKelvin);
  double vcs_nondimMult_TP(int mu_units, double TKelvin);
  void   vcs_nondim_TP(void);
  void   vcs_redim_TP(void);
  void   vcs_printChemPotUnits(int unitsFormat);

  void vcs_elab(void);
  int vcs_elabcheck(int ibound);
  void vcs_elabPhase(int iphase, double * const elemAbundPhase);
  int vcs_elcorr(double aa[], double x[]);

  int vcs_inest_TP(void);

#ifdef ALTLINPROG
  //! Extimate the initial mole numbers by constrained linear programming
  /*!
   *   This is done by running
   *   each reaction as far forward or backward as possible, subject
   *   to the constraint that all mole numbers remain
   *   non-negative. Reactions for which \f$ \Delta \mu^0 \f$ are
   *   positive are run in reverse, and ones for which it is negative
   *   are run in the forward direction. The end result is equivalent
   *   to solving the linear programming problem of minimizing the
   *   linear Gibbs function subject to the element and
   *   non-negativity constraints.
   */
  int vcs_setMolesLinProg();
#endif

  double vcs_Total_Gibbs(double *w, double *fe, double *tPhMoles);
  double vcs_GibbsPhase(int iphase, double *w, double *fe);

  double vcs_Gxs_phase_calc(vcs_VolPhase *Vphase, double *mf_PO);
  double vcs_Gxs_calc(int iphase);

  //! Transfer the results of the equilibrium calculation back to VCS_PROB
  /*!
   *   The VCS_PUB structure is returned to the user.
   *
   *  @param pub  Pointer to VCS_PROB object that will get the results of the
   *              equilibrium calculation transfered to it. 
   */
  int vcs_prob_update(VCS_PROB *pub);

  //! Fully specify the problem to be solved using VCS_PROB
  /*!
   *  Use the contents of the VCS_PROB to specify the contents of the
   *  private data, VCS_SOLVE.
   *
   *  @param pub  Pointer to VCS_PROB that will be used to
   *              initialize the current equilibrium problem
   */
  int vcs_prob_specifyFully(const VCS_PROB *pub);

  //! Specify the problem to be solved using VCS_PROB, incrementally
  /*!
   *  Use the contents of the VCS_PROB to specify the contents of the
   *  private data, VCS_SOLVE.
   *
   *  It's assumed we are solving the same problem.
   *
   *  @param pub  Pointer to VCS_PROB that will be used to
   *              initialize the current equilibrium problem
   */
  int VCS_SOLVE::vcs_prob_specify(const VCS_PROB *pub);
  
private:
  int zero_species(int kspec);
  int delete_species(int kspec);
  void delete_multiphase(int iph);
  int delta_species(int kspec, double *delta_ptr);
  void add_deleted(void);
  int recheck_deleted(void);

  //! Alternative treatment for the update of a minor species
  /*!
   * @param kspec Species index of the minor species
   * @param irxn  Rxn index of the same minor species
   * @param do_delete
   */
  double minor_alt_calc(int kspec, int irxn, int *do_delete
#ifdef DEBUG_MODE
			, char *ANOTE  
#endif
			);

  int force(int iti);
  int globStepDamp(int iti);
  void vcs_switch2D(double * const * const Jac, int k1, int k2);
  double l2normdg(double dg[]);
#ifdef DEBUG_MODE
  void prneav(void);
  void checkDelta1(double * const ds, double * const delTPhMoles, int kspec);
#endif
  void inest(double *aw, double *sa, double *sm, 
	     double *ss, double test);
  void vcs_SSPhase(void); 
  double deltaG_Recalc_Rxn(int irxn, const double *const molNum,
			   double * const ac, double * const mu_i);
  void delete_memory();

  //! Initialize the internal counters
  /*!
   * Initialize the internal counters containing the subroutine call
   * values and times spent in the subroutines.
   *
   *  ifunc = 0     Initialize only those counters appropriate for the top of
   *                vcs_solve_TP().
   *        = 1     Initialize all counters.
   */
  void vcs_counters_init(int ifunc);

  void vcs_TCounters_report(int timing_print_lvl = 1);

public:
  //! value of the number of species  used to malloc data structures
  int NSPECIES0;

  //! value of the number of phases  used to malloc data structures
  int NPHASE0;
  
  //!  Total number of species in the problems
  int m_numSpeciesTot;

  //! Number of element constraints in the problem
  /*! 
   * This is typically equal to the number of elements in the problem
   */
  int m_numElemConstraints;

  //! Number of components calculated for the problem
  int m_numComponents;

  //! Total number of non-component species in  the problem
  int m_numRxnTot;

  //! Current number of species in the problems
  /*!
   * Species can be deleted if they aren't
   * stable under the current conditions
   */
  int m_numSpeciesRdc;

  //! Current number of non-component species in the  problem 
  /*!
   * Species can be deleted if they aren't
   * stable under the current conditions
   */
  int m_numRxnRdc;

  //!  Number of active species which are  currently either zeroed out or
  //!  are minor species
  int m_numRxnMinorZeroed; 

  //! Number of Phases in the problem
  int NPhase;

  //! Formula matrix for the problem
  /*!
   *  FormulaMatrix[j][kspec] =  Number of elements, j, in the kspec species
   *
   *  Both element and species indecies are swapped.
   */
  DoubleStarStar FormulaMatrix;

  //! Stoichiometric coefficient matrix for the reaction mechanism 
  //! expressed in Reduced Canonical Form.
  /*!
   *   This is the stoichiometric coefficient matrix for the 
   *   reaction which forms species K from the component species. A
   *   stoichiometric coefficient of one is assumed for the 
   *   species K in this mechanism. 
   *
   *              NOTE: K = IRXN + NC
   *
   *   sc[irxn][j] :
   *     j refers to the component number, and irxn 
   *     refers to the irxn_th non-component species.
   *  
   *
   *   length = [nspecies0][nelements0]
   */
  //DoubleStarStar sc;
  DoubleStarStar m_stoichCoeffRxnMatrix;

  //! Absolute size of the stoichiometric coefficients
  /*!
   *  scSize[irxn] = abs(Size) of the stoichiometric
   *               coefficients. These are used to determine
   *               whether a given species should be
   *               handled by the alt_min treatment or 
   *               should be handled as a major species. 
   */
  std::vector<double> scSize;

  //! total size of the species
  /*!
   *  This is used as a multiplier to the mole number in figuring out which
   *  species should be components.
   */
  std::vector<double> m_spSize;

  //!  Standard state chemical potentials for species K at the current
  //!  temperature and pressure.
  /*!
   *  The first NC entries are for components. The following NR entries are
   *  for the current non-component species in the mechanism.
   */
  std::vector<double> m_SSfeSpecies;

  //! Dimensionless/Dimensional free energy for all the species in the mechanism at the 
  //! current T, P, and mole numbers. 
  /*!
   *   The first NC entries are for components. The following
   *   NR entries are for the current  non-component species in the mechanism.
   *  The dimension of this vector is specified by the m_VCS_UnitsFormat variable.
   *  Length = number of species
   */
  std::vector<double> m_feSpecies_curr;

  //! Free energy vector from the start of the current iteration
  /*!
   *  The free energies are saved at the start of the current iteration.
   *  Length = number of species  
   */
  std::vector<double> m_feSpecies_old;

  //! Dimensionless new free energy for all the species in the mechanism 
  //! at the new tentatite T, P, and mole numbers. 
  /*!
   *   The first NC entries are for components. The following
   *   NR entries are for the current  non-component species in the mechanism. 
   *  Length = number of species  
   */
  std::vector<double> m_feSpecies_new; 

  //! Setting for the initial estimate
  /*!
   *  Initial estimate: 0 user estimate
   *                   -1 machine estimate
   */
  int iest;

  //! Total moles of the species
  /*!
   *    Total number of moles of the kth species. 
   *    Length = Total number of species = m 
   */
  std::vector<double> m_molNumSpecies_old;

  //! Specifies the species unknown type
  /*!
   *   There are two types. One is the straightforward
   *   species, with the mole number w[k], as the
   *   unknown. The second is the an interfacial
   *   voltage where w[k] refers to the interfacial
   *   voltage in volts.
   *   These species types correspond to metalic
   *   electrons corresponding to electrodes.
   *   The voltage and other interfacial conditions
   *   sets up an interfacial current, which is
   *   set to zero in this initial treatment.
   *   Later we may have non-zero interfacial currents.
   */
  std::vector<int> SpeciesUnknownType; 

  //!  Change in the number of moles of phase,   iphase, due to the noncomponent formation
  //!  reaction, irxn, for species, k:
  /*!
   *    DnPhase[irxn][iphase]       =  k = nc + irxn 
   */
  DoubleStarStar DnPhase; 

  //!  This is 1 if the phase, iphase,  participates in the formation reaction
  //!   irxn, and zero otherwise.  PhaseParticipation[irxn][iphase]
  IntStarStar PhaseParticipation; 
 
  //! electric potential of the iph phase
  std::vector<double> phasePhi;

  //! Tentative value of the mole number vector. It's also used to store the
  //!     mole fraction vector.
  //std::vector<double> wt;
  std::vector<double> m_molNumSpecies_new;

  //! Delta G(I) for the noncomponent species  in  the mechanism.
  /*!
   *   Computed by the   subroutine  DELTAG. DG is the free
   *    energy change for the reaction which
   *              forms species K from the
   *              component species. This vector has length 
   *              equal to the number of noncomponent
   *              species in the mechanism. It starts with
   *              the first  current noncomponent species 
   *              in the mechanism. 
   */
  std::vector<double> m_deltaGRxn_new;   

  //!  Last deltag[irxn] from the previous step 
  std::vector<double> m_deltaGRxn_old;

  std::vector<double> m_deltaGRxn_tmp;

  //! Reaction Adjustments for each species during the current step
  /*!
   *  delta Moles for each species during the current step.
   *  Length = number of species
   */
  std::vector<double> m_deltaMolNumSpecies;


  std::vector<double> ga;  /* ga[j]      = Element abundances for jth element from 
			    *              estimate
			    *           -> this is calculated from the current mole 
			    *              fraction vector and BM, the formula 
			    *              vector. 
			    *              units = gmoles  */
  std::vector<double> gai; /* gai[j]     = Element abundances for jth element 
			    *              -> corrected
			    *              -> this is input from the input file and 
			    *                 is considered a constant from thereon.
			    *              units = gmoles  */
  double  TMoles;   /* TMoles      = Total number of moles in all phases
		     *               This number includes the inerts.
		     *            -> Don't use this except for scaling
		     *               purposes only                          */

  //! total gmols of species in each phase
  /*!
   *  This contains the total number of moles of species in each phase
   *
   *  Length = number of phases
   */
  std::vector<double> TPhMoles;

  //! total gmols of species in each phase in the tentative soln vector
  /*!
   *  This contains the total number of moles of species in each phase
   *  in the tentative solution vector
   *
   *  Length = number of phases
   */
  std::vector<double> TPhMoles1;

  //! Temporary vector of length NPhase 
  std::vector<double> TmpPhase;

  //! Temporary vector of length NPhase 
  std::vector<double> TmpPhase2;

  //! Change in the total moles in each phase
  /*!
   *  Length number of phases.
   */
  std::vector<double> DelTPhMoles;

  //! Temperature (Kelvin)
  double   T;

  //! Pressure (units are Pascals)
  double   Pres;

  //!  Total kmoles of inert to add to each phase 
  /*!
   *  TPhInertMoles[iph] = Total gmoles of  inert to add to each phase
   *  length = number of phases
   */
  std::vector<double> TPhInertMoles; 

  double   tolmaj;  /* tolmaj     = Tolerance requirement for major species */
  double   tolmin;  /* tolmin     = Tolerance requirement for minor species */
  double   tolmaj2; /* tolmaj2    = Below this, major species aren't refined
		     *              any more */
  double   tolmin2; /* tolmin2    = Below this, minor species aren't refined 
		     *              any more */
  std::vector<int> ind;    /* ind[k]     = Index vector that keeps track of the 
			    *              rearrangement
			    *              of the species vector within the problem.
			    *            -> At the end of each run, the species 
			    *              vector and associated data gets put back
			    *              in the original order. */

  //! Index that keeps track of the index of the species according
  //!  to the phase
  /*!
   *   indPhSp[k] = Index that keeps track of the index of the species according
   *   to the phase
   *  Length = number of species
   */
  std::vector<int> indPhSp;

  //! Index vector that keeps track of the rearrangement of the elements 
  /*!
   *  IndEl[j]
   */
  std::vector<int> IndEl;

  //!  Mapping between the species index for noncomponent species and the
  //!  full species  index.
  /*!
   * ir[irxn]   = Mapping between the species index for
   *              noncomponent species and the full species 
   *              index.
   *            - Initially set to a value of K = NC + I
   *              This vector has length equal to number 
   *              of noncomponent species in the mechanism. 
   *              It starts with the first current 
   *              noncomponent species in the mechanism.
   */
  std::vector<int> ir;

  //! Major - Minor status Vector for the noncomponent
  /*!
   *              species  irxn : 1 -> Major player  VCS_SPECIES_MAJOR
   *                              0 -> Minor player  VCS_SPECIES_MINOR
   *                             -1 -> Mole number is zero 
   *                                   in inactive phase VCS_SPECIES_ZEROEDPHASE
   *                             -2 -> Deleted species in an
   *                                   active multicom phase VCS_SPECIES_ZEROEDMS
   *                             -3 -> Mole number is zero
   *                                   in a stoich phase - VCS_SPECIES_ZEREODSS
   *                             -4 -> Species is deleted
   *                                   - VCS_SPECIES_DELETED
   *            -> Length equal to number of 
   *               non-components*/
  std::vector<int> spStatus;

  //!  Mapping from the species number to the phase number 
  std::vector<int> PhaseID;

  //!  Boolean indicating whether a species belongs to a single-species phase
  std::vector<int> SSPhase;


  //! Species string name for the kth species
  /*!
   *  SpName[k] = Species string name for the kth species
   */
  std::vector<std::string> SpName; 

  //! Vector of strings containing the element names
  /*!
   *   ElName[j]  = String containing element names      
   */
  std::vector<std::string> ElName; 

  //! Type of the element constraint
  /*!
   *  m_elType[j] = type of the element
   *             0  VCS_ELEM_TYPE_ABSPOS Normal element that is positive
   *                                     or zero in all species.
   *             1  VCS_ELEM_TPYE_ELECTRONCHARGE element dof that corresponds
   *                                        to the electronic charge DOF.
   *             2  VCS_ELEM_TYPE_CHARGENEUTRALITY element dof that
   *                              corresponds to a required charge
   *                              neutrality constraint on the phase.
   *                              The element abundance is always exactly zero.
   *             3  VCS_ELEM_TYPE_OTHERCONSTRAINT Other constraint which may
   *                              mean that a species has neg 0 or pos value 
   *                              of that constraint (other than charge)
   */
  std::vector<int> m_elType;

  //! Specifies whether an element constraint is active
  /*!
   * The default is true
   * Length = nelements
   */
  std::vector<int> ElActive;

  //! Array of Phase Structures
  /*!
   *  Length = number of phases
   */
  std::vector<vcs_VolPhase *> VPhaseList;

  //! String containing the title of the run
  std::string Title;

  //!   This specifies the current state of units for the Gibbs free energy
  //!   properties in the program.
  /*!
   *. The default is to have this unitless
   */
  char  UnitsState; 
 
  //! specifies the activity  convention of the phase containing the species
  /*!
   * SpecActConvention[kspec]
   *                 0 = molar based
   *                 1 = molality based                  
   *               length = number of species            
   */
  std::vector<int> SpecActConvention;

  //! specifies the activity convention of the phase.
  /*!
   *                 0 = molar based
   *                 1 = molality based                  
   *               length = number of phases          
   */
  std::vector<int> PhaseActConvention;
 
  //!  specifies the ln(Mnaught) used to   calculate the chemical potentials
  /*!
   *  For molar based activity conventions
   *  this will be equal to 0.0
   *   length = number of species            
   */
  std::vector<double> SpecLnMnaught;

  //! Activity Coefficients for Species
  /*!
   *
   * Length = number of species
   */
  std::vector<double> ActCoeff;

  //! Activity Coefficients for Species
  /*!
   *
   * Length = number of species
   */
  std::vector<double> ActCoeff0; 

  //! Change in activity coefficient with mole number
  /*!
   * length = [nspecies][nspecies]
   *
   *               (This is a temporary array that        
   *                gets regenerated every time it's      
   *                needed. It is not swapped wrt species 
   *  (unused atm)
   */
  DoubleStarStar dLnActCoeffdMolNum;

  //! This boolean indicates whether the   activity coefficients for a phase
  //!   are current.  
  std::vector<int> CurrPhAC;

  //! Molecular weight of each species
  /*!
   *  units = gm/gmol
   *  length = number of species
   */
  std::vector<double> WtSpecies;

  //! Charge of each species
  /*!
   * Length = number of species
   */
  std::vector<double> Charge;

  //! Vector of pointers to thermostructures which identify the model
  //! and parameters for evaluating the  thermodynamic functions for that 
  //! particular species.
  /*!
   * SpeciesThermo[k] pointer to the thermo information for the kth species
   */    
  std::vector<VCS_SPECIES_THERMO *> SpeciesThermo;

  //! Choice of Hessians
  /*!
   *  If this is true, then we will use a better approximation to the 
   *  Hessian based on Jacobian of the  ln(ActCoeff) with respect to mole 
   *  numbers                            
   */
  int UseActCoeffJac;

  double   Vol;     /* Vol        = Volume (cm^3) */
  
  //! Partialm molar volumes of the species
  /*!
   *  units = mks (m^3/kmol) -determined by m_VCS_UnitsFormat
   *  Length = number of species
   */
  std::vector<double> VolPM;

  //! dimensionless value of Faraday's constant 
  /*!
   *       F / RT  (1/volt)
   */
  double Faraday_dim;


  VCS_COUNTERS *m_VCount;

  int vcs_debug_print_lvl;

  //! printing level of timing information
  /*!
   *  1 allowing printing of timing
   *  0 do not allow printing of timing -> everything is printed
   *    as a NA.
   */
  int m_timing_print_lvl;

  //! Units for the chemical potential data:
  /*!
   *  VCS_UnitsFormat    = Units for the chemical potential data:
   *                          -1:  kcal/mol 
   *                           0:  MU/RT  
   *                           1:  kJ/mol 
   *                           2:  Kelvin 
   *                           3:  J / kmol 
   *                    and pressure data:
   *                          -1:  atm
   *                           0:  atm
   *                           1:  atm
   *                           2:  atm
   *                           3:  Pa
   */
  int m_VCS_UnitsFormat;

}; 

#ifdef ALTLINPROG
#else
int linprogmax(double *, double *, double *, double *, int, int, int);
#endif

}
#endif

