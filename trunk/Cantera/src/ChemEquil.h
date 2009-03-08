/**
 *  @file ChemEquil.h
 *
 *  Chemical equilibrium.
 *
 *  $Author: hkmoffa $
 *  $Date: 2006/10/20 21:34:00 $
 *  $Revision: 1.19 $
 *
 *  Copyright 2001 California Institute of Technology
 *
 */


#ifndef CT_CHEM_EQUIL_H
#define CT_CHEM_EQUIL_H


// Cantera includes
#include "ct_defs.h"
#include "vec_functions.h"
#include "ctexceptions.h"
#include "ThermoPhase.h"
#include "DenseMatrix.h"

#include "MultiPhaseEquil.h"

namespace Cantera {

  int _equilflag(const char* xy);

  /**
   *  Chemical equilibrium options. Used internally by class ChemEquil.
   */
  class EquilOpt {
  public:
    EquilOpt() : relTolerance(1.e-8), absElemTol(1.0E-70),maxIterations(1000), 
		 iterations(0), 
		 maxStepSize(10.0), propertyPair(TP), contin(false) {}
        
    doublereal relTolerance;      ///< Relative tolerance
    doublereal absElemTol;        ///< Abs Tol in element number
    int maxIterations;            ///< Maximum number of iterations
    int iterations;               ///< Iteration counter

    /**
     * Maximum step size. Largest change in any element potential or
     * in log(T) allowed in one Newton step. Default: 10.0
     */
    doublereal maxStepSize;       

    /** 
     * Property pair flag. Determines which two thermodynamic properties
     * are fixed.
     */
    int propertyPair;

    /** 
     * Continuation flag. Set true if the calculation should be
     * initialized from the last calculation. Otherwise, the
     * calculation will be started from scratch and the initial
     * composition and element potentials estimated.
     */
    bool contin;
  };

  template<class M>
  class PropertyCalculator;

  /**
   * @defgroup equil Chemical Equilibrium
   * 
   */

  /**
   *  Class ChemEquil implements a chemical equilibrium solver for
   *  single-phase solutions. It is a "non-stoichiometric" solver in
   *  the terminology of Smith and Missen, meaning that every
   *  intermediate state is a valid chemical equilibrium state, but
   *  does not necessarily satisfy the element constraints. In
   *  contrast, the solver implemented in class MultiPhaseEquil uses
   *  a "stoichiometric" algorithm, in which each intermediate state
   *  satisfies the element constraints but is not a state of
   *  chemical equilibrium. Non-stoichiometric methods are faster
   *  when they converge, but stoichiometric ones tend to be more
   *  robust and can be used also for problems with multiple
   *  condensed phases.  As expected, the ChemEquil solver is faster
   *  than MultiPhaseEquil for many single-phase equilibrium
   *  problems (particularly if there are only a few elements but
   *  vvery many species), but can be less stable. Problem
   *  situations include low temperatures where only a few species
   *  have non-zero mole fractions, precisely stoichiometric
   *  compositions (e.g. 2 H2 + O2). In general, if speed is
   *  important, this solver should be tried first, and if it fails
   *  then use MultiPhaseEquil.
   * @ingroup equil
   */
  class ChemEquil {

  public:
    ChemEquil();
    virtual ~ChemEquil();

    int equilibrate(thermo_t& s, const char* XY,
		    bool useThermoPhaseElementPotentials = false);
    int equilibrate(thermo_t& s, const char* XY, vector_fp& elMoles,
		    bool useThermoPhaseElementPotentials = false);
    const vector_fp& elementPotentials() const { return m_lambda; }

    /**
     * Options controlling how the calculation is carried out. 
     * @see EquilOptions
     */
    EquilOpt options;


  protected:

    thermo_t*  m_phase;
    thermo_t* m_thermo;

    /// number of atoms of element m in species k.
    doublereal nAtoms(int k, int m) const { return m_comp[k*m_mm + m]; }

    void initialize(thermo_t& s);

    void setToEquilState(thermo_t& s, 
			 const vector_fp& x, doublereal t);

    int setInitialMoles(thermo_t& s);

    int estimateElementPotentials(thermo_t& s,  vector_fp& lambda);
        
    int estimateEP_Brinkley(thermo_t&s, vector_fp& lambda, vector_fp& elMoles);

    int dampStep(thermo_t& s, vector_fp& oldx, 
		 double oldf, vector_fp& grad, vector_fp& step, vector_fp& x, 
		 double& f, vector_fp& elmols, double xval, double yval );

    void equilResidual(thermo_t& s, const vector_fp& x, 
		       const vector_fp& elmtotal, vector_fp& resid, 
		       double xval, double yval);

    void equilJacobian(thermo_t& s, vector_fp& x,  
		       const vector_fp& elmols, DenseMatrix& jac, 
		       double xval, double yval);

    void update(const thermo_t& s);

    int m_mm;
    int m_kk;
    int m_skip;

    PropertyCalculator<thermo_t> *m_p1, *m_p2;

    /**
     * Current value of the mole fractions in the single phase.
     * -> length = m_kk.
     */
    vector_fp m_molefractions;
    /**
     * Current value of the dimensional element potentials
     * -> length = m_mm
     */
    vector_fp m_lambda;

    /*
     * Current value of the sum of the element abundances given the
     * current element potentials.
     */
    doublereal m_elementTotalSum;
    /*
     * Current value of the element mole fractions. Note these aren't
     * the goal element mole fractions.
     */
    vector_fp m_elementmolefracs;
    vector_fp m_reswork;
    vector_fp m_jwork1;
    vector_fp m_jwork2;
    /*
     * Storage of the element compositions
     *      natom(k,m) = m_comp[k*m_mm+ m];
     */
    vector_fp m_comp;
    doublereal m_temp, m_dens;
    doublereal m_p0;
    int m_eloc;
    doublereal m_abscharge;

    doublereal m_startTemp, m_startDens;
    vector_fp m_startSoln;

    vector_fp m_grt;
    vector_fp m_mu_RT;
    vector_int m_component;

    /*
     * element fractional cutoff, below which the element will be
     * zeroed. 
     */
    double m_elemFracCutoff;
    bool m_doResPerturb;
  };

}


#endif
