/**
 *  @file ChemEquil.h
 *
 *  Chemical equilibrium.
 */
/*
 *  Copyright 2001 California Institute of Technology
 */


#ifndef CT_CHEM_EQUIL_H
#define CT_CHEM_EQUIL_H


// Cantera includes
#include "cantera/base/ct_defs.h"
#include "cantera/base/vec_functions.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/numerics/DenseMatrix.h"

#include "MultiPhaseEquil.h"

#include <memory>

namespace Cantera
{

int _equilflag(const char* xy);

/**
 *  Chemical equilibrium options. Used internally by class ChemEquil.
 */
class EquilOpt
{
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
 *  very many species), but can be less stable. Problem
 *  situations include low temperatures where only a few species
 *  have non-zero mole fractions, precisely stoichiometric
 *  compositions (e.g. 2 H2 + O2). In general, if speed is
 *  important, this solver should be tried first, and if it fails
 *  then use MultiPhaseEquil.
 * @ingroup equil
 */
class ChemEquil
{

public:
    //! Default Constructor
    ChemEquil();

    //! Constructor combined with the initialization function
    /*!
     *  This constructor initializes the ChemEquil object with everything it
     *  needs to start solving equilibrium problems.
     *   @param s ThermoPhase object that will be used in the equilibrium calls.
     */
    ChemEquil(thermo_t& s);

    virtual ~ChemEquil();

    int equilibrate(thermo_t& s, const char* XY,
                    bool useThermoPhaseElementPotentials = false, int loglevel = 0);
    int equilibrate(thermo_t& s, const char* XY, vector_fp& elMoles,
                    bool useThermoPhaseElementPotentials = false, int loglevel = 0);
    const vector_fp& elementPotentials() const {
        return m_lambda;
    }

    /**
     * Options controlling how the calculation is carried out.
     * @see EquilOptions
     */
    EquilOpt options;


protected:

    //! Pointer to the %ThermoPhase object used to initialize this object.

    /*!
     *  This %ThermoPhase object must be compatible with the %ThermoPhase
     *  objects input from the equilibrate function. Currently, this
     *  means that the 2 %ThermoPhases have to have consist of the same
     *  species and elements.
     */
    thermo_t*  m_phase;

    /// number of atoms of element m in species k.
    doublereal nAtoms(size_t k, size_t m) const {
        return m_comp[k*m_mm + m];
    }

    void initialize(thermo_t& s);

    void setToEquilState(thermo_t& s,
                         const vector_fp& x, doublereal t);

    int setInitialMoles(thermo_t& s, vector_fp& elMoleGoal, int loglevel = 0);

    int estimateElementPotentials(thermo_t& s,  vector_fp& lambda,
                                  vector_fp& elMolesGoal, int loglevel = 0);

    int estimateEP_Brinkley(thermo_t& s, vector_fp& lambda, vector_fp& elMoles);

    int dampStep(thermo_t& s, vector_fp& oldx,
                 double oldf, vector_fp& grad, vector_fp& step, vector_fp& x,
                 double& f, vector_fp& elmols, double xval, double yval);

    void equilResidual(thermo_t& s, const vector_fp& x,
                       const vector_fp& elmtotal, vector_fp& resid,
                       double xval, double yval, int loglevel = 0);

    void equilJacobian(thermo_t& s, vector_fp& x,
                       const vector_fp& elmols, DenseMatrix& jac,
                       double xval, double yval, int loglevel = 0);

    void adjustEloc(thermo_t& s, vector_fp& elMolesGoal);

    void update(const thermo_t& s);

    double calcEmoles(thermo_t& s, vector_fp& x,
                      const double& n_t, const vector_fp& Xmol_i_calc,
                      vector_fp& eMolesCalc, vector_fp& n_i_calc,
                      double pressureConst);

    size_t m_mm;
    size_t m_kk;
    size_t m_skip;

    /**
     * This is equal to the rank of the stoichiometric coefficient
     * matrix when it is computed. It's initialized to m_mm.
     */
    size_t m_nComponents;

    std::auto_ptr<PropertyCalculator<thermo_t> > m_p1, m_p2;

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
    /**
     * Index of the element id corresponding to the electric charge of each
     * species. Equal to -1 if there is no such element id.
     */
    size_t m_eloc;

    vector_fp m_startSoln;

    vector_fp m_grt;
    vector_fp m_mu_RT;
    /**
     * Dimensionless values of the gibbs free energy for the
     * standard state of each species, at the temperature and
     * pressure of the solution (the star standard state).
     */
    vector_fp m_muSS_RT;
    std::vector<size_t> m_component;

    /*
     * element fractional cutoff, below which the element will be
     * zeroed.
     */
    double m_elemFracCutoff;
    bool m_doResPerturb;


    std::vector<size_t> m_orderVectorElements;
    std::vector<size_t> m_orderVectorSpecies;


};

extern int ChemEquil_print_lvl;

}

#endif
