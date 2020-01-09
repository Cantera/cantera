/**
 *  @file ImplicitSurfChem.h
 * Declarations for the implicit integration of surface site density equations
 *  (see \ref  kineticsmgr and class
 *  \link Cantera::ImplicitSurfChem ImplicitSurfChem\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IMPSURFCHEM_H
#define CT_IMPSURFCHEM_H

#include "cantera/numerics/Integrator.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/kinetics/solveSP.h"

namespace Cantera
{

//! Advances the surface coverages of the associated set of SurfacePhase
//! objects in time
/*!
 * This function advances a set of SurfPhase objects, each associated with one
 * InterfaceKinetics object, in time. The following equation is used for each
 * surface phase, *i*.
 *
 *   \f[
 *        \dot \theta_k = \dot s_k (\sigma_k / s_0)
 *   \f]
 *
 * In this equation,
 * - \f$ \theta_k \f$ is the site coverage for the kth species.
 * - \f$ \dot s_k \f$ is the source term for the kth species
 * - \f$ \sigma_k \f$ is the number of surface sites covered by each species k.
 * - \f$ s_0 \f$ is the total site density of the interfacial phase.
 *
 * Additionally, the 0'th equation in the set is discarded. Instead the
 * alternate equation is solved for
 *
 * \f[
 *     \sum_{k=0}^{N-1}  \dot \theta_k = 0
 * \f]
 *
 * This last equation serves to ensure that sum of the \f$ \theta_k \f$ values
 * stays constant.
 *
 * The object uses the CVODE software to advance the surface equations.
 *
 * The solution vector used by this object is as follows: For each surface
 * phase with \f$ N_s \f$ surface sites, it consists of the surface coverages
 * \f$ \theta_k \f$ for \f$ k = 0, N_s - 1 \f$
 *
 * @ingroup  kineticsmgr
 */
class ImplicitSurfChem : public FuncEval
{
public:
    //! Constructor for multiple surfaces.
    /*!
     * @param k  Vector of pointers to InterfaceKinetics objects Each object
     *           consists of a surface or an edge containing internal degrees of
     *           freedom representing the concentration of surface adsorbates.
     * @param rtol   The relative tolerance for the integrator
     * @param atol   The absolute tolerance for the integrator
     * @param maxStepSize   The maximum step-size the integrator is allowed to take.
     *                      If zero, this option is disabled.
     * @param maxSteps   The maximum number of time-steps the integrator can take.
     *                   If not supplied, uses the default value in the CVodesIntegrator (20000).
     * @param maxErrTestFails   the maximum permissible number of error test failures
     *                           If not supplied, uses the default value in CVODES (7).
     */
    ImplicitSurfChem(std::vector<InterfaceKinetics*> k,
                     double rtol=1.e-7, double atol=1.e-14,
                     double maxStepSize=0, size_t maxSteps=20000,
                     size_t maxErrTestFails=7);

    virtual ~ImplicitSurfChem() {};

    /*!
     *  Must be called before calling method 'advance'
     */
    virtual void initialize(doublereal t0 = 0.0);

    /*!
     *  Set the maximum integration step-size.  Note, setting this value to zero
     *  disables this option
     */
    virtual void setMaxStepSize(double maxstep = 0.0);

    /*!
     *  Set the relative and absolute integration tolerances.
     */
    virtual void setTolerances(double rtol=1.e-7, double atol=1.e-14);

    /*!
     *  Set the maximum number of CVODES integration steps.
     */
    virtual void setMaxSteps(size_t maxsteps = 20000);

    /*!
     *  Set the maximum number of CVODES error test failures
     */
    virtual void setMaxErrTestFails(size_t maxErrTestFails = 7);

    //! Integrate from t0 to t1. The integrator is reinitialized first.
    /*!
     *   This routine does a time accurate solve from t = t0 to t = t1.
     *   of the surface problem.
     *
     *  @param t0  Initial Time -> this is an input
     *  @param t1  Final Time -> This is an input
     */
    void integrate(doublereal t0, doublereal t1);

    //! Integrate from t0 to t1 without reinitializing the integrator.
    /*!
     *  Use when the coverages have not changed from their values on return
     *  from the last call to integrate or integrate0.
     *
     *  @param t0  Initial Time -> this is an input
     *  @param t1  Final Time -> This is an input
     */
    void integrate0(doublereal t0, doublereal t1);

    //! Solve for the pseudo steady-state of the surface problem
    /*!
     * Solve for the steady state of the surface problem.
     * This is the same thing as the advanceCoverages() function,
     * but at infinite times.
     *
     * Note, a direct solve is carried out under the hood here,
     * to reduce the computational time.
     *
     * @param ifuncOverride One of the values defined in @ref solvesp_methods.
     *     The default is -1, which means that the program will decide.
     * @param timeScaleOverride When a pseudo transient is selected this value
     *             can be used to override the default time scale for
     *             integration which is one. When SFLUX_TRANSIENT is used, this
     *             is equal to the time over which the equations are integrated.
     *             When SFLUX_INITIALIZE is used, this is equal to the time used
     *             in the initial transient algorithm, before the equation
     *             system is solved directly.
     */
    void solvePseudoSteadyStateProblem(int ifuncOverride = -1,
                                       doublereal timeScaleOverride = 1.0);

    // overloaded methods of class FuncEval

    //! Return the number of equations
    virtual size_t neq() {
        return m_nv;
    }

    //! Evaluate the value of ydot[k] at the current conditions
    /*!
     *  @param t   Time (seconds)
     *  @param y   Vector containing the current solution vector
     *  @param ydot   Output vector containing the value of the
     *                derivative of the surface coverages.
     *  @param p   Unused parameter pass-through parameter vector
     */
    virtual void eval(doublereal t, doublereal* y, doublereal* ydot,
                      doublereal* p);

    //! Get the current state of the solution vector
    /*!
     *  @param y   Value of the solution vector to be used.
     *            On output, this contains the initial value
     *           of the solution.
     */
    virtual void getState(doublereal* y);

    /*!
     * Get the specifications for the problem from the values
     * in the ThermoPhase objects for all phases.
     *
     *  1. concentrations of all species in all phases, #m_concSpecies
     *  2. Temperature and pressure
     *
     *  @param vecConcSpecies Vector of concentrations. The phase concentration
     *                  vectors are contiguous within the object, in the same
     *                  order as the unknown vector.
     */
    void getConcSpecies(doublereal* const vecConcSpecies) const;

    //! Sets the concentrations within phases that are unknowns in
    //! the surface problem
    /*!
     * Fills the local concentration vector for all of the species in all of
     * the phases that are unknowns in the surface problem.
     *
     *  @param vecConcSpecies Vector of concentrations. The phase concentration
     *                  vectors are contiguous within the object, in the same
     *                  order as the unknown vector.
     */
    void setConcSpecies(const doublereal* const vecConcSpecies);

    //! Sets the state variable in all thermodynamic phases (surface and
    //! surrounding bulk phases) to the input temperature and pressure
    /*!
     *  @param TKelvin input temperature (kelvin)
     *  @param PresPa   input pressure in pascal.
     */
    void setCommonState_TP(doublereal TKelvin, doublereal PresPa);

    //! Returns a reference to the vector of pointers to the
    //! InterfaceKinetics objects
    /*!
     * This should probably go away in the future, as it opens up the class.
     */
    std::vector<InterfaceKinetics*> & getObjects() {
        return m_vecKinPtrs;
    }

    int checkMatch(std::vector<ThermoPhase*> m_vec, ThermoPhase* thPtr);

    void setIOFlag(int ioFlag) {
        m_ioFlag = ioFlag;
    }

protected:
    //! Set the mixture to a state consistent with solution
    //! vector y.
    /*!
     *  This function will set the surface site factions in the underlying
     *  SurfPhase objects to the current value of the solution vector.
     *
     * @param y Current value of the solution vector. The length is equal to
     *     the sum of the number of surface sites in all the surface phases.
     */
    void updateState(doublereal* y);

    //! vector of pointers to surface phases.
    std::vector<SurfPhase*> m_surf;

    //! Vector of pointers to bulk phases
    std::vector<ThermoPhase*> m_bulkPhases;

    //! vector of pointers to InterfaceKinetics objects
    std::vector<InterfaceKinetics*> m_vecKinPtrs;

    //! Vector of number of species in each Surface Phase
    std::vector<size_t> m_nsp;

    //! index of the surface phase in each InterfaceKinetics object
    std::vector<size_t> m_surfindex;

    std::vector<size_t> m_specStartIndex;

    //! Total number of surface species in all surface phases
    /*!
     * This is the total number of unknowns in m_mode 0 problem
     */
    size_t m_nv;

    size_t m_numTotalBulkSpecies;
    size_t m_numTotalSpecies;

    std::vector<vector_int> pLocVec;
    //! Pointer to the CVODE integrator
    std::unique_ptr<Integrator> m_integ;
    doublereal m_atol, m_rtol; // tolerances
    doublereal m_maxstep; //!< max step size
    size_t m_nmax; //!< maximum number of steps allowed
    size_t m_maxErrTestFails; //!< maximum number of error test failures allowed
    vector_fp m_work;

    /**
     * Temporary vector - length num species in the Kinetics object. This is
     * the sum of the number of species in each phase included in the kinetics
     * object.
     */
    vector_fp m_concSpecies;
    vector_fp m_concSpeciesSave;

    /**
     * Index into the species vector of the kinetics manager,
     * pointing to the first species from the surrounding medium.
     */
    int m_mediumSpeciesStart;
    /**
     * Index into the species vector of the kinetics manager, pointing to the
     * first species from the condensed phase of the particles.
     */
    int m_bulkSpeciesStart;
    /**
     * Index into the species vector of the kinetics manager, pointing to the
     * first species from the surface of the particles
     */
    int m_surfSpeciesStart;
    /**
     * Pointer to the helper method, Placid, which solves the surface problem.
     */
    std::unique_ptr<solveSP> m_surfSolver;

    //! If true, a common temperature and pressure for all surface and bulk
    //! phases associated with the surface problem is imposed
    bool m_commonTempPressForPhases;

    //! We make the solveSS class a friend because we need to access all of
    //! the above information directly. Adding the members into the class is
    //! also a possibility.
    friend class solveSS;

private:
    //! Controls the amount of printing from this routine
    //! and underlying routines.
    int m_ioFlag;
};
}

#endif
