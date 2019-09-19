/**
 * @file MultiPhase.h
 * Headers for the \link Cantera::MultiPhase MultiPhase\endlink
 * object that is used to set up multiphase equilibrium problems (see \ref equilfunctions).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MULTIPHASE_H
#define CT_MULTIPHASE_H

#include "cantera/numerics/DenseMatrix.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{

//! @defgroup equilfunctions

//! A class for multiphase mixtures. The mixture can contain any
//! number of phases of any type.
/*!
 * This object is the basic tool used by Cantera for use in Multiphase
 * equilibrium calculations.
 *
 * It is a container for a set of phases. Each phase has a given number of
 * kmoles. Therefore, MultiPhase may be considered an "extrinsic"
 * thermodynamic object, in contrast to the ThermoPhase object, which is an
 * "intrinsic" thermodynamic object.
 *
 * MultiPhase may be considered to be "upstream" of the ThermoPhase objects in
 * the sense that setting a property within MultiPhase, such as temperature,
 * pressure, or species mole number, affects the underlying ThermoPhase
 * object, but not the other way around.
 *
 * All phases have the same temperature and pressure, and a specified number
 * of moles for each phase. The phases do not need to have the same elements.
 * For example, a mixture might consist of a gaseous phase with elements (H,
 * C, O, N), a solid carbon phase containing only element C, etc. A master
 * element set will be constructed for the mixture that is the intersection of
 * the elements of each phase.
 *
 * Below, reference is made to global species and global elements. These refer
 * to the collective species and elements encompassing all of the phases
 * tracked by the object.
 *
 * The global element list kept by this object is an intersection of the
 * element lists of all the phases that comprise the MultiPhase.
 *
 * The global species list kept by this object is a concatenated list of all
 * of the species in all the phases that comprise the MultiPhase. The ordering
 * of species is contiguous with respect to the phase id.
 *
 * @ingroup equilfunctions
 */
class MultiPhase
{
public:
    //! Constructor.
    /*!
     *  The constructor takes no arguments, since phases are added using
     *  method addPhase().
     */
    MultiPhase();

    //! Destructor. Does nothing. Class MultiPhase does not take "ownership"
    //! (i.e. responsibility for destroying) the phase objects.
    virtual ~MultiPhase() {}

    //! Add a vector of phases to the mixture
    /*!
     * See the single addPhases command. This just does a bunch of phases
     * at one time
     *   @param phases Vector of pointers to phases
     *   @param phaseMoles Vector of mole numbers in each phase (kmol)
     */
    void addPhases(std::vector<ThermoPhase*>& phases, const vector_fp& phaseMoles);

    //! Add all phases present in 'mix' to this mixture.
    /*!
     *  @param mix  Add all of the phases in another MultiPhase
     *              object to the current object.
     */
    void addPhases(MultiPhase& mix);

    //! Add a phase to the mixture.
    /*!
     *  This function must be called before the init() function is called,
     *  which serves to freeze the MultiPhase.
     *
     *  @param p pointer to the phase object
     *  @param moles total number of moles of all species in this phase
     */
    void addPhase(ThermoPhase* p, doublereal moles);

    //! Number of elements.
    size_t nElements() const {
        return m_nel;
    }

    //! Check that the specified element index is in range.
    //! Throws an exception if m is greater than nElements()-1
    void checkElementIndex(size_t m) const;

    //! Check that an array size is at least nElements().
    //! Throws an exception if mm is less than nElements(). Used before calls
    //! which take an array pointer.
    void checkElementArraySize(size_t mm) const;

    //! Returns the name of the global element *m*.
    /*!
     *  @param m index of the global element
     */
    std::string elementName(size_t m) const;

    //! Returns the index of the element with name \a name.
    /*!
     * @param name   String name of the global element
     */
    size_t elementIndex(const std::string& name) const;

    //! Number of species, summed over all phases.
    size_t nSpecies() const {
        return m_nsp;
    }

    //! Check that the specified species index is in range.
    //! Throws an exception if k is greater than nSpecies()-1
    void checkSpeciesIndex(size_t k) const;

    //! Check that an array size is at least nSpecies().
    //! Throws an exception if kk is less than nSpecies(). Used before calls
    //! which take an array pointer.
    void checkSpeciesArraySize(size_t kk) const;

    //! Name of species with global index \a kGlob
    /*!
     * @param kGlob   global species index
     */
    std::string speciesName(const size_t kGlob) const;

    //! Returns the Number of atoms of global element \a mGlob in
    //! global species \a kGlob.
    /*!
     * @param kGlob   global species index
     * @param mGlob   global element index
     * @returns the number of atoms.
     */
    doublereal nAtoms(const size_t kGlob, const size_t mGlob) const;

    /// Returns the global Species mole fractions.
    /*!
     *   Write the array of species mole
     *   fractions into array \c x. The mole fractions are
     *   normalized to sum to one in each phase.
     *
     *    @param x  vector of mole fractions. Length = number of global species.
     */
    void getMoleFractions(doublereal* const x) const;

    //! Process phases and build atomic composition array.
    /*!
     *  This method must be called after all phases are added, before doing
     *  anything else with the mixture. After init() has been called, no more
     *  phases may be added.
     */
    void init();

    //! Returns the name of the n'th phase
    /*!
     *   @param iph  phase Index
     */
    std::string phaseName(const size_t iph) const;

    //! Returns the index, given the phase name
    /*!
     * @param pName Name of the phase
     * @returns the index. A value of -1 means the phase isn't in the object.
     */
    int phaseIndex(const std::string& pName) const;

    //! Return the number of moles in phase n.
    /*!
     * @param n  Index of the phase.
     */
    doublereal phaseMoles(const size_t n) const;

    //! Set the number of moles of phase with index n.
    /*!
     * @param n     Index of the phase
     * @param moles Number of moles in the phase (kmol)
     */
    void setPhaseMoles(const size_t n, const doublereal moles);

    /// Return a reference to phase n.
    /*!
     *  The state of phase n is also updated to match the state stored locally
     *  in the mixture object.
     *
     * @param n  Phase Index
     * @return   Reference to the ThermoPhase object for the phase
     */
    thermo_t& phase(size_t n);

    //! Check that the specified phase index is in range
    //! Throws an exception if m is greater than nPhases()
    void checkPhaseIndex(size_t m) const;

    //! Check that an array size is at least nPhases()
    //! Throws an exception if mm is less than nPhases(). Used before calls
    //! which take an array pointer.
    void checkPhaseArraySize(size_t mm) const;

    //! Returns the moles of global species \c k. units = kmol
    /*!
     * @param kGlob   Global species index k
     */
    doublereal speciesMoles(size_t kGlob) const;

    //! Return the global index of the species belonging to phase number \c p
    //! with local index \c k within the phase.
    /*!
     * @param k local index of the species within the phase
     * @param p index of the phase
     */
    size_t speciesIndex(size_t k, size_t p) const {
        return m_spstart[p] + k;
    }

    //! Return the global index of the species belonging to phase name \c phaseName
    //! with species name \c speciesName
    /*!
     * @param speciesName    Species Name
     * @param phaseName      Phase Name
     *
     * @returns the global index
     *
     * If the species or phase name is not recognized, this routine throws a
     * CanteraError.
     */
    size_t speciesIndex(const std::string& speciesName, const std::string& phaseName);

    /// Minimum temperature for which all solution phases have valid thermo
    /// data. Stoichiometric phases are not considered, since they may have
    /// thermo data only valid for conditions for which they are stable.
    doublereal minTemp() const {
        return m_Tmin;
    }

    /// Maximum temperature for which all solution phases have valid thermo
    /// data. Stoichiometric phases are not considered, since they may have
    /// thermo data only valid for conditions for which they are stable.
    doublereal maxTemp() const {
        return m_Tmax;
    }

    //! Total charge summed over all phases (Coulombs).
    doublereal charge() const;

    /// Charge (Coulombs) of phase with index \a p.
    /*!
     *  The net charge is computed as \f[ Q_p = N_p \sum_k F z_k X_k \f]
     *  where the sum runs only over species in phase \a p.
     *  @param p index of the phase for which the charge is desired.
     */
    doublereal phaseCharge(size_t p) const;

    //! Total moles of global element \a m, summed over all phases.
    /*!
     * @param m   Index of the global element
     */
    doublereal elementMoles(size_t m) const;

    //! Returns a vector of Chemical potentials.
    /*!
     * Write into array \a mu the chemical potentials of all species
     * [J/kmol]. The chemical potentials are related to the activities by
     *
     * \f$
     *          \mu_k = \mu_k^0(T, P) + RT \ln a_k.
     * \f$.
     *
     * @param mu Chemical potential vector. Length = num global species. Units
     *           = J/kmol.
     */
    void getChemPotentials(doublereal* mu) const;

    /// Returns a vector of Valid chemical potentials.
    /*!
     * Write into array \a mu the chemical potentials of all species with
     * thermo data valid for the current temperature [J/kmol]. For other
     * species, set the chemical potential to the value \a not_mu. If \a
     * standard is set to true, then the values returned are standard chemical
     * potentials.
     *
     * This method is designed for use in computing chemical equilibrium by
     * Gibbs minimization. For solution phases (more than one species), this
     * does the same thing as getChemPotentials. But for stoichiometric
     * phases, this writes into array \a mu the user-specified value \a not_mu
     * instead of the chemical potential if the temperature is outside the
     * range for which the thermo data for the one species in the phase are
     * valid. The need for this arises since many condensed phases have thermo
     * data fit only for the temperature range for which they are stable. For
     * example, in the NASA database, the fits for H2O(s) are only done up to
     * 0 C, the fits for H2O(L) are only done from 0 C to 100 C, etc. Using
     * the polynomial fits outside the range for which the fits were done can
     * result in spurious chemical potentials, and can lead to condensed
     * phases appearing when in fact they should be absent.
     *
     * By setting \a not_mu to a large positive value, it is possible to force
     * routines which seek to minimize the Gibbs free energy of the mixture to
     * zero out any phases outside the temperature range for which their
     * thermo data are valid.
     *
     * @param not_mu Value of the chemical potential to set species in phases,
     *               for which the thermo data is not valid
     * @param mu    Vector of chemical potentials. length = Global species,
     *              units = J kmol-1
     * @param standard  If this method is called with \a standard set to true,
     *                  then the composition-independent standard chemical
     *                  potentials are returned instead of the composition-
     *                  dependent chemical potentials.
     */
    void getValidChemPotentials(doublereal not_mu, doublereal* mu,
                                bool standard = false) const;

    //! Temperature [K].
    doublereal temperature() const {
        return m_temp;
    }

    //! Equilibrate a MultiPhase object
    /*!
     *  Set this mixture to chemical equilibrium by calling one of Cantera's
     *  equilibrium solvers. The XY parameter indicates what two thermodynamic
     *  quantities are to be held constant during the equilibration process.
     *
     *  @param XY      String representation of what two properties are being
     *                 held constant
     *  @param solver  Name of the solver to be used to equilibrate the phase.
     *      If solver = 'vcs', the vcs_MultiPhaseEquil solver will be used. If
     *      solver = 'gibbs', the MultiPhaseEquil solver will be used. If solver
     *      = 'auto', the 'vcs' solver will be tried first, followed by the
     *      'gibbs' solver if the first one fails.
     *  @param rtol      Relative tolerance
     *  @param max_steps Maximum number of steps to take to find the solution
     *  @param max_iter  The maximum number of outer temperature or pressure
     *      iterations to take when T and/or P is not held fixed.
     *  @param estimate_equil integer indicating whether the solver should
     *      estimate its own initial condition. If 0, the initial mole fraction
     *      vector in the ThermoPhase object is used as the initial condition.
     *      If 1, the initial mole fraction vector is used if the element
     *      abundances are satisfied. If -1, the initial mole fraction vector is
     *      thrown out, and an estimate is formulated.
     *  @param log_level  loglevel Controls amount of diagnostic output.
     *      log_level=0 suppresses diagnostics, and increasingly-verbose
     *      messages are written as loglevel increases.
     *
     * @ingroup equilfunctions
     */
    void equilibrate(const std::string& XY, const std::string& solver="auto",
                     double rtol=1e-9, int max_steps=50000, int max_iter=100,
                     int estimate_equil=0, int log_level=0);

    /// Set the temperature [K].
    /*!
     * @param T   value of the temperature (Kelvin)
     */
    void setTemperature(const doublereal T);

    //! Set the state of the underlying ThermoPhase objects in one call
    /*!
     * @param T      Temperature of the system (kelvin)
     * @param Pres   pressure of the system (pascal)
     */
    void setState_TP(const doublereal T, const doublereal Pres);

    //! Set the state of the underlying ThermoPhase objects in one call
    /*!
     * @param T      Temperature of the system (kelvin)
     * @param Pres   pressure of the system (pascal)
     * @param Moles  Vector of mole numbers of all the species in all the phases
     *               (kmol)
     */
    void setState_TPMoles(const doublereal T, const doublereal Pres, const doublereal* Moles);

    /// Pressure [Pa].
    doublereal pressure() const {
        return m_press;
    }

    /// The total mixture volume [m^3].
    /*!
     * Returns the cumulative sum of the volumes of all the phases in the
     * mixture.
     */
    doublereal volume() const;

    //! Set the pressure [Pa].
    /*!
     * @param P Set the pressure in the MultiPhase object (Pa)
     */
    void setPressure(doublereal P) {
        m_press = P;
        updatePhases();
    }

    //! The enthalpy of the mixture [J].
    doublereal enthalpy() const;

    //! The internal energy of the mixture [J].
    doublereal IntEnergy() const;

    //! The entropy of the mixture [J/K].
    doublereal entropy() const;

    //! The Gibbs function of the mixture [J].
    doublereal gibbs() const;

    //! Heat capacity at constant pressure [J/K]. Note that this does not
    //! account for changes in composition of the mixture with temperature.
    doublereal cp() const;

    //! Number of phases.
    size_t nPhases() const {
        return m_phase.size();
    }

    //! Return true is species \a kGlob is a species in a multicomponent
    //! solution phase.
    /*!
     * @param kGlob   index of the global species
     */
    bool solutionSpecies(size_t kGlob) const;

    //! Returns the phase index of the Kth "global" species
    /*!
     * @param kGlob Global species index.
     * @returns the index of the owning phase.
     */
    size_t speciesPhaseIndex(const size_t kGlob) const;

    //! Returns the mole fraction of global species k
    /*!
     * @param kGlob Index of the global species.
     */
    doublereal moleFraction(const size_t kGlob) const;

    //! Set the Mole fractions of the nth phase
    /*!
     * This function sets the mole fractions of the nth phase. Note, the mole
     * number of the phase stays constant
     *
     * @param n    index of the phase
     * @param x    Vector of input mole fractions.
     */
    void setPhaseMoleFractions(const size_t n, const doublereal* const x);

    //! Set the number of moles of species in the mixture
    /*!
     * @param xMap  CompositionMap of the species with nonzero mole numbers.
     *              Mole numbers that are less than or equal to zero will be
     *              set to zero. units = kmol.
     */
    void setMolesByName(const compositionMap& xMap);

    //! Set the moles via a string containing their names.
    /*!
     * The string x is in the form of a composition map. Species which are not
     * listed are set to zero.
     *
     * @param x string x in the form of a composition map
     *             where values are the moles of the species.
     */
    void setMolesByName(const std::string& x);

    //! Get the mole numbers of all species in the multiphase object
    /*!
     * @param[out] molNum Vector of doubles of length nSpecies containing the
     *               global mole numbers (kmol).
     */
    void getMoles(doublereal* molNum) const;

    //! Sets all of the global species mole numbers
    /*!
     * The state of each phase object is also updated to have the specified
     * composition and the mixture temperature and pressure.
     *
     * @param n    Vector of doubles of length nSpecies containing the global
     *             mole numbers (kmol).
     */
    void setMoles(const doublereal* n);

    //! Adds moles of a certain species to the mixture
    /*!
     * @param indexS      Index of the species in the MultiPhase object
     * @param addedMoles  Value of the moles that are added to the species.
     */
    void addSpeciesMoles(const int indexS, const doublereal addedMoles);

    //! Retrieves a vector of element abundances
    /*!
     * @param elemAbundances  Vector of element abundances
     * Length = number of elements in the MultiPhase object.
     * Index is the global element index. Units is in kmol.
     */
    void getElemAbundances(doublereal* elemAbundances) const;

    //! Return true if the phase \a p has valid thermo data for the current
    //! temperature.
    /*!
     * @param p  Index of the phase.
     */
    bool tempOK(size_t p) const;

    // These methods are meant for internal use.

    //! Update the locally-stored composition within this object to match the
    //! current compositions of the phase objects.
    /*!
     * Query the underlying ThermoPhase objects for their mole fractions and
     * fill in the mole fraction vector of this current object. Adjust element
     * compositions within this object to match.
     *
     * This is an upload operation in the sense that we are taking downstream
     * information (ThermoPhase object info) and applying it to an upstream
     * object (MultiPhase object).
     */
    void uploadMoleFractionsFromPhases();

    //! Set the states of the phase objects to the locally-stored
    //! state within this MultiPhase object.
    /*!
     * This method sets each phase to the mixture temperature and pressure,
     * and sets the phase mole fractions based on the mixture mole numbers.
     *
     * This is an download operation in the sense that we are taking upstream
     * object information (MultiPhase object) and applying it to downstream
     * objects (ThermoPhase object information)
     *
     * Therefore, the term, "update", is appropriate for a downstream operation.
     */
    void updatePhases() const;

private:
    //! Calculate the element abundance vector
    void calcElemAbundances() const;

    //! Set the mixture to a state of chemical equilibrium using the
    //! MultiPhaseEquil solver.
    /*!
     * @param XY   Integer flag specifying properties to hold fixed.
     * @param err  Error tolerance for \f$\Delta \mu/RT \f$ for all reactions.
     *             Also used as the relative error tolerance for the outer loop.
     * @param maxsteps Maximum number of steps to take in solving the fixed TP
     *                 problem.
     * @param maxiter Maximum number of "outer" iterations for problems holding
     *                fixed something other than (T,P).
     * @param loglevel Level of diagnostic output
     */
    double equilibrate_MultiPhaseEquil(int XY, doublereal err, int maxsteps,
                                       int maxiter, int loglevel);

    //! Vector of the number of moles in each phase.
    /*!
     * Length = m_np, number of phases.
     */
    vector_fp m_moles;

    //! Vector of the ThermoPhase pointers.
    std::vector<ThermoPhase*> m_phase;

    //! Global Stoichiometric Coefficient array
    /*!
     * This is a two dimensional array m_atoms(m, k). The first index is the
     * global element index. The second index, k, is the global species index.
     * The value is the number of atoms of type m in species k.
     */
    DenseMatrix m_atoms;

    //! Locally stored vector of mole fractions of all species comprising the
    //! MultiPhase object.
    vector_fp m_moleFractions;

    //! Mapping between the global species number and the phase ID
    /*!
     *  m_spphase[kGlobal] = iPhase
     *  Length = number of global species
     */
    std::vector<size_t> m_spphase;

    //! Vector of ints containing of first species index in the global list of
    //! species for each phase
    /*!
     *  kfirst = m_spstart[ip], kfirst is the index of the first species in
     *                          the ip'th phase.
     */
    std::vector<size_t> m_spstart;

    //! String names of the global elements. This has a length equal to the
    //! number of global elements.
    std::vector<std::string> m_enames;

    //! Atomic number of each global element.
    vector_int m_atomicNumber;

    //! Vector of species names in the problem. Vector is over all species
    //! defined in the object, the global species index.
    std::vector<std::string> m_snames;

    //! Returns the global element index, given the element string name
    /*!
     * -> used in the construction. However, wonder if it needs to be global.
     */
    std::map<std::string, size_t> m_enamemap;

    //! Current value of the temperature (kelvin)
    doublereal m_temp;

    //! Current value of the pressure (Pa)
    doublereal m_press;

    //! Number of distinct elements in all of the phases
    size_t m_nel;

    //! Number of distinct species in all of the phases
    size_t m_nsp;

    //! True if the init() routine has been called, and the MultiPhase frozen
    bool m_init;

    //! Global ID of the element corresponding to the electronic charge. If
    //! there is none, then this is equal to -1
    size_t m_eloc;

    //! Vector of bools indicating whether temperatures are ok for phases.
    /*!
     * If the current temperature is outside the range of valid temperatures
     * for the phase thermodynamics, the phase flag is set to false.
     */
    mutable std::vector<bool> m_temp_OK;

    //! Minimum temperature for which thermo parameterizations are valid.
    //! Stoichiometric phases are ignored in this determination. units Kelvin
    doublereal m_Tmin;

    //! Minimum temperature for which thermo parameterizations are valid.
    //! Stoichiometric phases are ignored in this determination. units Kelvin
    doublereal m_Tmax;

    //! Vector of element abundances
    /*!
     *  m_elemAbundances[mGlobal] = kmol of element mGlobal summed over all
     *      species in all phases.
     */
    mutable vector_fp m_elemAbundances;
};

//! Function to output a MultiPhase description to a stream
/*!
 *  Writes out a description of the contents of each phase of the
 *  MultiPhase using the report function.
 *
 *  @param s ostream
 *  @param x  Reference to a MultiPhase
 *  @returns a reference to the ostream
 */
inline std::ostream& operator<<(std::ostream& s, MultiPhase& x)
{
    x.updatePhases();
    for (size_t ip = 0; ip < x.nPhases(); ip++) {
        if (x.phase(ip).name() != "") {
            s << "*************** " << x.phase(ip).name() << " *****************" << std::endl;
        } else {
            s << "*************** Phase " << ip << " *****************" << std::endl;
        }
        s << "Moles: " << x.phaseMoles(ip) << std::endl;

        s << x.phase(ip).report() << std::endl;
    }
    return s;
}

//! Choose the optimum basis of species for the equilibrium calculations.
/*!
 * This is done by choosing the species with the largest mole fraction not
 * currently a linear combination of the previous components. Then, calculate
 * the stoichiometric coefficient matrix for that basis.
 *
 * Calculates the identity of the component species in the mechanism. Rearranges
 * the solution data to put the component data at the front of the species list.
 *
 * Then, calculates SC(J,I) the formation reactions for all noncomponent
 * species in the mechanism.
 *
 * @param[in] mphase   Pointer to the multiphase object. Contains the species
 *         mole fractions, which are used to pick the current optimal species
 *         component basis.
 * @param[in] orderVectorElements Order vector for the elements. The element
 *         rows in the formula matrix are rearranged according to this vector.
 * @param[in] orderVectorSpecies Order vector for the species. The species are
 *         rearranged according to this formula. The first nCompoments of this
 *         vector contain the calculated species components on exit.
 * @param[in] doFormRxn  If true, the routine calculates the formation
 *         reaction matrix based on the calculated component species. If
 *         false, this step is skipped.
 * @param[out] usedZeroedSpecies = If true, then a species with a zero
 *         concentration was used as a component. The problem may be converged.
 * @param[out] formRxnMatrix
 * @return      The number of components.
 *
 * @ingroup equilfunctions
 */
size_t BasisOptimize(int* usedZeroedSpecies, bool doFormRxn,
                     MultiPhase* mphase, std::vector<size_t>& orderVectorSpecies,
                     std::vector<size_t>& orderVectorElements,
                     vector_fp& formRxnMatrix);

//! Handles the potential rearrangement of the constraint equations
//! represented by the Formula Matrix.
/*!
 * Rearrangement is only necessary when the number of components is less
 * than the number of elements. For this case, some constraints can never
 * be satisfied exactly, because the range space represented by the Formula
 * Matrix of the components can't span the extra space. These constraints,
 * which are out of the range space of the component Formula matrix
 * entries, are migrated to the back of the Formula matrix.
 *
 * A prototypical example is an extra element column in FormulaMatrix[], which
 * is identically zero. For example, let's say that argon is has an element
 * column in FormulaMatrix[], but no species in the mechanism actually
 * contains argon. Then, nc < ne. Unless the entry for desired element
 * abundance vector for Ar is zero, then this element abundance constraint can
 * never be satisfied. The constraint vector is not in the range space of the
 * formula matrix.
 *
 * Also, without perturbation of FormulaMatrix[], BasisOptimize[] would
 * produce a zero pivot because the matrix would be singular (unless the argon
 * element column was already the last column of FormulaMatrix[].
 *
 * This routine borrows heavily from BasisOptimize algorithm. It finds nc
 * constraints which span the range space of the Component Formula matrix, and
 * assigns them as the first nc components in the formula matrix. This
 * guarantees that BasisOptimize has a nonsingular matrix to invert.
 *
 * @param[in] nComponents  Number of components calculated previously.
 * @param[in] elementAbundances  Current value of the element abundances
 * @param[in] mphase  Input pointer to a MultiPhase object
 * @param[in] orderVectorSpecies input vector containing the ordering of the
 *         global species in mphase. This is used to extract the component
 *         basis of the mphase object.
 * @param[out] orderVectorElements Output vector containing the order of the
 *         elements that is necessary for calculation of the formula matrix.
 *
 * @ingroup equilfunctions
 */
void ElemRearrange(size_t nComponents, const vector_fp& elementAbundances,
                   MultiPhase* mphase,
                   std::vector<size_t>& orderVectorSpecies,
                   std::vector<size_t>& orderVectorElements);

//! External int that is used to turn on debug printing for the
//! BasisOptimze program.
/*!
 *   Set this to 1 if you want debug printing from BasisOptimize.
 */
extern int BasisOptimize_print_lvl;
}

#endif
