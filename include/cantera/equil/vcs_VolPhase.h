/**
 * @file vcs_VolPhase.h
 *   Header for the object representing each phase within vcs
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef VCS_VOLPHASE_H
#define VCS_VOLPHASE_H

#include "cantera/equil/vcs_SpeciesProperties.h"
#include "cantera/base/Array.h"

namespace Cantera
{

class ThermoPhase;

//! Models for the standard state volume of each species
#define VCS_SSVOL_IDEALGAS 0
#define VCS_SSVOL_CONSTANT 1

/*
 * DEFINITIONS FOR THE vcs_VolPhase structure
 *
 * Equation of State Types
 * - Permissible values for the EqnState variable in CPC_PHASE structure
 */
#define VCS_EOS_CONSTANT 0
#define VCS_EOS_IDEAL_GAS 1
#define VCS_EOS_STOICH_SUB 5
#define VCS_EOS_IDEAL_SOLN 22
#define VCS_EOS_DEBEYE_HUCKEL 23
#define VCS_EOS_REDLICH_KWONG 24
#define VCS_EOS_REGULAR_SOLN 25
#define VCS_EOS_UNK_CANTERA -1

struct VCS_SPECIES;
class vcs_SpeciesProperties;
class VCS_SOLVE;

//!  Phase information and Phase calculations for vcs.
/*!
 * Each phase in a vcs calculation has a vcs_VolPhase object associated with it.
 * This object helps to coordinate property evaluations for species within the
 * phase. Usually these evaluations must be carried out on a per phase basis.
 * However, vcs frequently needs per species quantities. Therefore, we need an
 * interface layer between vcs and Cantera's ThermoPhase.
 *
 * The species stay in the same ordering within this structure. The vcs
 * algorithm will change the ordering of species in the global species list.
 * However, the indexing of species in this list stays the same. This structure
 * contains structures that point to the species belonging to this phase in the
 * global vcs species list.
 *
 * This object is considered not to own the underlying Cantera ThermoPhase
 * object for the phase.
 *
 * This object contains an idea of the temperature and pressure. It checks to
 * see if if the temperature and pressure has changed before calling underlying
 * property evaluation routines.
 *
 * The object contains values for the electric potential of a phase. It
 * coordinates the evaluation of properties wrt when the electric potential of a
 * phase has changed.
 *
 * The object knows about the mole fractions of the phase. It controls the
 * values of mole fractions, and coordinates the property evaluation wrt to
 * changes in the mole fractions. It also will keep track of the likely values
 * of mole fractions in multicomponent phases even when the phase doesn't
 * actually exist within the thermo program.
 *
 * The object knows about the total moles of a phase. It checks to see if the
 * phase currently exists or not, and modifies its behavior accordingly.
 *
 * Activity coefficients and volume calculations are lagged. They are only
 * called when they are needed (and when the state has changed so that they need
 * to be recalculated).
 */
class vcs_VolPhase
{
public:
    vcs_VolPhase(VCS_SOLVE* owningSolverObject = 0);

    vcs_VolPhase(const vcs_VolPhase& b) = delete;
    vcs_VolPhase& operator=(const vcs_VolPhase& b) = delete;
    ~vcs_VolPhase();

    //! The resize() function fills in all of the initial information if it
    //! is not given in the constructor.
    /*!
     * @param phaseNum    index of the phase in the vcs problem
     * @param numSpecies  Number of species in the phase
     * @param phaseName   String name for the phase
     * @param molesInert  kmoles of inert in the phase (defaults to zero)
     */
    void resize(const size_t phaseNum, const size_t numSpecies,
                const size_t numElem, const char* const phaseName,
                const double molesInert = 0.0);

    void elemResize(const size_t numElemConstraints);

    //! Set the moles and/or mole fractions within the phase
    /*!
     * @param molNum           total moles in the phase
     * @param moleFracVec      Vector of input mole fractions
     * @param vcsStateStatus   Status flag for this update
     */
    void setMoleFractionsState(const double molNum, const double* const moleFracVec,
                               const int vcsStateStatus);

    //! Set the moles within the phase
    /*!
     * This function takes as input the mole numbers in vcs format, and then
     * updates this object with their values. This is essentially a gather
     * routine.
     *
     * @param molesSpeciesVCS  Array of mole numbers. Note, the indices for
     *     species in this array may not be contiguous. IndSpecies[] is needed
     *     to gather the species into the local contiguous vector format.
     */
    void setMolesFromVCS(const int stateCalc,
                         const double* molesSpeciesVCS = 0);

    //! Set the moles within the phase
    /*!
     * This function takes as input the mole numbers in vcs format, and then
     * updates this object with their values. This is essentially a gather
     * routine.
     *
     * Additionally it checks to see that the total moles value in
     * TPhMoles[iplace] is equal to the internally computed value. If this isn't
     * the case, an error exit is carried out.
     *
     * @param vcsStateStatus  State calc value either `VCS_STATECALC_OLD` or
     *     `VCS_STATECALC_NEW`. With any other value nothing is done.
     * @param molesSpeciesVCS  array of mole numbers. Note, the indices for
     *     species in this array may not be contiguous. IndSpecies[] is needed
     *     to gather the species into the local contiguous vector format.
     * @param TPhMoles  VCS's array containing the number of moles in each phase.
     */
    void setMolesFromVCSCheck(const int vcsStateStatus,
                              const double* molesSpeciesVCS,
                              const double* const TPhMoles);

    //! Update the moles within the phase, if necessary
    /*!
     * This function takes as input the stateCalc value, which determines where
     * within VCS_SOLVE to fetch the mole numbers. It then updates this object
     * with their values. This is essentially a gather routine.
     *
     * @param stateCalc  State calc value either VCS_STATECALC_OLD or
     *     VCS_STATECALC_NEW. With any other value nothing is done.
     */
    void updateFromVCS_MoleNumbers(const int stateCalc);

    //! Fill in an activity coefficients vector within a VCS_SOLVE object
    /*!
     * This routine will calculate the activity coefficients for the current
     * phase, and fill in the corresponding entries in the VCS activity
     * coefficients vector.
     *
     * @param AC  vector of activity coefficients for all of the species in all
     *     of the phases in a VCS problem. Only the entries for the current
     *     phase are filled in.
     */
    void sendToVCS_ActCoeff(const int stateCalc, double* const AC);

    //! set the electric potential of the phase
    /*!
     * @param phi electric potential (volts)
     */
    void setElectricPotential(const double phi);

    //! Returns the electric field of the phase
    /*!
     *  Units are potential
     */
    double electricPotential() const;

    //! Gibbs free energy calculation for standard state of one species
    /*!
     * Calculate the Gibbs free energies for the standard state of the kth
     * species. The results are held internally within the object.
     *
     * @param kspec   Species number (within the phase)
     * @returns the Gibbs free energy for the standard state of the kth species.
     */
    double GStar_calc_one(size_t kspec) const;

    //! Gibbs free energy calculation at a temperature for the reference state
    //! of a species, return a value for one species
    /*!
     *  @param kspec   species index
     *  @returns       value of the Gibbs free energy
     */
    double G0_calc_one(size_t kspec) const;

    //! Molar volume calculation for standard state of one species
    /*!
     * Calculate the molar volume for the standard states. The results are held
     * internally within the object.
     *
     * @param kspec Species number (within the phase)
     * @return molar volume of the kspec species's standard state (m**3/kmol)
     */
    double VolStar_calc_one(size_t kspec) const;

    //! Fill in the partial molar volume vector for VCS
    /*!
     * This routine will calculate the partial molar volumes for the current
     * phase (if needed), and fill in the corresponding entries in the VCS
     * partial molar volumes vector.
     *
     * @param VolPM  vector of partial molar volumes for all of the species in
     *     all of the phases in a VCS problem. Only the entries for the current
     *     phase are filled in.
     */
    double sendToVCS_VolPM(double* const VolPM) const;

    //! Fill in the partial molar volume vector for VCS
    /*!
     * This routine will calculate the partial molar volumes for the
     * current phase (if needed), and fill in the corresponding entries in the
     * VCS partial molar volumes vector.
     *
     * @param VolPM  vector of partial molar volumes for all of the species in
     *     all of the phases in a VCS problem. Only the entries for the current
     *     phase are filled in.
     *
     * @todo This function's documentation is incorrect.
     */
    void sendToVCS_GStar(double* const gstar) const;

    //! Sets the temperature and pressure in this object and underlying
    //! ThermoPhase objects
    /*!
     * @param temperature_Kelvin    (Kelvin)
     * @param pressure_PA  Pressure (MKS units - Pascal)
     */
    void setState_TP(const double temperature_Kelvin, const double pressure_PA);

    //! Sets the temperature in this object and underlying ThermoPhase objects
    /*!
     * @param temperature_Kelvin    (Kelvin)
     */
    void setState_T(const double temperature_Kelvin);

    // Downloads the ln ActCoeff Jacobian into the VCS version of the
    // ln ActCoeff Jacobian.
    /*
     * This is essentially a scatter operation.
     *
     * @param LnAcJac_VCS Jacobian parameter
     *     The Jacobians are actually d( lnActCoeff) / d (MolNumber);
     *     dLnActCoeffdMolNumber(k,j)
     *
     *      j = id of the species mole number
     *      k = id of the species activity coefficient
     */
    void sendToVCS_LnActCoeffJac(Array2D& LnACJac_VCS);

    //! Set the pointer for Cantera's ThermoPhase parameter
    /*!
     * When we first initialize the ThermoPhase object, we read the state of the
     * ThermoPhase into vcs_VolPhase object.
     *
     * @param tp_ptr Pointer to the ThermoPhase object corresponding
     *               to this phase.
     */
    void setPtrThermoPhase(ThermoPhase* tp_ptr);

    //! Return a const ThermoPhase pointer corresponding to this phase
    /*!
     *  @return pointer to the ThermoPhase.
     */
    const ThermoPhase* ptrThermoPhase() const;

    //! Return the total moles in the phase
    double totalMoles() const;

    //! Returns the mole fraction of the kspec species
    /*!
     * @param kspec    Index of the species in the phase
     *
     * @return  Value of the mole fraction
     */
    double molefraction(size_t kspec) const;

    //! Sets the total moles in the phase
    /*!
     * We don't have to flag the internal state as changing here because we have
     * just changed the total moles.
     *
     * @param totalMols   Total moles in the phase (kmol)
     */
    void setTotalMoles(const double totalMols);

    //! Sets the mole flag within the object to out of date
    /*!
     * This will trigger the object to go get the current mole numbers when it
     * needs it.
     */
    void setMolesOutOfDate(int stateCalc = -1);

    //! Sets the mole flag within the object to be current
    void setMolesCurrent(int vcsStateStatus);

private:
    //! Set the mole fractions from a conventional mole fraction vector
    /*!
     * @param xmol Value of the mole fractions for the species in the phase.
     *             These are contiguous.
     */
    void setMoleFractions(const double* const xmol);

public:
    //! Return a const reference to the mole fractions stored in the object.
    const vector_fp & moleFractions() const;

    double moleFraction(size_t klocal) const;

    //! Sets the creationMoleNum's within the phase object
    /*!
     * @param F_k Pointer to a vector of n_k's
     */
    void setCreationMoleNumbers(const double* const n_k, const std::vector<size_t> &creationGlobalRxnNumbers);

    //! Return a const reference to the creationMoleNumbers stored in the object.
    /*!
     * @returns a const reference to the vector of creationMoleNumbers
     */
    const vector_fp & creationMoleNumbers(std::vector<size_t> &creationGlobalRxnNumbers) const;

    //! Returns whether the phase is an ideal solution phase
    bool isIdealSoln() const;

    //! Return the index of the species that represents the the voltage of the
    //! phase
    size_t phiVarIndex() const;

    void setPhiVarIndex(size_t phiVarIndex);

    //! Retrieve the kth Species structure for the species belonging to this
    //! phase
    /*!
     * The index into this vector is the species index within the phase.
     *
     * @param kindex kth species index.
     */
    vcs_SpeciesProperties* speciesProperty(const size_t kindex);

    //! int indicating whether the phase exists or not
    /*!
     * returns the m_existence int for the phase
     *
     * - VCS_PHASE_EXIST_ZEROEDPHASE = -6: Set to not exist by fiat from a
     *   higher level. This is used in phase stability boundary calculations
     * - VCS_PHASE_EXIST_NO = 0:   Doesn't exist currently
     * - VCS_PHASE_EXIST_MINORCONC = 1:  Exists, but the concentration is so low
     *   that an alternate method is used to calculate the total phase
     *   concentrations.
     * - VCS_PHASE_EXIST_YES = 2 : Does exist currently
     * - VCS_PHASE_EXIST_ALWAYS = 3: Always exists because it contains inerts
     *   which can't exist in any other phase. Or, the phase exists always
     *   because it consists of a single species, which is identified with the
     *   voltage, i.e., it's an electron metal phase.
     */
    int exists() const;

    //! Set the existence flag in the object
    /*!
     * Note the total moles of the phase must have been set appropriately before
     * calling this routine.
     *
     * @param existence Phase existence flag
     *
     * @note try to eliminate this routine
     */
    void setExistence(const int existence);

    //! Return the Global VCS index of the kth species in the phase
    /*!
     * @param spIndex local species index (0 to the number of species in the
     *                phase)
     *
     * @returns the VCS_SOLVE species index of the species. This changes as
     *         rearrangements are carried out.
     */
    size_t spGlobalIndexVCS(const size_t spIndex) const;


    //! set the Global VCS index of the kth species in the phase
    /*!
     * @param spIndex local species index (0 to the number of species
     *                in the phase)
     *
     * @returns the VCS_SOLVE species index of the that species This changes as
     *         rearrangements are carried out.
     */
    void setSpGlobalIndexVCS(const size_t spIndex, const size_t spGlobalIndex);

    //! Sets the total moles of inert in the phase
    /*!
     * @param tMolesInert Value of the total kmols of inert species in the
     *     phase.
     */
    void setTotalMolesInert(const double tMolesInert);

    //! Returns the value of the total kmol of inert in the phase
    double totalMolesInert() const;

    //! Returns the global index of the local element index for the phase
    size_t elemGlobalIndex(const size_t e) const;

    //! sets a local phase element to a global index value
    /*!
     * @param eLocal Local phase element index
     * @param eGlobal Global phase element index
     */
    void setElemGlobalIndex(const size_t eLocal, const size_t eGlobal);

    //! Returns the number of element constraints
    size_t nElemConstraints() const;

    //! Name of the element constraint with index \c e.
    /*!
     * @param e Element index.
     */
    std::string elementName(const size_t e) const;

    //! Type of the element constraint with index \c e.
    /*!
     * @param e Element index.
     */
    int elementType(const size_t e) const;

    //! Transfer all of the element information from the ThermoPhase object to
    //! the vcs_VolPhase object.
    /*!
     * Also decide whether we need a new charge neutrality element in the phase
     * to enforce a charge neutrality constraint.
     *
     * @param tPhase Pointer to the ThermoPhase object
     */
    size_t transferElementsFM(const ThermoPhase* const tPhase);

    //! Get a constant form of the Species Formula Matrix
    /*!
     *  Returns a `double**` pointer such that `fm[e][f]` is the formula
     *  matrix entry for element `e` for species `k`
     */
    const Array2D& getFormulaMatrix() const;

    //! Returns the type of the species unknown
    /*!
     * @param k species index
     * @return the SpeciesUnknownType[k] = type of species
     *      - Normal -> VCS_SPECIES_TYPE_MOLUNK (unknown is the mole number in
     *        the phase)
     *      - metal electron -> VCS_SPECIES_INTERFACIALVOLTAGE (unknown is the
     *        interfacial voltage (volts))
     */
    int speciesUnknownType(const size_t k) const;

    int elementActive(const size_t e) const;

    //! Return the number of species in the phase
    size_t nSpecies() const;

    //! Return the name corresponding to the equation of state
    std::string eos_name() const;

private:
    //! Evaluate the activity coefficients at the current conditions
    /*!
     * We carry out a calculation whenever #m_UpToDate_AC is false. Specifically
     * whenever a phase goes zero, we do not carry out calculations on it.
     */
    void _updateActCoeff() const;

    //! Gibbs free energy calculation for standard states
    /*!
     * Calculate the Gibbs free energies for the standard states The results are
     * held internally within the object.
     */
    void _updateGStar() const;

    //! Gibbs free energy calculation at a temperature for the reference state
    //! of each species
    void _updateG0() const;

    //! Molar volume calculation for standard states
    /*!
     * Calculate the molar volume for the standard states. The results are held
     * internally within the object. Units are in m**3/kmol.
     */
    void _updateVolStar() const;

    //! Calculate the partial molar volumes of all species and return the
    //! total volume
    /*!
     * Calculates these quantities internally and then stores them
     *
     * @return total volume [m^3]
     */
    double _updateVolPM() const;

    //! Evaluation of Activity Coefficient Jacobians
    /*!
     * This is the derivative of the ln of the activity coefficient with respect
     * to mole number of jth species. (temp, pressure, and other mole numbers
     * held constant)
     *
     * We employ a finite difference derivative approach here. Because we have
     * to change the mole numbers, this is not a const function, even though the
     * paradigm would say that it should be.
     */
    void _updateLnActCoeffJac();

    //! Updates the mole fraction dependencies
    /*!
     * Whenever the mole fractions change, this routine should be called.
     */
    void _updateMoleFractionDependencies();

private:
    //! Backtrack value of VCS_SOLVE *
    VCS_SOLVE* m_owningSolverObject;

public:
    //! Original ID of the phase in the problem.
    /*!
     * If a non-ideal phase splits into two due to a miscibility gap, these
     * numbers will stay the same after the split.
     */
    size_t VP_ID_;

    //! If true, this phase consists of a single species
    bool m_singleSpecies;

    //! If true, this phase is a gas-phase like phase
    /*!
     * A RTlog(p/1atm) term is added onto the chemical potential for inert
     * species if this is true.
     */
    bool m_gasPhase;

    //! Type of the equation of state
    /*!
     * The known types are listed at the top of this file.
     */
    int m_eqnState;

    //! This is the element number for the charge neutrality condition of the
    //! phase
    /*!
     *  If it has one.  If it does not have a charge neutrality
     * constraint, then this value is equal to -1
     */
    size_t ChargeNeutralityElement;

    //! Convention for the activity formulation
    /*!
     *  * 0 = molar based activities (default)
     *  * 1 = Molality based activities, mu = mu_0 + ln a_molality. Standard
     *    state is based on unity molality
     */
    int p_activityConvention;

private:
    //! Number of element constraints within the problem
    /*!
     *  This is usually equal to the number of elements.
     */
    size_t m_numElemConstraints;

    //! vector of strings containing the element constraint names
    /*!
     * Length =  nElemConstraints
     */
    std::vector<std::string> m_elementNames;

    //! boolean indicating whether an element constraint is active
    //! for the current  problem
    vector_int m_elementActive;

    //! Type of the element constraint
    /*!
     * m_elType[j] = type of the element:
     * * 0  VCS_ELEM_TYPE_ABSPOS Normal element that is positive or zero in
     *   all species.
     * * 1  VCS_ELEM_TYPE_ELECTRONCHARGE element dof that corresponds to the
     *   charge DOF.
     * * 2  VCS_ELEM_TYPE_OTHERCONSTRAINT Other constraint which may mean that
     *   a species has neg 0 or pos value of that constraint (other than
     *   charge)
     */
    vector_int m_elementType;

    //! Formula Matrix for the phase
    /*!
     *  FormulaMatrix(kspec,j) = Formula Matrix for the species
     *  Number of elements, j, in the kspec species
     */
    Array2D m_formulaMatrix;

    //! Type of the species unknown
    /*!
     *  SpeciesUnknownType[k] = type of species
     *  - Normal -> VCS_SPECIES_TYPE_MOLUNK.
     *    (unknown is the mole number in the phase)
     *  - metal electron -> VCS_SPECIES_INTERFACIALVOLTAGE.
     *    (unknown is the interfacial voltage (volts))
     */
    vector_int m_speciesUnknownType;

    //! Index of the element number in the global list of elements stored in VCS_SOLVE
    std::vector<size_t> m_elemGlobalIndex;

    //! Number of species in the phase
    size_t m_numSpecies;

public:
    //! String name for the phase
    std::string PhaseName;

private:
    //!  Total moles of inert in the phase
    double m_totalMolesInert;

    //! Boolean indicating whether the phase is an ideal solution
    //! and therefore its molar-based activity coefficients are
    //! uniformly equal to one.
    bool m_isIdealSoln;

    //! Current state of existence:
    /*!
     * - VCS_PHASE_EXIST_ZEROEDPHASE = -6: Set to not exist by fiat from a
     *   higher level. This is used in phase stability boundary calculations
     * - VCS_PHASE_EXIST_NO = 0:   Doesn't exist currently
     * - VCS_PHASE_EXIST_MINORCONC = 1:  Exists, but the concentration is so
     *   low that an alternate method is used to calculate the total phase
     *   concentrations.
     * - VCS_PHASE_EXIST_YES = 2 : Does exist currently
     * - VCS_PHASE_EXIST_ALWAYS = 3: Always exists because it contains inerts
     *   which can't exist in any other phase. Or, the phase exists always
     *   because it consists of a single species, which is identified with the
     *   voltage, i.e., its an electron metal phase.
     */
    int m_existence;

    // Index of the first MF species in the list of unknowns for this phase
    /*!
     *  This is always equal to zero.
     *  Am anticipating the case where the phase potential is species # 0,
     *  for multiphase phases. Right now we have the phase potential equal
     *  to 0 for single species phases, where we set by hand the mole fraction
     *  of species 0 to one.
     */
    int m_MFStartIndex;

    //! Index into the species vectors
    /*!
     *  Maps the phase species number into the global species number.
     *  Note, as part of the vcs algorithm, the order of the species
     *  vector is changed during the algorithm
     */
    std::vector<size_t> IndSpecies;

    //! Vector of Species structures for the species belonging to this phase
    /*!
     * The index into this vector is the species index within the phase.
     */
    std::vector<vcs_SpeciesProperties*> ListSpeciesPtr;

    /**
     *  If we are using Cantera, this is the pointer to the ThermoPhase
     *  object. If not, this is null.
     */
    ThermoPhase* TP_ptr;

    //!  Total mols in the phase. units are kmol
    double v_totalMoles;

    //! Vector of the current mole fractions for species in the phase
    vector_fp Xmol_;

    //! Vector of current creationMoleNumbers_
    /*!
     *  These are the actual unknowns in the phase stability problem
     */
    vector_fp creationMoleNumbers_;

    //! Vector of creation global reaction numbers for the phase stability problem
    /*!
     * The phase stability problem requires a global reaction number for each
     * species in the phase. Usually this is the krxn = kglob - M for species in
     * the phase that are not components. For component species, the choice of
     * the reaction is one which maximizes the chance that the phase pops into
     * (or remains in) existence.
     *
     * The index here is the local phase species index. the value of the
     * variable is the global vcs reaction number. Note, that the global
     * reaction number will go out of order when the species positions are
     * swapped. So, this number has to be recalculated.
     *
     * Length = number of species in phase
     */
    std::vector<size_t> creationGlobalRxnNumbers_;

    //! If the potential is a solution variable in VCS, it acts as a species.
    //! This is the species index in the phase for the potential
    size_t m_phiVarIndex;

    //! Total Volume of the phase. Units are m**3.
    mutable double m_totalVol;

    //! Vector of calculated SS0 chemical potentials for the
    //! current Temperature.
    /*!
     * Note, This is the chemical potential derived strictly from the polynomial
     * in temperature. Pressure effects have to be added in to get to the
     * standard state. Units are J/kmol.
     */
    mutable vector_fp SS0ChemicalPotential;

    //! Vector of calculated Star chemical potentials for the
    //! current Temperature and pressure.
    /*!
     * Note, This is the chemical potential at unit activity. Thus, we can call
     * it the standard state chemical potential as well. Units are J/kmol.
     */
    mutable vector_fp StarChemicalPotential;

    //! Vector of the Star molar Volumes of the species. units  m3 / kmol
    mutable vector_fp StarMolarVol;

    //! Vector of the Partial molar Volumes of the species. units  m3 / kmol
    mutable vector_fp PartialMolarVol;

    //! Vector of calculated activity coefficients for the current state
    /*!
     * Whether or not this vector is current is determined by the bool
     * #m_UpToDate_AC.
     */
    mutable vector_fp ActCoeff;

    //! Vector of the derivatives of the ln activity coefficient wrt to the
    //! current mole number multiplied by the current phase moles
    /*!
     * np_dLnActCoeffdMolNumber(k,j);
     * - j = id of the species mole number
     * - k = id of the species activity coefficient
     */
    mutable Array2D np_dLnActCoeffdMolNumber;

    //! Status
    /*!
     *  valid values are
     *  - VCS_STATECALC_OLD
     *  - VCS_STATECALC_NEW
     *  - VCS_STATECALC_TMP
     */
    int m_vcsStateStatus;

    //! Value of the potential for the phase (Volts)
    double m_phi;

    //! Boolean indicating whether the object has an up-to-date mole number vector
    //! and potential with respect to the current vcs state calc status
    bool m_UpToDate;

    //! Boolean indicating whether activity coefficients are up to date.
    /*!
     * Activity coefficients and volume calculations are lagged. They are only
     * called when they are needed (and when the state has changed so that they
     * need to be recalculated).
     */
    mutable bool m_UpToDate_AC;

    //! Boolean indicating whether Star volumes are up to date.
    /*!
     * Activity coefficients and volume calculations are lagged. They are only
     * called when they are needed (and when the state has changed so that they
     * need to be recalculated). Star volumes are sensitive to temperature and
     * pressure
     */
    mutable bool m_UpToDate_VolStar;

    //! Boolean indicating whether partial molar volumes are up to date.
    /*!
     * Activity coefficients and volume calculations are lagged. They are only
     * called when they are needed (and when the state has changed so that they
     * need to be recalculated). partial molar volumes are sensitive to
     * everything
     */
    mutable bool m_UpToDate_VolPM;

    //! Boolean indicating whether GStar is up to date.
    /*!
     * GStar is sensitive to the temperature and the pressure, only
     */
    mutable bool m_UpToDate_GStar;

    //! Boolean indicating whether G0 is up to date.
    /*!
     * G0 is sensitive to the temperature and the pressure, only
     */
    mutable bool m_UpToDate_G0;

    //! Current value of the temperature for this object, and underlying objects
    double Temp_;

    //! Current value of the pressure for this object, and underlying objects
    double Pres_;
};

}

#endif
