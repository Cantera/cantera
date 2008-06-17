/**
 * @file vcs_VolPhase.h  
 *   Header for the object representing each phase within vcs
 */
/*
 * $Id$ 
 */
/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef VCS_VOLPHASE_H
#define VCS_VOLPHASE_H

#include "vcs_DoubleStarStar.h"

#include <vector>
#include <string>

/*
 * Forward references
 */
// Forward reference for ThermoPhase object within the Cantera namespace
namespace Cantera {
  class ThermoPhase;
}

namespace VCSnonideal {
  /*
   *     Models for the species activity coefficients
   *
   */
#define VCS_AC_CONSTANT       0
  //#define VCS_AC_DEBYE_HUCKEL   23
  //#define VCS_AC_REGULAR_SOLN   25
  //#define VCS_AC_MARGULES       300
#define VCS_AC_UNK_CANTERA    -1
#define VCS_AC_UNK            -2
  /*
   *
   *    Models for the standard state volume of each species
   */
#define VCS_SSVOL_IDEALGAS    0
#define VCS_SSVOL_CONSTANT    1

  /*
   *            DEFINITIONS FOR THE vcs_VolPhase structure
   *
   *
   *   Equation of State Types
   *        - Permissible values for the EqnState variable in CPC_PHASE structure
   */
#define VCS_EOS_CONSTANT      0
#define VCS_EOS_IDEAL_GAS     1
#define VCS_EOS_STOICH_SUB    5
#define VCS_EOS_IDEAL_SOLN    22
#define VCS_EOS_DEBEYE_HUCKEL 23
#define VCS_EOS_REDLICK_KWONG 24
#define VCS_EOS_REGULAR_SOLN  25
#define VCS_EOS_UNK_CANTERA   -1


  struct VCS_SPECIES;
  class vcs_SpeciesProperties;
  class VCS_SOLVE;


  //!  Phase information and Phase calculations for vcs.
  /*!
   *  Each phase in a vcs calculation has a vcs_VolPhase object associated
   *  with it. This object helps to coordinate property evaluations for
   *  species within the phase. Usually these evaluations must be carried
   *  out on a per phase basis. However, vcs frequently needs per species
   *  quantitites. Therefore, we need an interface layer between vcs
   *  and Cantera's ThermoPhase. 
   *
   *  The species stay in the same ordering within this structure.
   *  The vcs algorithm will change the ordering of species in 
   *  the global species list. However, the indexing of species in this
   *  list stays the same. This structure contains structures that
   *  point to the species belonging to this phase in the global
   *  vcs species list.
   *
   * This object is considered not to own the underlying Cantera ThermoPhase
   * object for the phase.
   *
   * This object contains an idea of the temperature and pressure.
   * It checks to see if if the temperature and pressure has changed before calling
   * underlying property evalulation routines.
   *
   * The object contains values for the electric potential of a phase.
   * It coordinates the evalulation of properties wrt when the electric
   * potential of a phase has changed.
   *
   * The object knows about the mole fractions of the phase. It controls
   * the values of mole fractions, and coordinates the property evalulation
   * wrt to changes in the mole fractions. It also will keep track of the
   * likely values of mole fractions in multicomponent phases even when
   * the phase doesn't actually exist within the thermo program.
   *
   * The object knows about the total moles of a phase. It checkes to 
   * see if the phase currently exists or not, and modifies its behavior
   * accordingly.
   * 
   *
   * Activity coefficients and volume calculations are lagged. They are only
   * called when they are needed (and when the state has changed so that they
   * need to be recalculated).
   */
  class vcs_VolPhase {
  public:

    /*************************************************************************
     *              FUNCTIONS                                                *
     ************************************************************************/

    //! Base constructor for the class
    vcs_VolPhase(VCS_SOLVE * owningSolverObject = 0);

    //! Copy constructor
    /*!
     *  @param b object to be copied
     */
    vcs_VolPhase(const vcs_VolPhase& b);

    //! Assignment operator
    /*!
     *  @param b object to be copied
     */
    vcs_VolPhase& operator=(const vcs_VolPhase& b);

    //! Destructor
    ~vcs_VolPhase();

  
    //! The resize() function fills in all of the initial information if it
    //! is not given in the constructor.
    /*!
     *  @param phaseNum    index of the phase in the vcs problem
     *  @param numSpecies  Number of species in the phase
     *  @param phaseName   String name for the phase
     *  @param molesInert  kmoles of inert in the phase (defaults to zero)
     */
    void resize(int phaseNum, int numSpecies, const char *phaseName,
		double molesInert = 0.0);

  private:
    //! Evaluate activity coefficients
    /*!
     *   We carry out a calculation whenever UpTODate_AC is false. Specifically
     *   whenever a phase goes zero, we do not carry out calculations on it.
     */
    void evaluateActCoeff() const;

  public:
    //! Evaluate activity coefficients and return the kspec coefficient
    /*!
     *   We carry out a calculation whenever UpTODate_AC is false. Specifically
     *   whenever a phase goes zero, we do not carry out calculations on it.
     *
     * @param kspec species number
     */
    double AC_calc_one(int kspec) const;
  
    //! Set the moles within the phase
    /*!
     *  This function takes as input the mole numbers in vcs format, and
     *  then updates this object with their values. This is essentially
     *  a gather routine.
     *
     *  @param molesSpeciesVCS  array of mole numbers. Note, the indecises for species in 
     *            this array may not be contiguous. IndSpecies[] is needed
     *            to gather the species into the local contiguous vector
     *            format. 
     */
    void setMolesFromVCS(const int stateCalc, 
			 const double * const molesSpeciesVCS = 0);

    //! Set the moles within the phase
    /*!
     *  This function takes as input the mole numbers in vcs format, and
     *  then updates this object with their values. This is essentially
     *  a gather routine.
     *  Additionally it checks to see that the total moles value in 
     *  TPhMoles[iplace] is equal to the internally computed value.
     *  If this isn't the case, an error exit is carried out.
     *
     *  
     *  @param molesSpeciesVCS  array of mole numbers. Note, the indecises
     *            for species in 
     *            this array may not be contiguous. IndSpecies[] is needed
     *            to gather the species into the local contiguous vector
     *            format. 
     *  @param TPhMoles   VCS's array containing the number of moles
     *                    in each phase.
     *  @param iphase     index of the current phase.
     *
     */
    void setMolesFromVCSCheck(const int stateCalc,
			      const double * molesSpeciesVCS,
			      const double * const TPhMoles);

    //! Update the moles within the phase, if necessary
    /*!
     *  This function takes as input the stateCalc value, which 
     *  determines where within VCS_SOLVE to fetch the mole numbers.
     *  It then updates this object with their values. This is essentially
     *  a gather routine.
     *
     *  @param stateCalc    State calc value either VCS_STATECALC_OLD 
     *                      or  VCS_STATECALC_NEW. With any other value
     *                      nothing is done.
     *
     */
    void updateFromVCS_MoleNumbers(const int stateCalc);

    //! Fill in an activity coefficients vector within a VCS_SOLVE object
    /*!
     *  This routine will calculate the activity coefficients for the
     *  current phase, and fill in the corresponding entries in the
     *  VCS activity coefficients vector.
     *  
     * @param AC  vector of activity coefficients for all of the species
     *            in all of the phases in a VCS problem. Only the
     *            entries for the current phase are filled in.
     */
    void sendToVCS_ActCoeff(const int stateCalc, double * const AC);

    //! set the electric potential of the phase
    /*!
     * @param phi electric potential (volts)
     */
    void setElectricPotential(double phi);

    //! Returns the electric field of the phase
    /*!
     *  Units are potential
     */
    double electricPotential() const;
  
  private:
    //! Gibbs free energy calculation for standard states
    /*!
     * Calculate the Gibbs free energies for the standard states
     * The results are held internally within the object.
     *
     * @param TKelvin Current temperature
     * @param pres    Current pressure
     */
    void GStar_calc() const;

  public:
    //! Gibbs free energy calculation for standard state of one species
    /*!
     * Calculate the Gibbs free energies for the standard state
     * of the kth species.
     * The results are held internally within the object.
     * The kth species standard state G is returned
     *
     * @param kspec   Species number (within the phase)
     * @param TKelvin Current temperature
     * @param pres    Current pressure
     *
     * @return Gstar[kspec] returns the gibbs free energy for the
     *         standard state of the kth species.
     */
    double GStar_calc_one(int kspec);

    //! Gibbs free energy calculation at a temperature for the reference state
    //! of each species
    /*!
     *  @param TKelvin temperature
     */
    void G0_calc(double TKelvin);

    //! Gibbs free energy calculation at a temperature for the reference state
    //! of a species, return a value for one species
    /*!
     *  @param kspec   species index
     *  @param TKelvin temperature
     *
     *  @return return value of the gibbs free energy
     */
    double G0_calc_one(int kspec, double TKelvin);

  
  private:
    //! Molar volume calculation for standard states
    /*!
     * Calculate the molar volume for the standard states
     * The results are held internally within the object.
     *
     * @param TKelvin Current temperature
     * @param pres    Current pressure
     *
     *  Units are in m**3/kmol
     */
    void VolStar_calc() const;

  public:
    //! Molar volume calculation for standard state of one species
    /*!
     * Calculate the molar volume for the standard states
     * The results are held internally within the object.
     * Return the molar volume for one species
     *
     * @param kspec Species number (within the phase)
     * @param TKelvin Current temperature
     * @param pres    Current pressure
     *
     * @return molar volume of the kspec species's standard
     *         state (m**3/kmol)
     */
    double VolStar_calc_one(int kglob, double TKelvin, double pres);

  private:
    //! Calculate the partial molar volumes of all species and return the
    //! total volume
    /*!
     *  Calculates these quantitites internally
     *
     * @return total volume
     */
    double VolPM_calc() const;

  public:
    //! Fill in the partial molar volume vector for VCS
    /*!
     *  This routine will calculate the partial molar volumes for the
     *  current phase (if needed), and fill in the corresponding entries in the
     *  VCS partial molar volumes vector.
     *  
     * @param VolPM  vector of partial molar volumes for all of the species
     *            in all of the phases in a VCS problem. Only the
     *            entries for the current phase are filled in.
     */
    double sendToVCS_VolPM(double * const VolPM) const;

    //! Fill in the partial molar volume vector for VCS
    /*!
     *  This routine will calculate the partial molar volumes for the
     *  current phase (if needed), and fill in the corresponding entries in the
     *  VCS partial molar volumes vector.
     *  
     * @param VolPM  vector of partial molar volumes for all of the species
     *            in all of the phases in a VCS problem. Only the
     *            entries for the current phase are filled in.
     */
    void sendToVCS_GStar(double * const gstar);

    //! Sets the temperature and pressure in this object and
    //! underlying objects
    /*!
     *  Sets the temperature and pressure in this object and
     *  underlying objects. The underlying objects refers to the
     *  Cantera's ThermoPhase object for this phase.
     *
     *  @param temperature_Kelvin    (Kelvin)
     *  @param pressure_PA  Pressure (MKS units - Pascal)
     */
    void setState_TP(double temperature_Kelvin, double pressure_PA);

  private:
    //! Evaluation of Activity Coefficient Jacobians
    /*!
     *  This is the derivative of the ln of the activity coefficient
     *  with respect to mole number of jth species.
     *  (temp, pressure, and other mole numbers held constant)
     *
     *  We employ a finite difference derivative approach here.
     *  Because we have to change the mole numbers, this is not
     *  a const function, even though the paradigm would say that
     *  it should be.
     *
     *  @param moleNumbers Mole numbers are input.
     */
    void updateLnActCoeffJac();
 
  public:
    // Downloads the ln ActCoeff jacobian into the VCS version of the
    // ln ActCoeff jacobian.
    /*
     *
     *   This is essentially a scatter operation.
     *
     *  @param LnAcJac_VCS jacobian parameter
     *   The Jacobians are actually d( lnActCoeff) / d (MolNumber);
     *   dLnActCoeffdMolNumber[j][k]
     * 
     *      j = id of the species mole number
     *      k = id of the species activity coefficient
     */
    void sendToVCS_LnActCoeffJac(double * const * const LnACJac_VCS);

    //! Set the pointer for Cantera's ThermoPhase parameter
    /*!
     *  When we first initialize the ThermoPhase object, we read the
     *  state of the ThermoPhase into vcs_VolPhase object.
     *
     * @param tp_ptr Pointer to the ThermoPhase object corresponding
     *               to this phase.
     */
    void setPtrThermoPhase(Cantera::ThermoPhase *tp_ptr);

    //! Return a const ThermoPhase pointer corresponding to this phase
    /*!
     *  @return pointer to the ThermoPhase.
     */
    const Cantera::ThermoPhase *ptrThermoPhase() const;

    //! Return the total moles in the phase
    /*! 
     *
     *  Units -> depends on VCS_UnitsFormat variable
     *             Cantera -> J/kmol
     */
    double TotalMoles() const;

    //! Returns the mole fraction of the kspec species
    /*!
     *  @param kspec    Index of the species in the phase
     *
     * @return  Value of the mole fraction
     */
    double molefraction(int kspec) const;

    //! Sets the total moles in the phase
    /*!
     *  @param totalMols   Total moles in the phase (kmol)
     */
    void setTotalMoles(double totalMols);

    //! Set the mole fractions from a conventional mole fraction vector
    /*!
     *
     * @param xmol Value of the mole fractions for the species
     *             in the phase. These are contiguous. 
     */
    void setMoleFractions (const double * const xmol);

    //! Return a const reference to the mole fractions
    const std::vector<double> & moleFractions() const;

    //! Returns whether the phase is an ideal solution phase
    bool isIdealSoln() const;

    //! Returns whether the object is using cantera calls.
    bool usingCanteraCalls() const;

  private:

    //! Updates the mole fraction depenpencies
    /*!
     *  Whenever the mole fractions change, this routine
     *  should be called.
     */
    void _updateMoleFractionDependencies();

    /*************************************************************************
     *         MEMBER  DATA                                                  *
     ************************************************************************/

  public:
    //! Backtrack value of VCS_SOLVE *
    /*!
     *  Note the default for this is 0. That's a valid value too, since
     *  VCS_PROB also uses vcs_VolPhase objects.
     */
    VCS_SOLVE *m_owningSolverObject;

    //! Original ID of the phase in the problem. 
    /*!
     * If a non-ideal phase splits into two due to a
     *  miscibility gap, these numbers will stay the
     * same after the split.                          
     */
    int VP_ID;

    //! ID of the surface or volume domain in which the
    //!  this phase exists
    /*!
     * This ventures into the idea of installing a physical location
     * into a thermodynamics program. This unknown is currently not
     * being used.
     */
    int Domain_ID; 

    //! If true, this phase consists of a single species
    int SingleSpecies;

    //! If true, this phase is a gas-phase like phase
    /*!
     * A RTlog(p/1atm) term is added onto the chemical potential for inert
     * species if this is true.
     */
    bool m_gasPhase;        

    //! Type of the equation of state
    /*!
     *  The known types are listed at the top of this file.
     */
    int EqnState;

    //! Number of element constraints within the problem
    /*!
     *  This is usually equal to the number of elements.
     */
    int nElemConstraints;

    //!  This is the element number for the  charge neutrality 
    //!  condition of the phase
    /*!
     *  If it has one.  If it does not have a charge neutrality 
     * constraint, then this value is equal to -1    
     */
    int ChargeNeutralityElement;

    //! vector of strings containing the element names
    /*!
     * Length =  nElemConstraints
     */
    std::vector<std::string> ElName;

    //! boolean indicating whether element  constraint is active
    //! for the current  problem
    std::vector<int> ElActive;

    //! Type of the element
    /*!
     * m_elType[j] = type of the element
     *             0  VCS_ELEM_TYPE_ABSPOS Normal element that is positive
     *                                     or zero in all species.
     *             1  VCS_ELEM_TPYE_ELECTRONCHARGE element dof that corresponds
     *                                        to the charge DOF.
     *             2  VCS_ELEM_TYPE_OTHERCONSTRAINT Other constraint which may
     *                              mean that a species has neg 0 or pos value 
     *                              of that constraint (other than charge)
     */
    std::vector<int> m_elType;

    //! Formula Matrix for the phase
    /*!
     *  FormulaMatrix[j][kspec]
     *            = Formula Matrix for the species
     *              Number of elements, j,
     *              in the kspec  species  
     */
    DoubleStarStar FormulaMatrix;  

    //! Type of the species unknown
    /*!
     *  SpeciesUnknownType[k] = type of species 
     *            Normal -> VCS_SPECIES_TYPE_MOLUNK
     *                 ( unknown is the mole number in the phase)
     *            metal electron -> VCS_SPECIES_INTERFACIALVOLTAGE
     *                 ( unknown is the interfacial voltage (volts) 
     */
    std::vector<int> SpeciesUnknownType;

    //! Index of the element number in the global list of elements
    //!  storred in VCS_PROB or VCS_SOLVE       
    std::vector<int> ElGlobalIndex;

    //! Number of species in the phase       
    int NVolSpecies;

    //! String name for the phase
    std::string PhaseName;

    //!  Total moles of inert in the phase 
    double TMolesInert;

    //! molar volume of the inert species if present
    /*!
     *  units m**3 / kg
     */
    double m_molarVolInert;

    //! Convention for the activity formulation
    /*!
     *  0 = molar based activities (default)
     *  1 = Molality based activities
     *          mu = mu_0 + ln a_molality
     *          standard state is based on unity molality
     */
    int ActivityConvention;

    //! Boolean indicating whether the phase is an ideal solution
    //! and therefore it's molar-based activity coefficients are
    //! uniformly equal to one.
    bool m_isIdealSoln;

    //! Current state of existence:
    /*!
     *      0 : Doesn't exist currently
     *      1 : Does    exist currently
     *      2 : Always exists because it contains
     *          inerts which can't exist in any other
     *          phase  
     */
    int Existence; 

    //! Index of the species which is special in 
    //! with respect to the thermo treatment. 
    /*!
     * For water models this index will point to 
     * the index for water.
     *   defaults to 0
     */
    int IndexSpecialSpecies;

    //! Integer representing the activity coefficient model
    /*!
     *  The known models are listed at the top of this page
     */
    int Activity_Coeff_Model;

    //! General pointer for hanging stuff off of
    /*!
     *  Currently, not implemented very well
     */
    void *Activity_Coeff_Params;

    //! Index into the species vectors
    /*!
     *  Maps the phase species number into the global species number. 
     *  Note, as part of the vcs algorithm, the order of the species 
     *  vector is changed during the algorithm      
     */
    std::vector<int> IndSpecies;
 
    //!  Boolean indicating whether IndSpecies is contiguous 
    bool IndSpeciesContig;  

    //! Vector of Species structures for the species belonging to this phase
    /*!
     * The index into this vector is the species index within the phase.
     */
    std::vector<vcs_SpeciesProperties *> ListSpeciesPtr;

    //! Units for the chemical potential data, pressure data, volume,
    //! and species amounts
    /*!
     *  All internally storred quantities will have these units. Also, printed
     *  quantitities will display in these units. Input quantities are expected
     *  in these units.
     *
     *                           Chem_Pot                 Pres      vol   moles
     * ----------------------------------------------------------------------
     * -1  VCS_UNITS_KCALMOL  = kcal/gmol                 Pa     m**3   kmol
     *  0  VCS_UNITS_UNITLESS = MU / RT -> no units       Pa     m**3   kmol
     *  1  VCS_UNITS_KJMOL    = kJ / gmol                 Pa     m**3   kmol
     *  2  VCS_UNITS_KELVIN   = KELVIN -> MU / R          Pa     m**3   kmol
     *  3  VCS_UNITS_MKS      = Joules / Kmol (Cantera)   Pa     m**3   kmol
     * ----------------------------------------------------------------------
     *
     *  see vcs_defs.h for more information.
     *
     *  Currently, this value should be the same as the owning VCS_PROB or
     *  VCS_SOLVE object. There is no code for handling anything else atm.
     */
    int m_VCS_UnitsFormat;

  private:
    //!  If this is true, then calculations are actually performed within
    //!  Cantera
    bool m_useCanteraCalls;
    /**
     *  If we are using Cantera, this is the
     *	pointer to the ThermoPhase object. If not, this is null. 
     */
    Cantera::ThermoPhase *TP_ptr;
 
    //!  Total mols in the phase 
    /*!
     *      units are kmol
     */
    double v_totalMoles;

    //! Vector of the current mole fractions for species 
    //! in the phase
    std::vector<double> Xmol;

  public:
  
    //! If the potential is a solution variable in VCS, it acts as a species.
    //!  This is the species index in the phase for the potential
    int m_phiVarIndex;
 
    //! Total Volume of the phase 
    /*!
     *  units are m**3
     */
    mutable double m_totalVol;
  
    //! Vector of calculated SS0 chemical potentials for the
    //! current Temperature. 
    /*!
     * Note, This is the chemical potential derived strictly from the polynomial
     * in temperature. Pressure effects have to be added in to
     * get to the standard state.
     *
     * Units -> depends on VCS_UnitsFormat variable
     *             Cantera -> J/kmol
     */
    mutable std::vector<double> SS0ChemicalPotential;

  private:
    //! Vector of calculated Star chemical potentials for the
    //! current Temperature and pressure.
    /*!
     * Note, This is the chemical potential at unit activity. Thus, we can call
     * it the standard state chemical potential as well.
     *
     * Units -> depends on VCS_UnitsFormat variable
     *             Cantera -> J/kmol
     */
    mutable std::vector<double> StarChemicalPotential;

  public:
    //! Vector of the Star molar Volumes of the species.
    /*!
     * units  m3 / kmol
     */
    mutable std::vector<double> StarMolarVol;

    //! Vector of the Partial molar Volumes of the species.
    /*!
     * units  m3 / kmol
     */
    mutable std::vector<double> PartialMolarVol;

    //! Vector of calculated activity coefficients for the current state
    /*!
     *  Whether or not this vector is current is determined by
     *  the bool m_UpToDate_AC.
     */
    mutable std::vector<double> ActCoeff;

    //! Vector of the derivatives of the ln activity coefficient wrt to the
    //! current mole number
    /*!
     * dLnActCoeffdMolNumber[j][k];
     *      j = id of the species mole number
     *      k = id of the species activity coefficient
     */
    mutable DoubleStarStar dLnActCoeffdMolNumber;


    //! Status 
    /*!
     *  valid values are 
     *       VCS_STATECALC_OLD
     *       VCS_STATECALC_NEW
     */
    int  m_vcsStateStatus;

  private:

    //! Value of the potential for the phase (Volts)
    double m_phi;

  public:
    //! Boolean indicating whether the object has an uptodate mole number vector
    //! and potential with respect to the current vcs state calc status
    bool m_UpToDate;

  private:
    //! Boolean indicating whether activity coefficients are uptodate.
    /*!
     * Activity coefficients and volume calculations are lagged. They are only
     * called when they are needed (and when the state has changed so that they
     * need to be recalculated).
     */
    mutable bool m_UpToDate_AC;

    //! Boolean indicating whether Star volumes are uptodate.
    /*!
     * Activity coefficients and volume calculations are lagged. They are only
     * called when they are needed (and when the state has changed so that they
     * need to be recalculated).
     *  Star volumes are sensitive to temperature and pressure
     */
    mutable bool m_UpToDate_VolStar;

    //! Boolean indicating whether partial molar volumes are uptodate.
    /*!
     * Activity coefficients and volume calculations are lagged. They are only
     * called when they are needed (and when the state has changed so that they
     * need to be recalculated).
     *  partial molar volumes are sensitive to everything
     */
    mutable bool m_UpToDate_VolPM;

    //! Boolean indicating whether GStar is uptodate.
    /*!
     * GStar is sensitive to the temperature and the pressure, only
     */
    mutable bool m_UpToDate_GStar;

    //! Current value of the temperature for this object, and underlying objects
    double Temp;

    //! Current value of the pressure for this object, and underlying objects
    double Pres;

  public:
    //! Reference pressure for the phase
    double RefPres;

  };

  //! Return a string representing the equation of state
  /*!
   *  @param EOSType : integer value of the equation of state
   *
   * @return returns a string representing the EOS
   */
  std::string string16_EOSType(int EOSType);

}

#endif
