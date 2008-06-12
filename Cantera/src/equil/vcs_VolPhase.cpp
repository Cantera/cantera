/**
 * @file vcs_VolPhase.cpp
 */

/* $Id$ */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#include "vcs_VolPhase.h"
#include "vcs_internal.h"
#include "vcs_SpeciesProperties.h"
#include "vcs_species_thermo.h"
#include "vcs_solve.h"

#include "ThermoPhase.h"
#include "mix_defs.h"

#include <cstdio>
#include <cstdlib>

namespace VCSnonideal {

  /*
   * 
   *  vcs_VolPhase():
   *
   *    Constructor for the VolPhase object.
   */
  vcs_VolPhase::vcs_VolPhase(VCS_SOLVE * owningSolverObject) :
    m_owningSolverObject(0),
    VP_ID(-1),
    Domain_ID(-1),
    SingleSpecies(true),
    m_gasPhase(false),
    EqnState(VCS_EOS_CONSTANT),
    nElemConstraints(0),
    ChargeNeutralityElement(-1),
    ElGlobalIndex(0),
    NVolSpecies(0),
    TMolesInert(0.0),
    m_molarVolInert(1000.),
    ActivityConvention(0),
    m_isIdealSoln(false),
    Existence(0),
    IndexSpecialSpecies(-1),
    Activity_Coeff_Model(VCS_AC_CONSTANT),
    Activity_Coeff_Params(0),
    IndSpecies(0),
    IndSpeciesContig(true),
    m_VCS_UnitsFormat(VCS_UNITS_MKS),
    m_useCanteraCalls(false),
    TP_ptr(0),
    TMoles(0.0),
    m_totalVol(0.0),
    m_vcsStateStatus(VCS_STATECALC_OLD),
    m_phi(0.0),
    m_UpToDate(false),
    m_UpToDate_AC(false),
    m_UpToDate_VolStar(false),
    m_UpToDate_VolPM(false),
    m_UpToDate_GStar(false),
    Temp(273.15),
    Pres(1.01325E5),
    RefPres(1.01325E5)
  {
    m_owningSolverObject = owningSolverObject;
  }
  /************************************************************************************/

  /*
   * 
   *  ~vcs_VolPhase():
   *
   *   Destructor for the VolPhase object.
   */
  vcs_VolPhase::~vcs_VolPhase() {
    for (int k = 0; k < NVolSpecies; k++) {
      vcs_SpeciesProperties *sp = ListSpeciesPtr[k];
      delete sp;
      sp = 0;
    }
  }
  /************************************************************************************/

  /*
   * 
   *  Copy Constructor():
   *
   *  Objects that are owned by this object are deep copied here, except
   *  for the ThermoPhase object.
   *  The assignment operator does most of the work.
   */
  vcs_VolPhase::vcs_VolPhase(const vcs_VolPhase& b) :
    m_owningSolverObject(b.m_owningSolverObject),
    VP_ID(b.VP_ID),
    Domain_ID(b.Domain_ID),
    SingleSpecies(b.SingleSpecies),
    m_gasPhase(b.m_gasPhase),
    EqnState(b.EqnState),
    nElemConstraints(b.nElemConstraints),
    ChargeNeutralityElement(b.ChargeNeutralityElement),
    NVolSpecies(b.NVolSpecies),
    TMolesInert(b.TMolesInert),
    ActivityConvention(b.ActivityConvention),
    m_isIdealSoln(b.m_isIdealSoln),
    Existence(b.Existence),
    IndexSpecialSpecies(b.IndexSpecialSpecies),
    Activity_Coeff_Model(b.Activity_Coeff_Model),
    Activity_Coeff_Params(b.Activity_Coeff_Params),
    IndSpeciesContig(b.IndSpeciesContig),
    m_VCS_UnitsFormat(b.m_VCS_UnitsFormat),
    m_useCanteraCalls(b.m_useCanteraCalls),
    TP_ptr(b.TP_ptr),
    TMoles(b.TMoles),
    m_phiVarIndex(-1),
    m_totalVol(b.m_totalVol),
    m_vcsStateStatus(VCS_STATECALC_OLD),
    m_phi(b.m_phi),
    m_UpToDate(false),
    m_UpToDate_AC(false),
    m_UpToDate_VolStar(false),
    m_UpToDate_VolPM(false),
    m_UpToDate_GStar(false),
    Temp(b.Temp),
    Pres(b.Pres)
  {
    /*
     * Call the Assignment operator to do the heavy
     * lifting.
     */
    *this = b;
  }
  /*****************************************************************************/
  
  /*
   * Assignment operator()
   *
   *   (note, this is used, so keep it current!)
   */
  vcs_VolPhase& vcs_VolPhase::operator=(const vcs_VolPhase& b)
  {
    int k;
    if (&b != this) {
      int old_num = NVolSpecies;

      m_owningSolverObject = b.m_owningSolverObject;
      VP_ID               = b.VP_ID;
      Domain_ID           = b.Domain_ID;
      SingleSpecies       = b.SingleSpecies;
      m_gasPhase            = b.m_gasPhase;
      EqnState            = b.EqnState;
 
      NVolSpecies         = b.NVolSpecies;
      nElemConstraints    = b.nElemConstraints;
      ChargeNeutralityElement = b.ChargeNeutralityElement;


      ElName.resize(b.nElemConstraints);
      for (int e = 0; e < b.nElemConstraints; e++) {
	ElName[e] = b.ElName[e];
      }
 
      ElActive = b.ElActive;
      m_elType = b.m_elType;
  
      FormulaMatrix.resize(nElemConstraints, NVolSpecies, 0.0);
      for (int e = 0; e < nElemConstraints; e++) {
	for (int k = 0; k < NVolSpecies; k++) {
	  FormulaMatrix[e][k] = b.FormulaMatrix[e][k];
	}
      }

      SpeciesUnknownType = b.SpeciesUnknownType;
      ElGlobalIndex = b.ElGlobalIndex;
      NVolSpecies         = b.NVolSpecies;
      PhaseName           = b.PhaseName;
      TMolesInert         = b.TMolesInert;
      ActivityConvention  = b.ActivityConvention;
      m_isIdealSoln       = b.m_isIdealSoln;
      Existence           = b.Existence;
      IndexSpecialSpecies = b.IndexSpecialSpecies;
      Activity_Coeff_Model = b.Activity_Coeff_Model;

      /*
       * Do a shallow copy because we haven' figured this out.
       */
      Activity_Coeff_Params = b.Activity_Coeff_Params;
      IndSpecies = b.IndSpecies;
      IndSpeciesContig = b.IndSpeciesContig;

      for (k = 0; k < old_num; k++) {
	if ( ListSpeciesPtr[k]) {
	  delete  ListSpeciesPtr[k];
	  ListSpeciesPtr[k] = 0;
	}
      }
      ListSpeciesPtr.resize(NVolSpecies, 0);
      for (k = 0; k < NVolSpecies; k++) {
	ListSpeciesPtr[k] = 
	  new vcs_SpeciesProperties(*(b.ListSpeciesPtr[k]));
      }
    
      m_VCS_UnitsFormat   = b.m_VCS_UnitsFormat;
      m_useCanteraCalls   = b.m_useCanteraCalls;
      /*
       * Do a shallow copy of the ThermoPhase object pointer.
       * We don't duplicate the object.
       *  Um, there is no reason we couldn't do a 
       *  duplicateMyselfAsThermoPhase() call here. This will
       *  have to be looked into.
       */
      TP_ptr              = b.TP_ptr;
      TMoles              = b.TMoles;
 
      Xmol = b.Xmol;

      m_phi               = b.m_phi;
      m_phiVarIndex       = b.m_phiVarIndex;
 
      SS0ChemicalPotential = b.SS0ChemicalPotential;
      StarChemicalPotential = b.StarChemicalPotential;

      StarMolarVol = b.StarMolarVol;
      PartialMolarVol = b.PartialMolarVol; 
      ActCoeff = b.ActCoeff;

      dLnActCoeffdMolNumber = b.dLnActCoeffdMolNumber;

      m_UpToDate            = false;
      m_vcsStateStatus      = b.m_vcsStateStatus;
      m_UpToDate_AC         = false;
      m_UpToDate_VolStar    = false;
      m_UpToDate_VolPM      = false;
      m_UpToDate_GStar      = false;
      Temp                = b.Temp;
      Pres                = b.Pres;
      setState_TP(Temp, Pres);
      _updateMoleFractionDependencies();
    }
    return *this;
  }
  /************************************************************************************/

  void vcs_VolPhase::resize(int phaseNum, int nspecies, const char *phaseName,
			    double molesInert) {
    if (nspecies <= 0) {
      plogf("nspecies Error\n");
      std::exit(-1);
    }
    if (phaseNum < 0) {
      plogf("phaseNum should be greater than 0\n");
      std::exit(-1);
    }

    TMolesInert = molesInert;
    if (TMolesInert > 0.0) {
      Existence = 2;
    } 

    m_phi = 0.0;
    m_phiVarIndex = -1;

    if (phaseNum == VP_ID) {
      if (strcmp(PhaseName.c_str(), phaseName)) {
	plogf("Strings are different: %s %s :unknown situation\n",
	      PhaseName.c_str(), phaseName);
	std::exit(-1);
      }
    } else {
      VP_ID = phaseNum;
      if (!phaseName) {
	char itmp[40];
	sprintf(itmp, "Phase_%d", VP_ID);
	PhaseName = itmp;
      } else {
	PhaseName = phaseName;
      }
    }
    if (nspecies > 1) {
      SingleSpecies = false;
    } else {
      SingleSpecies = true;
    }

    if (NVolSpecies == nspecies) {
      return;
    }
 
    NVolSpecies = nspecies;
    if (nspecies > 1) {
      SingleSpecies = false;
    }

    IndSpecies.resize(nspecies,-1);

    if ((int) ListSpeciesPtr.size() >=  NVolSpecies) {
      for (int i = 0; i < NVolSpecies; i++) {
	if (ListSpeciesPtr[i]) {
	  delete ListSpeciesPtr[i]; 
	  ListSpeciesPtr[i] = 0;
	}
      }
    }
    ListSpeciesPtr.resize(nspecies, 0);
    for (int i = 0; i < nspecies; i++) {
      ListSpeciesPtr[i] = new vcs_SpeciesProperties(phaseNum, i, this);
    }

    Xmol.resize(nspecies, 0.0);
    for (int i = 0; i < nspecies; i++) {
      Xmol[i] = 1.0/nspecies;
    }

    SS0ChemicalPotential.resize(nspecies, -1.0);
    StarChemicalPotential.resize(nspecies, -1.0);
    StarMolarVol.resize(nspecies, -1.0);
    PartialMolarVol.resize(nspecies, -1.0);
    ActCoeff.resize(nspecies, 1.0);
    dLnActCoeffdMolNumber.resize(nspecies, nspecies, 0.0);
 

    SpeciesUnknownType.resize(nspecies, VCS_SPECIES_TYPE_MOLNUM);
    m_UpToDate            = false;
    m_vcsStateStatus      = VCS_STATECALC_OLD;
    m_UpToDate_AC         = false;
    m_UpToDate_VolStar    = false;
    m_UpToDate_VolPM      = false;
    m_UpToDate_GStar      = false;
  }
  /*******************************************************************************/

  //! Evaluate activity coefficients
  /*!
   *   We carry out a calculation whenever UpTODate_AC is false. Specifically
   *   whenever a phase goes zero, we do not carry out calculations on it.
   */
  void vcs_VolPhase::evaluateActCoeff() const {
    if (m_UpToDate_AC == true) return;
    if (m_isIdealSoln) {
      m_UpToDate_AC = true;
      return;
    }
    if (m_useCanteraCalls) {
      TP_ptr->getActivityCoefficients(VCS_DATA_PTR(ActCoeff));
    } else {
      switch (Activity_Coeff_Model) {
      case VCS_AC_CONSTANT:
	/*
	 * Don't need to do anything since ActCoeff[] is initialized to
	 * the value of one, and never changed for this model.
	 */
	break;
      default:
	plogf("%sERROR: unknown model\n");
	std::exit(-1);
      }
    }
    m_UpToDate_AC = true;
  }
  /********************************************************************************/

  /*
   *
   * Evaluate one activity coefficients.
   *
   *   return one activity coefficient. Have to recalculate them all to get
   *   one.
   */
  double vcs_VolPhase::AC_calc_one(int kspec) const {
    evaluateActCoeff();
    return(ActCoeff[kspec]);
  }
  /************************************************************************************/

  // Gibbs free energy calculation at a temperature for the reference state
  // of each species
  /*
   *  @param TKelvin temperature
   */
  void vcs_VolPhase::G0_calc(double tkelvin) {
    bool lsame = false;
    if (Temp == tkelvin) {
      lsame = true;
    }

    bool doit = !lsame;
    setState_TP(tkelvin, Pres);
    if (SS0ChemicalPotential[0] == -1) doit = true;
    if (doit) {
      if (m_useCanteraCalls) {
	TP_ptr->getGibbs_ref(VCS_DATA_PTR(SS0ChemicalPotential));
      } else {
	double R = vcsUtil_gasConstant(m_VCS_UnitsFormat);
	for (int k = 0; k < NVolSpecies; k++) {
	  int kglob = IndSpecies[k];
	  vcs_SpeciesProperties *sProp = ListSpeciesPtr[k];
	  VCS_SPECIES_THERMO *sTherm = sProp->SpeciesThermo;
	  SS0ChemicalPotential[k] =
	    R * (sTherm->G0_R_calc(kglob, tkelvin));
	}
      }
    }
  }
  /***********************************************************************/

  // Gibbs free energy calculation at a temperature for the reference state
  // of a species, return a value for one species
  /*
   *  @param kspec   species index
   *  @param TKelvin temperature
   *
   *  @return return value of the gibbs free energy
   */
  double vcs_VolPhase::G0_calc_one(int kspec, double tkelvin) {
    G0_calc(tkelvin);
    return SS0ChemicalPotential[kspec];
  }
  /***********************************************************************/

  // Gibbs free energy calculation for standard states
  /*
   * Calculate the Gibbs free energies for the standard states
   * The results are held internally within the object.
   *
   * @param TKelvin Current temperature
   * @param pres    Current pressure (pascal)
   */
  void vcs_VolPhase::GStar_calc() const {
    if (!m_UpToDate_GStar) {
      if (m_useCanteraCalls) {
	TP_ptr->getStandardChemPotentials(VCS_DATA_PTR(StarChemicalPotential));
      } else {
	double R = vcsUtil_gasConstant(m_VCS_UnitsFormat);
	for (int k = 0; k < NVolSpecies; k++) {
	  int kglob = IndSpecies[k];
	  vcs_SpeciesProperties *sProp = ListSpeciesPtr[k];
	  VCS_SPECIES_THERMO *sTherm = sProp->SpeciesThermo;
	  StarChemicalPotential[k] =
	    R * (sTherm->GStar_R_calc(kglob, Temp, Pres));
	}
      }
      m_UpToDate_GStar = true;
    }
  }
  /***********************************************************************/

  // Gibbs free energy calculation for standard state of one species
  /*
   * Calculate the Gibbs free energies for the standard state
   * of the kth species.
   * The results are held internally within the object.
   * The kth species standard state G is returned
   *
   * @param kspec   Species number (within the phase)
   *
   * @return Gstar[kspec] returns the gibbs free energy for the
   *         standard state of the kspec species.
   */
  double vcs_VolPhase::GStar_calc_one(int kspec) {
    if (!m_UpToDate_GStar) {
      GStar_calc();
    }
    return StarChemicalPotential[kspec];
  }
  /***********************************************************************/

  // Set the mole fractions from a conventional mole fraction vector
  /*
   *
   * @param xmol Value of the mole fractions for the species
   *             in the phase. These are contiguous. 
   */
  void vcs_VolPhase::setMoleFractions(const double * const xmol) {
    double sum = -1.0;
    for (int k = 0; k < NVolSpecies; k++) {
      Xmol[k] = xmol[k];
      sum+= xmol[k];
    }
    if (std::fabs(sum) > 1.0E-13) {
      for (int k = 0; k < NVolSpecies; k++) {
	Xmol[k] /= sum;
      }
    }
    _updateMoleFractionDependencies();
    m_UpToDate = false;
    m_vcsStateStatus = VCS_STATECALC_TMP;
  }
  /***********************************************************************/

  // Updates the mole fractions in subobjects
  /*
   *  Whenever the mole fractions change, this routine
   *  should be called.
   */
  void vcs_VolPhase::_updateMoleFractionDependencies() {
    if (m_useCanteraCalls) {
      if (TP_ptr) {
	TP_ptr->setState_PX(Pres, VCS_DATA_PTR(Xmol));
      }
    }
    if (!m_isIdealSoln) {
      m_UpToDate_AC = false;
      m_UpToDate_VolPM = false;
    }
  }
  /************************************************************************/

  // Return a const reference to the mole fraction vector in the phase
  const std::vector<double> & vcs_VolPhase::moleFractions() const {
    return Xmol;
  }
  /***********************************************************************/

  // Set the moles within the phase
  /*
   *  This function takes as input the mole numbers in vcs format, and
   *  then updates this object with their values. This is essentially
   *  a gather routine.
   *
   *  
   *  @param molesSpeciesVCS  array of mole numbers. Note, the indecises 
   *            for species in 
   *            this array may not be contiguous. IndSpecies[] is needed
   *            to gather the species into the local contiguous vector
   *            format. 
   */
  void vcs_VolPhase::setMolesFromVCS(const int stateCalc, 
				     const double * molesSpeciesVCS) {
    int kglob;
    double tmp;
    TMoles = TMolesInert;

    if (molesSpeciesVCS == 0) {
#ifdef DEBUG_MODE
      if (m_owningSolverObject == 0) {
	printf("shouldn't be here\n");
	std::exit(-1);
      }
#endif
      if (stateCalc == VCS_STATECALC_OLD) {
	molesSpeciesVCS = VCS_DATA_PTR(m_owningSolverObject->m_molNumSpecies_old);
      } else if (stateCalc == VCS_STATECALC_NEW) {
	molesSpeciesVCS = VCS_DATA_PTR(m_owningSolverObject->m_molNumSpecies_new);
      }
#ifdef DEBUG_MODE
      else {
	printf("shouldn't be here\n");
	std::exit(-1);
      }
#endif
    }
#ifdef DEBUG_MODE
    else {
      if (m_owningSolverObject) {
        if (stateCalc == VCS_STATECALC_OLD) {
  	  if (molesSpeciesVCS != VCS_DATA_PTR(m_owningSolverObject->m_molNumSpecies_old)) {
	    printf("shouldn't be here\n");
	    std::exit(-1);
          }
        } else if (stateCalc == VCS_STATECALC_NEW) {
          if (molesSpeciesVCS != VCS_DATA_PTR(m_owningSolverObject->m_molNumSpecies_new)) {
	    printf("shouldn't be here\n");
	    std::exit(-1);
          }
        }
      }
    }
#endif

    for (int k = 0; k < NVolSpecies; k++) {
      if (SpeciesUnknownType[k] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	kglob = IndSpecies[k];
	tmp = MAX(0.0, molesSpeciesVCS[kglob]);
	Xmol[k] = tmp;
	TMoles += tmp;
      }
    }
    if (TMoles > 0.0) {
      for (int k = 0; k < NVolSpecies; k++) {
	Xmol[k] /= TMoles;
      }
      Existence = 1;
    } else {
      // This is where we will start to store a better approximation 
      // for the mole fractions, when the phase doesn't exist.
      // This is currently unimplemented.
      for (int k = 0; k < NVolSpecies; k++) {
	Xmol[k] = 1.0 / NVolSpecies;
      }
      Existence = 0;
    }
    /*
     * Update the electric potential if it is a solution variable
     * in the equation system
     */
    if (m_phiVarIndex >= 0) {
      kglob = IndSpecies[m_phiVarIndex];
      if (NVolSpecies == 1) {
	Xmol[m_phiVarIndex] = 1.0;
      } else {
	Xmol[m_phiVarIndex] = 0.0;
      }
      double phi = molesSpeciesVCS[kglob];
      setElectricPotential(phi);
      if (NVolSpecies == 1) {
	Existence = 1;
      }
    }
    _updateMoleFractionDependencies();
    if (TMolesInert > 0.0) {
      Existence = 2;
    }
    /*
     * Set flags indicating we are up to date with the VCS state vector.
     */
    m_UpToDate = true;
    m_vcsStateStatus = stateCalc; 
 
  }
  /***********************************************************************/

  // Set the moles within the phase
  /*
   *  This function takes as input the mole numbers in vcs format, and
   *  then updates this object with their values. This is essentially
   *  a gather routine.
   *
   *  
   *  @param molesSpeciesVCS  array of mole numbers. Note, the indecises for species in 
   *            this array may not be contiguous. IndSpecies[] is needed
   *            to gather the species into the local contiguous vector
   *            format. 
   */
  void vcs_VolPhase::setMolesFromVCSCheck(const int stateCalc,
					  const double * molesSpeciesVCS, 
					  const double * const TPhMoles) {
    setMolesFromVCS(stateCalc, molesSpeciesVCS);
    /*
     * Check for consistency with TPhMoles[]
     */
    double Tcheck = TPhMoles[VP_ID];
    if (Tcheck != TMoles) {
      if (vcs_doubleEqual(Tcheck, TMoles)) {
	Tcheck = TMoles;
      } else {
	plogf("vcs_VolPhase::setMolesFromVCSCheck: "
	      "We have a consistency problem: %21.16g %21.16g\n",
	      Tcheck, TMoles);
	std::exit(-1);
      }
    }
  }
  /***********************************************************************/

  // Update the moles within the phase, if necessary
  /*
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
  void vcs_VolPhase::updateFromVCS_MoleNumbers(const int stateCalc) {
    if (!m_UpToDate || (stateCalc != m_vcsStateStatus)) {
      if (stateCalc == VCS_STATECALC_OLD || stateCalc == VCS_STATECALC_NEW) {
	if (m_owningSolverObject) {
	  setMolesFromVCS(stateCalc);
	}
      }
    }
  }
  /***********************************************************************/

  // Fill in an activity coefficients vector within a VCS_SOLVE object
  /*
   *  This routine will calculate the activity coefficients for the
   *  current phase, and fill in the corresponding entries in the
   *  VCS activity coefficients vector.
   *  
   * @param AC  vector of activity coefficients for all of the species
   *            in all of the phases in a VCS problem. Only the
   *            entries for the current phase are filled in.
   */
  void vcs_VolPhase::sendToVCS_ActCoeff(double * const AC) const {
    if (!m_UpToDate_AC) {
      evaluateActCoeff();
    }
    int kglob;
    for (int k = 0; k < NVolSpecies; k++) {
      kglob = IndSpecies[k];
      AC[kglob] = ActCoeff[k];
    }
  }
  /***********************************************************************/

  // Fill in the partial molar volume vector for VCS
  /*
   *  This routine will calculate the partial molar volumes for the
   *  current phase (if needed), and fill in the corresponding entries in the
   *  VCS partial molar volumes vector.
   *  
   * @param VolPM  vector of partial molar volumes for all of the species
   *            in all of the phases in a VCS problem. Only the
   *            entries for the current phase are filled in.
   */
  double vcs_VolPhase::sendToVCS_VolPM(double * const VolPM) const {  
    if (!m_UpToDate_VolPM) {
      (void) VolPM_calc();
    }
    int kglob;
    for (int k = 0; k < NVolSpecies; k++) {
      kglob = IndSpecies[k];
      VolPM[kglob] = PartialMolarVol[k];
    }
    return m_totalVol;
  }
  /***********************************************************************/

  // Fill in the partial molar volume vector for VCS
  /*
   *  This routine will calculate the partial molar volumes for the
   *  current phase (if needed), and fill in the corresponding entries in the
   *  VCS partial molar volumes vector.
   *  
   * @param VolPM  vector of partial molar volumes for all of the species
   *            in all of the phases in a VCS problem. Only the
   *            entries for the current phase are filled in.
   */
  void vcs_VolPhase::sendToVCS_GStar(double * const gstar){  
    if (!m_UpToDate_GStar) {
      GStar_calc();
    }
    int kglob;
    for (int k = 0; k < NVolSpecies; k++) {
      kglob = IndSpecies[k];
      gstar[kglob] = StarChemicalPotential[k];
    }
  }
 /***********************************************************************/


  void vcs_VolPhase::setElectricPotential(double phi) {
    m_phi = phi;
    if (m_useCanteraCalls) {
      TP_ptr->setElectricPotential(m_phi);
    }
    // We have changed the state variable. Set uptodate flags to false
    m_UpToDate_AC = false;
    m_UpToDate_VolStar = false;
    m_UpToDate_VolPM = false;
    m_UpToDate_GStar = false;
  }
  /***********************************************************************/

  double vcs_VolPhase::electricPotential() const {
    return m_phi;
  }
  /***********************************************************************/

  // Sets the temperature and pressure in this object and
  //  underlying objects
  /*
   *  Sets the temperature and pressure in this object and
   *  underlying objects. The underlying objects refers to the
   *  Cantera's ThermoPhase object for this phase.
   *
   *  @param temperature_Kelvin    (Kelvin)
   *  @param pressure_PA  Pressure (MKS units - Pascal)
   */
  void vcs_VolPhase::setState_TP(double temp, double pres)
  {
    if (Temp == temp) {
      if (Pres == pres) {
	return;
      }
    }
    if (m_useCanteraCalls) {
      TP_ptr->setElectricPotential(m_phi);
      TP_ptr->setState_TP(temp, pres);
    }
    Temp = temp;
    Pres = pres;
    m_UpToDate_AC      = false;
    m_UpToDate_VolStar = false;
    m_UpToDate_VolPM   = false;
    m_UpToDate_GStar   = false;
  }
  /***********************************************************************/

  // Molar volume calculation for standard states
  /*
   * Calculate the molar volume for the standard states
   * The results are held internally within the object.
   *
   * @param TKelvin Current temperature
   * @param pres    Current pressure (pascal)
   *
   *  Calculations are in m**3/kmol
   */
  void vcs_VolPhase::VolStar_calc() const {
    if (!m_UpToDate_VolStar) {     
      if (m_useCanteraCalls) {
	TP_ptr->getStandardVolumes(VCS_DATA_PTR(StarMolarVol));
      } else {
	for (int k = 0; k < NVolSpecies; k++) {
	  int kglob = IndSpecies[k];
	  vcs_SpeciesProperties *sProp = ListSpeciesPtr[k];
	  VCS_SPECIES_THERMO *sTherm = sProp->SpeciesThermo;
	  StarMolarVol[k] =
	    (sTherm->VolStar_calc(kglob, Temp, Pres));
	}
      }
      m_UpToDate_VolStar = true;
    }
  }
  /***********************************************************************/


  // Molar volume calculation for standard state of one species
  /*
   * Calculate the molar volume for the standard states
   * The results are held internally within the object.
   * Return the molar volume for one species
   *
   * @param kspec   Species number (within the phase)
   * @param TKelvin Current temperature
   * @param pres    Current pressure (pascal)
   *
   * @return molar volume of the kspec species's standard
   *         state
   */
  double vcs_VolPhase::VolStar_calc_one(int kspec, double tkelvin, 
					double pres) {
    setState_TP(tkelvin, pres);
    if (!m_UpToDate_VolStar) { 
      VolStar_calc();
    }
    return StarMolarVol[kspec];
  }
  /****************************************************************************/

  /*
   *
   * VolPM_calc
   */
  double vcs_VolPhase::VolPM_calc() const {
    int k, kglob;
    if (!m_UpToDate_VolPM) {
      if (m_useCanteraCalls) {
	TP_ptr->getPartialMolarVolumes(VCS_DATA_PTR(PartialMolarVol));
      } else {
	for (k = 0; k < NVolSpecies; k++) {
	  kglob = IndSpecies[k];
	  vcs_SpeciesProperties *sProp = ListSpeciesPtr[k];
	  VCS_SPECIES_THERMO *sTherm = sProp->SpeciesThermo;
	  StarMolarVol[k] = (sTherm->VolStar_calc(kglob, Temp, Pres));
	}
	for (k = 0; k < NVolSpecies; k++) {
	  PartialMolarVol[k] = StarMolarVol[k];
	}
      }

      m_totalVol = 0.0;
      for (k = 0; k < NVolSpecies; k++) {
	m_totalVol += PartialMolarVol[k] * Xmol[k];
      }
      m_totalVol *= TMoles;

      if (TMolesInert > 0.0) {
	if (m_gasPhase) {
	  double volI = TMolesInert * 8314.47215 * Temp / Pres;
	  m_totalVol += volI;
	} else {
	  printf("unknown situation\n");
	  std::exit(-1);
	}
      }
    }
    m_UpToDate_VolPM = true;
    return m_totalVol;
  }
  /************************************************************************************/

  /*
   * updateLnActCoeffJac():
   *
   */
  void vcs_VolPhase::updateLnActCoeffJac() {
    int k, j;
    double deltaMoles_j = 0.0;
  

    /*
     * Evaluate the current base activity coefficients.
     */  
    evaluateActCoeff();

    // Make copies of ActCoeff and Xmol for use in taking differences
    std::vector<double> ActCoeff_Base(ActCoeff);
    std::vector<double> Xmol_Base(Xmol);
    double TMoles_base = TMoles;

    /*
     *  Loop over the columns species to be deltad
     */
    for (j = 0; j < NVolSpecies; j++) {
      /*
       * Calculate a value for the delta moles of species j
       * -> NOte Xmol[] and Tmoles are always positive or zero
       *    quantities.
       */
      double moles_j_base = TMoles * Xmol_Base[j];
      deltaMoles_j = 1.0E-7 * moles_j_base + 1.0E-20 * TMoles + 1.0E-150;
      /*
       * Now, update the total moles in the phase and all of the
       * mole fractions based on this.
       */
      TMoles = TMoles_base + deltaMoles_j;      
      for (k = 0; k < NVolSpecies; k++) {
	Xmol[k] = Xmol_Base[k] * TMoles_base / TMoles;
      }
      Xmol[j] = (moles_j_base + deltaMoles_j) / TMoles;
 
      /*
       * Go get new values for the activity coefficients.
       * -> Note this calls setState_PX();
       */
      _updateMoleFractionDependencies();
      evaluateActCoeff();
      /*
       * Calculate the column of the matrix
       */
      double * const lnActCoeffCol = dLnActCoeffdMolNumber[j];
      for (k = 0; k < NVolSpecies; k++) {
	lnActCoeffCol[k] = (ActCoeff[k] - ActCoeff_Base[k]) /
	  ((ActCoeff[k] + ActCoeff_Base[k]) * 0.5 * deltaMoles_j);
      }
      /*
       * Revert to the base case Xmol, TMoles
       */
      TMoles = TMoles_base;
      vcs_vdcopy(Xmol, Xmol_Base, NVolSpecies);
    }
    /*
     * Go get base values for the activity coefficients.
     * -> Note this calls setState_TPX() again;
     * -> Just wanted to make sure that cantera is in sync
     *    with VolPhase after this call.
     */
    setMoleFractions(VCS_DATA_PTR(Xmol_Base));
    _updateMoleFractionDependencies();
    evaluateActCoeff();
  }
  /************************************************************************************/

  // Downloads the ln ActCoeff jacobian into the VCS version of the
  // ln ActCoeff jacobian.
  /*
   *
   *   This is essentially a scatter operation.
   *
   *   The Jacobians are actually d( lnActCoeff) / d (MolNumber);
   *   dLnActCoeffdMolNumber[j][k]
   * 
   *      j = id of the species mole number
   *      k = id of the species activity coefficient
   */
  void vcs_VolPhase::sendToVCS_LnActCoeffJac(double * const * const LnACJac_VCS) {
    /*
     * update the Ln Act Coeff jacobian entries with respect to the
     * mole number of species in the phase -> we always assume that
     * they are out of date.
     */
    updateLnActCoeffJac();
    /*
     *  Now copy over the values
     */
    int j, k, jglob, kglob;
    for (j = 0; j < NVolSpecies; j++) {
      jglob = IndSpecies[j];
      double * const lnACJacVCS_col = LnACJac_VCS[jglob];
      const double * const lnACJac_col = dLnActCoeffdMolNumber[j];
      for (k = 0; k < NVolSpecies; k++) {
	kglob = IndSpecies[k];
	lnACJacVCS_col[kglob] = lnACJac_col[k];
      }
    }
  }
  /************************************************************************************/

  // Set the pointer for Cantera's ThermoPhase parameter
  /*
   *  When we first initialize the ThermoPhase object, we read the
   *  state of the ThermoPhase into vcs_VolPhase object.
   *
   * @param tp_ptr Pointer to the ThermoPhase object corresponding
   *               to this phase.
   */
  void vcs_VolPhase::setPtrThermoPhase(Cantera::ThermoPhase *tp_ptr) {
    TP_ptr = tp_ptr;
    if (TP_ptr) {
      m_useCanteraCalls = true;
      Temp = TP_ptr->temperature();
      Pres = TP_ptr->pressure();
      setState_TP(Temp, Pres);
      m_VCS_UnitsFormat = VCS_UNITS_MKS;
      m_phi = TP_ptr->electricPotential();
      int nsp = TP_ptr->nSpecies();
      if (nsp !=  NVolSpecies) {
	if (NVolSpecies != 0) {
	  plogf("Warning Nsp != NVolSpeces: %d %d \n", nsp, NVolSpecies);
	}
	resize(VP_ID, nsp, PhaseName.c_str());
      }
      TP_ptr->getMoleFractions(VCS_DATA_PTR(Xmol));
      _updateMoleFractionDependencies();

      /*
       *  figure out ideal solution tag
       */
      if (nsp == 1) {
	m_isIdealSoln = true;
      } else {
	int eos = TP_ptr->eosType();
	switch (eos) {
	case Cantera::cIdealGas:
	case Cantera::cIncompressible:
	case Cantera::cSurf:
	case Cantera::cMetal:
	case Cantera::cStoichSubstance:
	case Cantera::cSemiconductor:
	case Cantera::cLatticeSolid:
	case Cantera::cLattice:
	case Cantera::cEdge:
	case Cantera::cIdealSolidSolnPhase:
	  m_isIdealSoln = true;
	  break;
	default:
	  m_isIdealSoln = false;
	};
      }
    } else {
      m_useCanteraCalls = false;
    }
  }
  /************************************************************************************/

  // Return a const ThermoPhase pointer corresponding to this phase
  /*
   *  @return pointer to the ThermoPhase.
   */
  const Cantera::ThermoPhase *vcs_VolPhase::ptrThermoPhase() const {
    return TP_ptr;
  }
  /************************************************************************************/

  double vcs_VolPhase::TotalMoles() const {
    return TMoles;
  }
  /************************************************************************************/

  double vcs_VolPhase::molefraction(int k) const {
    return Xmol[k];
  }
  /************************************************************************************/

  void vcs_VolPhase::setTotalMoles(double tmols)  {
    TMoles = tmols;
  }
  /************************************************************************************/

  // Return a string representing the equation of state
  /* 
   * The string is no more than 16 characters. 
   *  @param EOSType : integer value of the equation of state
   *
   * @return returns a string representing the EOS
   */
  std::string string16_EOSType(int EOSType) {
    char st[32];
    st[16] = '\0';
    switch (EOSType) {
    case VCS_EOS_CONSTANT:
      sprintf(st,"Constant        ");
      break;
    case VCS_EOS_IDEAL_GAS:
      sprintf(st,"Ideal Gas       ");
      break;
    case  VCS_EOS_STOICH_SUB:
      sprintf(st,"Stoich Sub      ");
      break;
    case VCS_EOS_IDEAL_SOLN:
      sprintf(st,"Ideal Soln      ");
      break;
    case VCS_EOS_DEBEYE_HUCKEL:
      sprintf(st,"Debeye Huckel   ");
      break;
    case VCS_EOS_REDLICK_KWONG:
      sprintf(st,"Redlick_Kwong   ");
      break;
    case VCS_EOS_REGULAR_SOLN:
      sprintf(st,"Regular Soln    ");
      break;
    default:
      sprintf(st,"UnkType: %-7d", EOSType);
      break;
    }
    st[16] = '\0';
    std::string sss=st;
    return sss;
  }
  /**********************************************************************/

  // Returns whether the phase is an ideal solution phase
  bool vcs_VolPhase::isIdealSoln() const {
    return m_isIdealSoln;
  }
  /**********************************************************************/

  // Returns whether the phase uses Cantera calls
  bool vcs_VolPhase::usingCanteraCalls() const {
    return m_useCanteraCalls;
  }
  /**********************************************************************/

}

