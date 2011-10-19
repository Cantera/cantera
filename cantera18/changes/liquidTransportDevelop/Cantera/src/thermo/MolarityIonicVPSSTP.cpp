/**
 *  @file MolarityIonicVPSSTP.cpp
 *   Definitions for intermediate ThermoPhase object for phases which
 *   employ excess gibbs free energy formulations
 *  (see \ref thermoprops 
 * and class \link Cantera::MolarityIonicVPSSTP MolarityIonicVPSSTP\endlink).
 *
 * Header file for a derived class of ThermoPhase that handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties that are further based upon expressions
 * for the excess gibbs free energy expressed as a function of
 * the mole fractions.
 */
/*
 * Copywrite (2009) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 *  $Date: 2009-11-09 16:36:49 -0700 (Mon, 09 Nov 2009) $
 *  $Revision: 255 $
 */


#include "MolarityIonicVPSSTP.h"

#include <cmath>

using namespace std;

namespace Cantera {
  //====================================================================================================================
  /*
   * Default constructor.
   *
   */
  MolarityIonicVPSSTP::MolarityIonicVPSSTP() :
    GibbsExcessVPSSTP(),
    PBType_(PBTYPE_PASSTHROUGH),
    numPBSpecies_(m_kk),
    indexSpecialSpecies_(-1),
    numCationSpecies_(0),
    numAnionSpecies_(0),
    numPassThroughSpecies_(0),
    neutralPBindexStart(0)
  {
  }
  //====================================================================================================================
  /*
   * Copy Constructor:
   *
   *  Note this stuff will not work until the underlying phase
   *  has a working copy constructor
   */
  MolarityIonicVPSSTP::MolarityIonicVPSSTP(const MolarityIonicVPSSTP &b) :
    GibbsExcessVPSSTP(),
    PBType_(PBTYPE_PASSTHROUGH),
    numPBSpecies_(m_kk),
    indexSpecialSpecies_(-1),
    numCationSpecies_(0),
    numAnionSpecies_(0),
    numPassThroughSpecies_(0),
    neutralPBindexStart(0)
  {
    *this = operator=(b);
  }
  //====================================================================================================================
  /*
   * operator=()
   *
   *  Note this stuff will not work until the underlying phase
   *  has a working assignment operator
   */
  MolarityIonicVPSSTP& MolarityIonicVPSSTP::
  operator=(const MolarityIonicVPSSTP &b) {
    if (&b != this) {
      GibbsExcessVPSSTP::operator=(b);
    }

    PBType_                     = b.PBType_;
    numPBSpecies_               = b.numPBSpecies_;
    indexSpecialSpecies_        = b.indexSpecialSpecies_;
    PBMoleFractions_            = b.PBMoleFractions_;
    cationList_                 = b.cationList_;
    numCationSpecies_           = b.numCationSpecies_;
    anionList_                  = b.anionList_;
    numAnionSpecies_            = b.numAnionSpecies_;
    passThroughList_            = b.passThroughList_;
    numPassThroughSpecies_      = b.numPassThroughSpecies_;
    neutralPBindexStart         = b.neutralPBindexStart;
    moleFractionsTmp_           = b.moleFractionsTmp_;

    return *this;
  }
  //====================================================================================================================
  /**
   *
   * ~MolarityIonicVPSSTP():   (virtual)
   *
   * Destructor: does nothing:
   *
   */
  MolarityIonicVPSSTP::~MolarityIonicVPSSTP() {
  }

  /*
   * This routine duplicates the current object and returns
   * a pointer to ThermoPhase.
   */
  ThermoPhase* 
  MolarityIonicVPSSTP::duplMyselfAsThermoPhase() const {
    MolarityIonicVPSSTP* mtp = new MolarityIonicVPSSTP(*this);
    return (ThermoPhase *) mtp;
  }

  /*
   *  -------------- Utilities -------------------------------
   */
  //====================================================================================================================
 
  // Equation of state type flag.
  /*
   * The ThermoPhase base class returns
   * zero. Subclasses should define this to return a unique
   * non-zero value. Known constants defined for this purpose are
   * listed in mix_defs.h. The MolarityIonicVPSSTP class also returns
   * zero, as it is a non-complete class.
   */
  int MolarityIonicVPSSTP::eosType() const { 
    return 0;
  }

 

  /*
   * ------------ Molar Thermodynamic Properties ----------------------
   */


  /*
   * - Activities, Standard States, Activity Concentrations -----------
   */

  //====================================================================================================================
  doublereal MolarityIonicVPSSTP::standardConcentration(int k) const {
    err("standardConcentration");
    return -1.0;
  }
  //====================================================================================================================
  doublereal MolarityIonicVPSSTP::logStandardConc(int k) const {
    err("logStandardConc");
    return -1.0;
  }
  //====================================================================================================================

  void MolarityIonicVPSSTP::getElectrochemPotentials(doublereal* mu) const {
    getChemPotentials(mu);
    double ve = Faraday * electricPotential();
    for (int k = 0; k < m_kk; k++) {
      mu[k] += ve*charge(k);
    }
  }
  //====================================================================================================================
  void MolarityIonicVPSSTP::calcPseudoBinaryMoleFractions() const {
    int k;
    int kCat;
    int kMax;
    doublereal sumCat; 
    doublereal sumAnion;
    doublereal chP, chM;
    doublereal sum = 0.0;
    doublereal sumMax;
    switch (PBType_) {
    case PBTYPE_PASSTHROUGH:
      for (k = 0; k < m_kk; k++) {
	PBMoleFractions_[k] = moleFractions_[k];
      }
      break;
    case PBTYPE_SINGLEANION:
      sumCat = 0.0;
      sumAnion = 0.0;
      for (k = 0; k < m_kk; k++) {
	moleFractionsTmp_[k] = moleFractions_[k];
      }
      kMax = -1;
      sumMax = 0.0;
      for (k = 0; k < (int) cationList_.size(); k++) {
	kCat = cationList_[k];
	chP = m_speciesCharge[kCat];
	if (moleFractions_[kCat] > sumMax) {
	  kMax = k;
	  sumMax = moleFractions_[kCat];
	}
	sumCat += chP * moleFractions_[kCat];
      }
      k = anionList_[0];
      chM = m_speciesCharge[k];
      sumAnion = moleFractions_[k] * chM;
      sum = sumCat - sumAnion;
      if (fabs(sum) > 1.0E-16) {
	moleFractionsTmp_[cationList_[kMax]] -= sum / m_speciesCharge[kMax];
	sum = 0.0;
	for (k = 0; k < numCationSpecies_; k++) {
	  sum +=  moleFractionsTmp_[k];
	}
	for (k = 0; k < numCationSpecies_; k++) {
	  moleFractionsTmp_[k]/= sum;
	}
      }

      for (k = 0; k < numCationSpecies_; k++) {
        PBMoleFractions_[k] = moleFractionsTmp_[cationList_[k]];
      }
      for (k = 0; k <  numPassThroughSpecies_; k++) {
	PBMoleFractions_[neutralPBindexStart + k] = moleFractions_[passThroughList_[k]];
      }
      
      sum = fmaxx(0.0, PBMoleFractions_[0]);
      for (k = 1; k < numPBSpecies_; k++) {
        sum += PBMoleFractions_[k];
	
      }
      for (k = 0; k < numPBSpecies_; k++) {
        PBMoleFractions_[k] /= sum;
      }

      break;
    case PBTYPE_SINGLECATION:
      throw CanteraError("eosType", "Unknown type");
     
      break;
     
    case PBTYPE_MULTICATIONANION:
      throw CanteraError("eosType", "Unknown type");
     
      break;
    default:
      throw CanteraError("eosType", "Unknown type");
      break;

    } 
  }
  //====================================================================================================================
  /*
   * ------------ Partial Molar Properties of the Solution ------------
   */
  //====================================================================================================================
  doublereal MolarityIonicVPSSTP::err(std::string msg) const {
    throw CanteraError("MolarityIonicVPSSTP","Base class method "
		       +msg+" called. Equation of state type: "+int2str(eosType()));
    return 0;
  }
  //====================================================================================================================
  /*
   * @internal Initialize. This method is provided to allow
   * subclasses to perform any initialization required after all
   * species have been added. For example, it might be used to
   * resize internal work arrays that must have an entry for
   * each species.  The base class implementation does nothing,
   * and subclasses that do not require initialization do not
   * need to overload this method.  When importing a CTML phase
   * description, this method is called just prior to returning
   * from function importPhase.
   *
   * @see importCTML.cpp
   */
  void MolarityIonicVPSSTP::initThermo() {
    GibbsExcessVPSSTP::initThermo();
    initLengths();
    /*
     *  Go find the list of cations and anions
     */
    double ch;
    numCationSpecies_ = 0.0;
    cationList_.clear();
    anionList_.clear();
    passThroughList_.clear();
    for (int k = 0; k < m_kk; k++) {
      ch = m_speciesCharge[k];
      if (ch > 0.0) {
	cationList_.push_back(k);
	numCationSpecies_++;
      } else if (ch < 0.0) {
	anionList_.push_back(k);
	numAnionSpecies_++;
      } else {
	passThroughList_.push_back(k);
	numPassThroughSpecies_++;
      }
    }
    numPBSpecies_ = numCationSpecies_ + numAnionSpecies_ - 1;
    neutralPBindexStart = numPBSpecies_;
    PBType_ = PBTYPE_MULTICATIONANION;
    if (numAnionSpecies_ == 1) {
      PBType_ = PBTYPE_SINGLEANION;
    } else if (numCationSpecies_ == 1) {
      PBType_ = PBTYPE_SINGLECATION;
    }
    if (numAnionSpecies_ == 0 && numCationSpecies_ == 0) {
      PBType_ = PBTYPE_PASSTHROUGH;
    }
  }
  //====================================================================================================================
  //   Initialize lengths of local variables after all species have been identified.
  void  MolarityIonicVPSSTP::initLengths() {
    m_kk = nSpecies();
    moleFractionsTmp_.resize(m_kk);    
  }
  //====================================================================================================================
  /*
   * initThermoXML()                (virtual from ThermoPhase)
   *   Import and initialize a ThermoPhase object
   *
   * @param phaseNode This object must be the phase node of a
   *             complete XML tree
   *             description of the phase, including all of the
   *             species data. In other words while "phase" must
   *             point to an XML phase object, it must have
   *             sibling nodes "speciesData" that describe
   *             the species in the phase.
   * @param id   ID of the phase. If nonnull, a check is done
   *             to see if phaseNode is pointing to the phase
   *             with the correct id. 
   */
  void MolarityIonicVPSSTP::initThermoXML(XML_Node& phaseNode, std::string id) {

 
    GibbsExcessVPSSTP::initThermoXML(phaseNode, id);
  }
  //====================================================================================================================
  /*
   * Format a summary of the mixture state for output.
   */           
  std::string MolarityIonicVPSSTP::report(bool show_thermo) const {
    char p[800];
    string s = "";
    try {
      if (name() != "") {
	sprintf(p, " \n  %s:\n", name().c_str());
	s += p;
      }
      sprintf(p, " \n       temperature    %12.6g  K\n", temperature());
      s += p;
      sprintf(p, "          pressure    %12.6g  Pa\n", pressure());
      s += p;
      sprintf(p, "           density    %12.6g  kg/m^3\n", density());
      s += p;
      sprintf(p, "  mean mol. weight    %12.6g  amu\n", meanMolecularWeight());
      s += p;

      doublereal phi = electricPotential();
      sprintf(p, "         potential    %12.6g  V\n", phi);
      s += p;

      int kk = nSpecies();
      array_fp x(kk);
      array_fp molal(kk);
      array_fp mu(kk);
      array_fp muss(kk);
      array_fp acMolal(kk);
      array_fp actMolal(kk);
      getMoleFractions(&x[0]);
   
      getChemPotentials(&mu[0]);
      getStandardChemPotentials(&muss[0]);
      getActivities(&actMolal[0]);
 

      if (show_thermo) {
        sprintf(p, " \n");
        s += p;
        sprintf(p, "                          1 kg            1 kmol\n");
        s += p;
        sprintf(p, "                       -----------      ------------\n");
        s += p;
        sprintf(p, "          enthalpy    %12.6g     %12.4g     J\n", 
		enthalpy_mass(), enthalpy_mole());
        s += p;
        sprintf(p, "   internal energy    %12.6g     %12.4g     J\n", 
		intEnergy_mass(), intEnergy_mole());
        s += p;
        sprintf(p, "           entropy    %12.6g     %12.4g     J/K\n", 
		entropy_mass(), entropy_mole());
        s += p;
        sprintf(p, "    Gibbs function    %12.6g     %12.4g     J\n", 
		gibbs_mass(), gibbs_mole());
        s += p;
        sprintf(p, " heat capacity c_p    %12.6g     %12.4g     J/K\n", 
		cp_mass(), cp_mole());
        s += p;
        try {
	  sprintf(p, " heat capacity c_v    %12.6g     %12.4g     J/K\n", 
		  cv_mass(), cv_mole());
	  s += p;
        }
        catch(CanteraError) {
	  sprintf(p, " heat capacity c_v    <not implemented>       \n");
	  s += p;
        }
      }
  
    } catch (CanteraError) {
      ;
    }
    return s;
  }
  //====================================================================================================================
}

