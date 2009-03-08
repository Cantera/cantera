/**
 *
 *  @file MolalityVPSSTP.cpp
 */
/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 *  $Author: hkmoffa $
 *  $Date: 2006/06/19 23:14:12 $
 *  $Revision: 1.5 $
 */
#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif


#include "MolalityVPSSTP.h"


namespace Cantera {

  /*
   * Default constructor.
   *
   * This doesn't do much more than initialize constants with
   * default values for water at 25C. Water molecular weight 
   * comes from the default elements.xml file. It actually
   * differs slightly from the IAPWS95 value of 18.015268. However,
   * density conservation and therefore element conservation
   * is the more important principle to follow.
   */
  MolalityVPSSTP::MolalityVPSSTP() :
    VPStandardStateTP(),
    m_indexSolvent(0),
    m_weightSolvent(18.01528),
    m_xmolSolventMIN(0.01),
    m_Mnaught(18.01528E-3)
  {
  }

  /**
   * Copy Constructor:
   *
   *  Note this stuff will not work until the underlying phase
   *  has a working copy constructor
   */
  MolalityVPSSTP::MolalityVPSSTP(const MolalityVPSSTP &b) :
    VPStandardStateTP(),
    m_indexSolvent(b.m_indexSolvent),
    m_xmolSolventMIN(b.m_xmolSolventMIN),
    m_Mnaught(b.m_Mnaught),
    m_molalities(b.m_molalities)
  {
    throw CanteraError("MolalityVPSSTP::operator=()",
		       "Not Implemented Fully");
    *this = operator=(b);
  }

  /*
   * operator=()
   *
   *  Note this stuff will not work until the underlying phase
   *  has a working assignment operator
   */
  MolalityVPSSTP& MolalityVPSSTP::
  operator=(const MolalityVPSSTP &b) {
    if (&b != this) {
      VPStandardStateTP::operator=(b);
      m_indexSolvent     = b.m_indexSolvent;
      m_weightSolvent    = b.m_weightSolvent;
      m_xmolSolventMIN   = b.m_xmolSolventMIN;
      m_Mnaught          = b.m_Mnaught;
      m_molalities       = b.m_molalities;
    }
    throw CanteraError("MolalityVPSSTP::operator=()",
		       "Not Implemented Fully");
    return *this;
  }

  /**
   *
   * ~MolalityVPSSTP():   (virtual)
   *
   * Destructor: does nothing:
   *
   */
  MolalityVPSSTP::~MolalityVPSSTP() {
  }

  /**
   * This routine duplicates the current object and returns
   * a pointer to ThermoPhase.
   */
  ThermoPhase* 
  MolalityVPSSTP::duplMyselfAsThermoPhase() {
    MolalityVPSSTP* mtp = new MolalityVPSSTP(*this);
    return (ThermoPhase *) mtp;
  }

  /*
   *  -------------- Utilities -------------------------------
   */

  /*
   * setSolvent():
   *  Utilities for Solvent ID and Molality 
   *  Here we also calculate and store the molecular weight
   *  of the solvent and the m_Mnaught parameter.
   *  @param k index of the solvent.
   */
  void MolalityVPSSTP::setSolvent(int k) {
    if (k < 0 || k >= m_kk) {
      throw CanteraError("MolalityVPSSTP::setSolute ", 
			 "bad value");
    }
    m_indexSolvent = k;
    m_weightSolvent = molecularWeight(k);
    m_Mnaught = m_weightSolvent / 1000.;
  }
    
  /*
   * return the solvent id index number.
   */
  int MolalityVPSSTP::solventIndex() const {
    return m_indexSolvent;
  }

  /**
   * Sets the minimum mole fraction in the molality formulation. The
   * minimum mole fraction must be in the range 0 to 0.9.
   */
  void  MolalityVPSSTP::
  setMoleFSolventMin(doublereal xmolSolventMIN) {
    if (xmolSolventMIN <= 0.0) {
      throw CanteraError("MolalityVPSSTP::setSolute ", "trouble");
    } else if (xmolSolventMIN > 0.9) {
      throw CanteraError("MolalityVPSSTP::setSolute ", "trouble");
    }
    m_xmolSolventMIN = xmolSolventMIN;
  }

  /**
   * Returns the minimum mole fraction in the molality formulation.
   */
  doublereal MolalityVPSSTP::moleFSolventMin() const {
    return m_xmolSolventMIN;
  }

  /**
   * calcMolalities():
   *   We calculate the vector of molalities of the species
   *   in the phase and store the result internally:
   * \f[
   *     m_i = (n_i) / (1000 * M_o * n_{o,p})
   * \f]
   *    where 
   *    - \f$ M_o \f$ is the molecular weight of the solvent
   *    - \f$ n_o \f$ is the mole fraction of the solvent
   *    - \f$ n_i \f$ is the mole fraction of the solute.
   *    - \f$ n_{o,p} = max (n_{o, min}, n_o) \f$
   *    - \f$ n_{o,min} \f$ = minimum mole fraction of solvent allowed
   *              in the denominator.
   */
  void MolalityVPSSTP::calcMolalities() const {
    getMoleFractions(DATA_PTR(m_molalities));
    double xmolSolvent = m_molalities[m_indexSolvent];
    if (xmolSolvent < m_xmolSolventMIN) {
      xmolSolvent = m_xmolSolventMIN;
    }
    double denomInv = 1.0/ (m_Mnaught * xmolSolvent);
    for (int k = 0; k < m_kk; k++) {
      m_molalities[k] *= denomInv;
    }
  }

  /**
   * getMolalities():
   *   We calculate the vector of molalities of the species
   *   in the phase
   * \f[
   *     m_i = (n_i) / (1000 * M_o * n_{o,p})
   * \f]
   *    where 
   *    - \f$ M_o \f$ is the molecular weight of the solvent
   *    - \f$ n_o \f$ is the mole fraction of the solvent
   *    - \f$ n_i \f$ is the mole fraction of the solute.
   *    - \f$ n_{o,p} = max (n_{o, min}, n_o) \f$
   *    - \f$ n_{o,min} \f$ = minimum mole fraction of solvent allowed
   *              in the denominator.
   */
  void MolalityVPSSTP::getMolalities(doublereal * const molal) const {
    calcMolalities();
    for (int k = 0; k < m_kk; k++) {
      molal[k] = m_molalities[k];
    }
  }

  /**
   * setMolalities():
   *   We are supplied with the molalities of all of the
   *   solute species. We then calculate the mole fractions of all
   *   species and update the ThermoPhase object.
   *
   *     m_i = (n_i) / (W_o/1000 * n_o_p)
   *
   *    where M_o is the molecular weight of the solvent
   *    n_o is the mole fraction of the solvent
   *    n_i is the mole fraction of the solute.
   *    n_o_p = max (n_o_min, n_o)
   *    n_o_min = minimum mole fraction of solvent allowed
   *              in the denominator.
   */
  void MolalityVPSSTP::setMolalities(const doublereal * const molal) {
	
    double Lsum = 1.0 / m_Mnaught;
    for (int k = 0; k < m_kk; k++) {
      if (k != m_indexSolvent) { 
	m_molalities[k] = molal[k];
	Lsum += molal[k];
      }
    }
    double tmp = 1.0 / Lsum;
    m_molalities[m_indexSolvent] = tmp / m_Mnaught;
    double sum = m_molalities[m_indexSolvent];
    for (int k = 0; k < m_kk; k++) {
      if (k != m_indexSolvent) { 
	m_molalities[k] = tmp * molal[k];
	sum += m_molalities[k];
      }
    }
    if (sum != 1.0) {
      tmp = 1.0 / sum;
      for (int k = 0; k < m_kk; k++) {
	m_molalities[k] *= tmp;
      }
    }
    setMoleFractions(DATA_PTR(m_molalities));
    /*
     * Essentially we don't trust the input: We calculate
     * the molalities from the mole fractions that we 
     * just obtained.
     */
    calcMolalities();
  }
    
  /*
   * setMolalitiesByName()
   *
   *  This routine sets the molalities by name
   *  HKM -> Might need to be more complicated here, setting
   *         neutrals so that the existing mole fractions are
   *         preserved.
   */
  void MolalityVPSSTP::setMolalitiesByName(compositionMap& mMap) {
    int kk = nSpecies();
    doublereal x;
    /*
     * Get a vector of mole fractions
     */
    vector_fp mf(kk, 0.0);
    getMoleFractions(DATA_PTR(mf));
    double xmolS = mf[m_indexSolvent];
    double xmolSmin = max(xmolS, m_xmolSolventMIN);
    compositionMap::iterator p;
    for (int k = 0; k < kk; k++) {
      p = mMap.find(speciesName(k));
      if (p != mMap.end()) {
	x = mMap[speciesName(k)];
	if (x > 0.0) {
	  mf[k] = x * m_Mnaught * xmolSmin;
	}
      }
    }
    /*
     * check charge neutrality
     */
    int largePos = -1;
    double cPos = 0.0;
    int largeNeg = -1;
    double cNeg = 0.0;
    double sum = 0.0;
    for (int k = 0; k < kk; k++) {
      double ch = charge(k);
      if (mf[k] > 0.0) {
	if (ch > 0.0) {
	  if (ch * mf[k] > cPos) {
	    largePos = k;
	    cPos = ch * mf[k];
	  }
	}
	if (ch < 0.0) {
	  if (fabs(ch) * mf[k] > cNeg) {
	    largeNeg = k;
	    cNeg = fabs(ch) * mf[k];
	  }
	}
      }
      sum += mf[k] * ch;
    }
    if (sum != 0.0) {
      if (sum > 0.0) {
	if (cPos > sum) {
	  mf[largePos] -= sum / charge(largePos);
	} else {
	  throw CanteraError("MolalityVPSSTP:setMolalitiesbyName",
			     "unbalanced charges");
	}
      } else {
	if (cNeg > (-sum)) {
	  mf[largeNeg] -= (-sum) / fabs(charge(largeNeg));
	} else {
	  throw CanteraError("MolalityVPSSTP:setMolalitiesbyName",
			     "unbalanced charges");
	}
      }
            
    }
    sum = 0.0;
    for (int k = 0; k < kk; k++) {
      sum += mf[k];
    }
    sum = 1.0/sum;
    for (int k = 0; k < kk; k++) {
      mf[k] *= sum;
    }
    setMoleFractions(DATA_PTR(mf));
    /*
     * After we formally set the mole fractions, we
     * calculate the molalities again and store it in
     * this object.
     */
    calcMolalities();
  }

  /*
   * setMolalitiesByNames()
   *
   *   Set the molalities of the solutes by name
   */
  void MolalityVPSSTP::setMolalitiesByName(const string& x) {
    compositionMap xx;
    int kk = nSpecies();
    for (int k = 0; k < kk; k++) {
      xx[speciesName(k)] = -1.0;
    }
    parseCompString(x, xx);
    setMolalitiesByName(xx);
  }


  /*
   * ------------ Molar Thermodynamic Properties ----------------------
   */


  /*
   * - Activities, Standard States, Activity Concentrations -----------
   */

  /**
   * This method returns the activity convention.
   * Currently, there are two activity conventions
   *  Molar-based activities
   *       Unit activity of species at either a hypothetical pure
   *       solution of the species or at a hypothetical
   *       pure ideal solution at infinite dilution
   *   cAC_CONVENTION_MOLAR 0
   *      - default
   *  
   *  Molality based activities
   *       (unit activity of solutes at a hypothetical 1 molal
   *        solution referenced to infinite dilution at all
   *        pressures and temperatures).
   *       (solvent is still on molar basis).
   *   cAC_CONVENTION_MOLALITY 1
   *
   *  We set the convention to molality here.
   */
  int MolalityVPSSTP::activityConvention() const {
    return cAC_CONVENTION_MOLALITY;
  }

  /**
   * Get the array of non-dimensional activity coefficients at
   * the current solution temperature, pressure, and
   * solution concentration.
   * These are mole fraction based activity coefficients. In this
   * object, their calculation is based on translating the values
   * of Molality based activity coefficients.
   *  See Denbigh p. 278 for a thorough discussion.
   *
   * Note, the solvent is treated differently. getMolalityActivityCoeff()
   * returns the molar based solvent activity coefficient already.
   * Therefore, we do not have to divide by x_s here.
   */
  void MolalityVPSSTP::getActivityCoefficients(doublereal* ac) const {
    getMolalityActivityCoefficients(ac);
    double xmolSolvent = moleFraction(m_indexSolvent);
    if (xmolSolvent < m_xmolSolventMIN) {
      xmolSolvent = m_xmolSolventMIN;
    }
    for (int k = 0; k < m_kk; k++) {
      if (k != m_indexSolvent) {
	ac[k] /= xmolSolvent;
      }
    }
  }
    
  /**
   * osmotic coefficient:
   * 
   *  Calculate the osmotic coefficient of the solvent. Note there
   *  are lots of definitions of the osmotic coefficient floating
   *  around. We use the one defined in the Pitzer's book:
   *  (Activity Coeff in Electrolyte Solutions, K. S. Pitzer
   *   CRC Press, Boca Raton, 1991, p. 85, Eqn. 28).
   *
   *        Definition:
   *         - sum(m_i) * Mnaught * oc = ln(activity_solvent)
   */
  doublereal MolalityVPSSTP::osmoticCoefficient() const {
    /*
     * First, we calculate the activities all over again
     */
    vector_fp act(m_kk);
    getActivities(DATA_PTR(act));
    /*
     * Then, we calculate the sum of the solvent molalities
     */
    double sum = 0;
    for (int k = 0; k < m_kk; k++) {
      if (k != m_indexSolvent) {
	sum += MAX(m_molalities[k], 0.0);
      }
    }
    double oc = 1.0;
    double lac = log(act[m_indexSolvent]);
    if (sum > 1.0E-200) {
      oc = - lac / (m_Mnaught * sum);
    }
    return oc;
  }

  /*
   * ------------ Partial Molar Properties of the Solution ------------
   */


  doublereal MolalityVPSSTP::err(string msg) const {
    throw CanteraError("MolalityVPSSTP","Base class method "
		       +msg+" called. Equation of state type: "+int2str(eosType()));
    return 0;
  }

  /**
   * Returns the units of the standard and general concentrations
   * Note they have the same units, as their divisor is 
   * defined to be equal to the activity of the kth species
   * in the solution, which is unitless.
   *
   * This routine is used in print out applications where the
   * units are needed. Usually, MKS units are assumed throughout
   * the program and in the XML input files. 
   *
   * On return uA contains the powers of the units (MKS assumed)
   * of the standard concentrations and generalized concentrations
   * for the kth species.
   *
   *  uA[0] = kmol units - default  = 1
   *  uA[1] = m    units - default  = -nDim(), the number of spatial
   *                                dimensions in the Phase class.
   *  uA[2] = kg   units - default  = 0;
   *  uA[3] = Pa(pressure) units - default = 0;
   *  uA[4] = Temperature units - default = 0;
   *  uA[5] = time units - default = 0
   */
  void MolalityVPSSTP::getUnitsStandardConc(double *uA, int k, int sizeUA) {
    for (int i = 0; i < sizeUA; i++) {
      if (i == 0) uA[0] = 1.0;
      if (i == 1) uA[1] = -nDim();
      if (i == 2) uA[2] = 0.0;
      if (i == 3) uA[3] = 0.0;
      if (i == 4) uA[4] = 0.0;
      if (i == 5) uA[5] = 0.0;
    }
  }
    

  /*
   * Set the thermodynamic state.
   */
  void MolalityVPSSTP::setStateFromXML(const XML_Node& state) {
    VPStandardStateTP::setStateFromXML(state);
    string comp = getString(state,"soluteMolalities");
    if (comp != "") {
      setMolalitiesByName(comp);
    }
    if (state.hasChild("pressure")) {
      double p = getFloat(state, "pressure", "pressure");
      setPressure(p);
    }
  }

  /** 
   * Set the temperature (K), pressure (Pa), and molalities
   * (gmol kg-1) of the solutes
   */
  void MolalityVPSSTP::setState_TPM(doublereal t, doublereal p, 
				    const doublereal * const molalities) {
    setMolalities(molalities);
    setTemperature(t);
    setPressure(p);
  }

  /** Set the temperature (K), pressure (Pa), and molalities.  */
  void MolalityVPSSTP::setState_TPM(doublereal t, doublereal p, compositionMap& m) {
    setMolalitiesByName(m);
    setTemperature(t);
    setPressure(p);
  }

  /** Set the temperature (K), pressure (Pa), and molality.  */
  void MolalityVPSSTP::setState_TPM(doublereal t, doublereal p, const string& m) {
    setMolalitiesByName(m);
    setTemperature(t);
    setPressure(p);
  }


  /**
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
  void MolalityVPSSTP::initThermo() {
    initLengths();
    VPStandardStateTP::initThermo();

    /*
     * The solvent defaults to species 0
     */
    setSolvent(0);
  }

  void  MolalityVPSSTP::initLengths() {
    int m_kk = nSpecies();
    m_molalities.resize(m_kk);
  }

  /**
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
  void MolalityVPSSTP::initThermoXML(XML_Node& phaseNode, string id) {

    initLengths();
    /*
     * The solvent defaults to species 0
     */
    setSolvent(0);

    VPStandardStateTP::initThermoXML(phaseNode, id);
  }
  

}




