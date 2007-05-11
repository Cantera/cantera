/**
 *  @file VPStandardStateTP.cpp
 * Definition file for a derived class of ThermoPhase that handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties (see \ref thermoprops and
 * class \link Cantera::VPStandardStateTP VPStandardStateTP\endlink).
 */
/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 *  $Author$
 *  $Date$
 *  $Revision$
 */

// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "VPStandardStateTP.h"

using namespace std;

namespace Cantera {

  /*
   * Default constructor
   */
  VPStandardStateTP::VPStandardStateTP() :
    ThermoPhase(),
    m_Pcurrent(OneAtm),
    m_tlast(-1.0),
    m_tlast_ref(-1.0),
    m_plast(-1.0),
    m_p0(OneAtm),
    m_useTmpRefStateStorage(false),
    m_useTmpStandardStateStorage(false)
  {
  }

  /*
   * Copy Constructor:
   *
   *  Note this stuff will not work until the underlying phase
   *  has a working copy constructor.
   *
   *  The copy constructor just calls the assignment operator
   *  to do the heavy lifting.
   */
  VPStandardStateTP::VPStandardStateTP(const VPStandardStateTP &b) :
    ThermoPhase(),
    m_Pcurrent(OneAtm),
    m_tlast(-1.0),
    m_tlast_ref(-1.0),
    m_plast(-1.0),
    m_p0(OneAtm),
    m_useTmpRefStateStorage(false),
    m_useTmpStandardStateStorage(false)
  {
    *this = b;
  }

  /*
   * operator=()
   *
   *  Note this stuff will not work until the underlying phase
   *  has a working assignment operator
   */
  VPStandardStateTP& VPStandardStateTP::
  operator=(const VPStandardStateTP &b) {
    if (&b != this) {
      /*
       * Mostly, this is a passthrough to the underlying
       * assignment operator for the ThermoPhae parent object.
       */
      ThermoPhase::operator=(b);
      /*
       * However, we have to handle data that we own.
       */
      m_Pcurrent  = b.m_Pcurrent;
      m_tlast     = b.m_tlast;
      m_tlast_ref = b.m_tlast_ref;
      m_plast     = b.m_plast;
      m_p0        = b.m_p0;
      m_useTmpRefStateStorage = b.m_useTmpRefStateStorage;
      m_h0_RT     = b.m_h0_RT;
      m_cp0_R     = b.m_cp0_R;
      m_g0_RT     = b.m_g0_RT;
      m_s0_R      = b.m_s0_R;
      m_useTmpStandardStateStorage = b.m_useTmpStandardStateStorage;
      m_hss_RT     = b.m_hss_RT;
      m_cpss_R     = b.m_cpss_R;
      m_gss_RT     = b.m_gss_RT;
      m_sss_R      = b.m_sss_R;
      m_Vss        = b.m_Vss;
    }
    return *this;
  }

  /*
   * ~VPStandardStateTP():   (virtual)
   *
   * This destructor does nothing. All of the owned objects
   * handle themselves.
   */
  VPStandardStateTP::~VPStandardStateTP() {
  }

  /*
   * Duplication function.
   *  This calls the copy constructor for this object.
   */
  ThermoPhase* VPStandardStateTP::duplMyselfAsThermoPhase() {
    VPStandardStateTP* vptp = new VPStandardStateTP(*this);
    return (ThermoPhase *) vptp;
  }
 
  
  /*
   * ------------Molar Thermodynamic Properties -------------------------
   */
  
  
  doublereal VPStandardStateTP::err(std::string msg) const {
    throw CanteraError("VPStandardStateTP","Base class method "
		       +msg+" called. Equation of state type: "+int2str(eosType()));
    return 0;
  }
  
  /*
   * ---- Partial Molar Properties of the Solution -----------------
   */
  
  /*
   * Get the array of non-dimensional species chemical potentials
   * These are partial molar Gibbs free energies.
   * \f$ \mu_k / \hat R T \f$.
   * Units: unitless
   *
   * We close the loop on this function, here, calling
   * getChemPotentials() and then dividing by RT.
   */
  void VPStandardStateTP::getChemPotentials_RT(doublereal* muRT) const{
    getChemPotentials(muRT);
    doublereal invRT = 1.0 / _RT();
    for (int k = 0; k < m_kk; k++) {
      muRT[k] *= invRT;
    }
  }
  
  /*
   * ----- Thermodynamic Values for the Species Standard States States ----
   */
  void VPStandardStateTP::getStandardChemPotentials(doublereal* g) const {
    getGibbs_RT(g);
    doublereal RT = _RT();
    for (int k = 0; k < m_kk; k++) {
      g[k] *= RT;
    }
  } 
  
  void VPStandardStateTP::getEnthalpy_RT(doublereal* hrt) const {
    if (m_useTmpStandardStateStorage) {
      _updateStandardStateThermo();
      copy(m_hss_RT.begin(), m_hss_RT.end(), hrt);
    } else {
      err("getEnthalpy_RT ERROR: Must be overwritten in child classes");
      _updateStandardStateThermo();
    }
  }

  void VPStandardStateTP::getEntropy_R(doublereal* srt) const {
    if (m_useTmpStandardStateStorage) {
      _updateStandardStateThermo();
      copy(m_sss_R.begin(), m_sss_R.end(), srt);
    } else {
      err("getEntropy_R ERROR: Must be overwritten in child classes");
      _updateStandardStateThermo();
    }
  }
  
  void VPStandardStateTP::getGibbs_RT(doublereal* grt) const {
    if (m_useTmpStandardStateStorage) {
      _updateStandardStateThermo();
      copy(m_gss_RT.begin(), m_gss_RT.end(), grt);
    } else {
      err("getGibbs_RT ERROR: Must be overwritten in child classes");
      _updateStandardStateThermo();
    }
  }
  
  void VPStandardStateTP::getPureGibbs(doublereal* g) const {
    getGibbs_RT(g);
    doublereal RT = _RT();
    for (int k = 0; k < m_kk; k++) {
      g[k] *= RT;
    }
  }

  void VPStandardStateTP::getIntEnergy_RT(doublereal* urt) const {
    if (m_useTmpStandardStateStorage) {
      _updateStandardStateThermo();
      copy(m_hss_RT.begin(), m_hss_RT.end(), urt);
      doublereal RT = _RT();
      doublereal tmp = pressure() / RT;
      for (int k = 0; k < m_kk; k++) {
	urt[k] -= tmp * m_Vss[k];
      }
    } else {
      err("getIntEnergy_RT ERROR: Must be overwritten in child classes");
      _updateStandardStateThermo();
    }
  }

  void VPStandardStateTP::getCp_R(doublereal* cpr) const {
    if (m_useTmpStandardStateStorage) {
      _updateStandardStateThermo();
      copy(m_cpss_R.begin(), m_cpss_R.end(), cpr);
    } else {
      err("getCp_R ERROR: Must be overwritten in child classes");
      _updateStandardStateThermo();
    }
  }

  void VPStandardStateTP::getStandardVolumes(doublereal *vol) const {
   if (m_useTmpStandardStateStorage) {
     _updateStandardStateThermo();
     copy(m_Vss.begin(), m_Vss.end(), vol);
   } else {
     err("getStandardVolumes ERROR: Must be overwritten in child classes");
     _updateStandardStateThermo();
   }
  }

  /*
   * ----- Thermodynamic Values for the Species Reference States ----
   */

  /*
   *  Returns the vector of nondimensional enthalpies of the
   *  reference state at the current temperature of the solution and
   *  the reference pressure for the species.
   */
  void VPStandardStateTP::getEnthalpy_RT_ref(doublereal *hrt) const {
    if (m_useTmpRefStateStorage) {
      /*
       * Call the function that makes sure the local copy of the
       * species reference thermo functions are up to date for the
       * current temperature.
       */
      _updateRefStateThermo();
      /*
       * Copy the enthalpy function into return vector.
       */
      copy(m_h0_RT.begin(), m_h0_RT.end(), hrt);
    } else if (m_useTmpStandardStateStorage) {
      _updateStandardStateThermo(m_p0);
      copy(m_hss_RT.begin(), m_hss_RT.end(), hrt);
    } else {
      err("getEnthalpy_RT_ref() ERROR: not handled");
    }
  }
    
  /*
   *  Returns the vector of nondimensional
   *  enthalpies of the reference state at the current temperature
   *  of the solution and the reference pressure for the species.
   */
  void VPStandardStateTP::getGibbs_RT_ref(doublereal *grt) const {
    if (m_useTmpRefStateStorage) {
      /*
       * Call the function that makes sure the local copy of 
       * the species reference thermo functions are up to date
       * for the current temperature.
       */
      _updateRefStateThermo();
      /*
       * Copy the gibbs function into return vector.
       */
      copy(m_g0_RT.begin(), m_g0_RT.end(), grt);
    } else if (m_useTmpStandardStateStorage) {
      _updateStandardStateThermo(m_p0);
      copy(m_gss_RT.begin(), m_gss_RT.end(), grt);
    } else {
      err("getGibbs_RT_ref() ERROR: not handled");
    }
  }    
        
  /*
   *  Returns the vector of the
   *  gibbs function of the reference state at the current temperature
   *  of the solution and the reference pressure for the species.
   *  units = J/kmol
   *
   *  This is filled in here so that derived classes don't have to
   *  take care of it.
   */
  void VPStandardStateTP::getGibbs_ref(doublereal *g) const {
    getGibbs_RT_ref(g);
    double RT = _RT();
    for (int k = 0; k < m_kk; k++) {
      g[k] *= RT;
    }
  }

  /*
   *  Returns the vector of nondimensional
   *  entropies of the reference state at the current temperature
   *  of the solution and the reference pressure for the species.
   */
  void VPStandardStateTP::getEntropy_R_ref(doublereal *er) const {
    if (m_useTmpRefStateStorage) {
      /*
       * Call the function that makes sure the local copy of 
       * the species reference thermo functions are up to date
       * for the current temperature.
       */
      _updateRefStateThermo();
      /*
       * Copy the gibbs function into return vector.
       */
      copy(m_s0_R.begin(), m_s0_R.end(), er);
    } else if (m_useTmpStandardStateStorage) {
      _updateStandardStateThermo(m_p0);
      copy(m_sss_R.begin(), m_sss_R.end(), er);
    } else {
      err("getEntropy_R_ref() ERROR: not handled");
    }
  }
     
  /*
   *  Returns the vector of nondimensional
   *  constant pressure heat capacities of the reference state
   *  at the current temperature of the solution
   *  and reference pressure for the species.
   */
  void VPStandardStateTP::getCp_R_ref(doublereal *cpr) const {
    if (m_useTmpRefStateStorage) {
      /*
       * Call the function that makes sure the local copy of 
       * the species reference thermo functions are up to date
       * for the current temperature.
       */
      _updateRefStateThermo();
      /*
       * Copy the gibbs function into return vector.
       */
      copy(m_cp0_R.begin(), m_cp0_R.end(), cpr);
    } else if (m_useTmpStandardStateStorage) {
      _updateStandardStateThermo(m_p0);
      copy(m_cpss_R.begin(), m_cpss_R.end(), cpr);
    } else {
      err("getCp_R_ref() ERROR: not handled");
    }
  }
  
  /*
   *  Get the molar volumes of the species reference states at the current
   *  <I>T</I> and <I>P_ref</I> of the solution.
   *
   * units = m^3 / kmol
   */
 void VPStandardStateTP::getStandardVolumes_ref(doublereal *vol) const {
   if (m_useTmpStandardStateStorage) {
     _updateStandardStateThermo(m_p0);
     copy(m_Vss.begin(), m_Vss.end(), vol);
   } else {
     err("getStandardVolumes_ref() ERROR: not handled");
   }
 }

  /*
   * Perform initializations after all species have been
   * added.
   */
  void VPStandardStateTP::initThermo() {
    initLengths();
    ThermoPhase::initThermo();
  }
  
  /*
   * Initialize the internal lengths.
   *       (this is not a virtual function)
   */
  void VPStandardStateTP::initLengths() {
    m_kk = nSpecies();
    int leng = m_kk;
    /*
     * malloc the storage for this even if
     * m_useTmpRefStateStorage is set to false.
     * So many functions need that temporary storage anyway.
     * However, that variable is still used to see if the
     * storage is used to supply the complete picture.
     */
    m_h0_RT.resize(leng);
    m_g0_RT.resize(leng);
    m_cp0_R.resize(leng);
    m_s0_R.resize(leng);
    
    if (m_useTmpStandardStateStorage) {
      m_hss_RT.resize(leng);
      m_gss_RT.resize(leng);
      m_cpss_R.resize(leng);
      m_sss_R.resize(leng);
      m_Vss.resize(leng);
    }
  }

  /*
   *   Import and initialize a ThermoPhase object
   *
   * param phaseNode This object must be the phase node of a
   *             complete XML tree
   *             description of the phase, including all of the
   *             species data. In other words while "phase" must
   *             point to an XML phase object, it must have
   *             sibling nodes "speciesData" that describe
   *             the species in the phase.
   * param id   ID of the phase. If nonnull, a check is done
   *             to see if phaseNode is pointing to the phase
   *             with the correct id.
   *
   * This routine initializes the lengths in the current object and
   * then calls the parent routine.
   */
  void VPStandardStateTP::initThermoXML(XML_Node& phaseNode, std::string id) {
    VPStandardStateTP::initLengths();
    ThermoPhase::initThermoXML(phaseNode, id);
  }
  
  /*
   * void _updateRefStateThermo()            (protected, virtual, const)
   *
   * This function checks to see whether the temperature has changed and
   * thus the reference thermodynamics functions for all of the species
   * must be recalculated.
   * It must be called for every reference state function evaluation,
   * if m_useTmpRefStateStorage is set to true.
   * If the temperature has changed, the species thermo manager is called
   * to recalculate the following internal arrays at the current temperature and at
   * the reference pressure:
   *
   *  - m_h0_RT
   *  - m_g0_RT
   *  - m_s0_R
   *  - m_cp0_R
   *
   * This function may be reimplemented in child objects. However, it doesn't
   * necessarily have to be, if the species thermo manager can carry
   * out the full calculation.
   */                    
  void VPStandardStateTP::_updateRefStateThermo() const {
    if (m_spthermo) {
      doublereal tnow = temperature();
      if (m_tlast_ref != tnow) {
	m_spthermo->update(tnow, DATA_PTR(m_cp0_R), DATA_PTR(m_h0_RT),
			   DATA_PTR(m_s0_R));
	m_tlast_ref = tnow;
	for (int k = 0; k < m_kk; k++) {
	  m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
	}
      }
    }
  }
  
  /*
   * void _updateStandardStateThermo()            (protected, virtual, const)
   *
   * If m_useTmpStandardStateStorage is true,
   * This function must be called for every call to functions in this
   * class that need standard state properties.
   * Child classes may require that it be called even if  m_useTmpStandardStateStorage
   * is not true.
   * It checks to see whether the temperature has changed and
   * thus the ss thermodynamics functions for all of the species
   * must be recalculated.
   *
   * This 
   */                    
  void VPStandardStateTP::_updateStandardStateThermo(doublereal pnow) const {
    _updateRefStateThermo();
    doublereal tnow = temperature();
    if (pnow == -1.0) {
      pnow = pressure();
    }
    if (m_tlast != tnow || m_plast != pnow) {
      err("_updateStandardStateThermo ERROR: Must be overwritten in child classes");
      /*
       * Redo objects that need reevaluation.
       */
      for (int k = 0; k < m_kk; k++) {
	m_g0_RT[k] = m_g0_RT[k];
      }
      m_tlast = tnow;
      m_plast = pnow;
    }
  }
}


