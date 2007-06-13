/**
 *
 *  @file StoichSubstance.cpp
 *
 */

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ct_defs.h"
#include "mix_defs.h"
#include "StoichSubstance.h"
#include "SpeciesThermo.h"

namespace Cantera {


  // Default empty constructor
  StoichSubstance::StoichSubstance() :
    m_kk(0),
    m_tmin(0.0),
    m_tmax(0.0),
    m_press(OneAtm),
    m_p0(OneAtm),
    m_tlast(-1.0)  {
  }
  
  // Copy Constructor
  /*
   * Copy constructor for the object. Constructed
   * object will be a clone of this object, but will
   * also own all of its data.
   * This is a wrapper around the assignment operator
   *
   * @param right Object to be copied.
   */
  StoichSubstance::StoichSubstance(const StoichSubstance &right) :
	m_kk(0),
	m_tmin(0.0),
	m_tmax(0.0),
	m_press(OneAtm),
	m_p0(OneAtm),
	m_tlast(-1.0)  { 
    *this = operator=(right);
  }

  // Asignment operator
  /*
   * Assignment operator for the object. Constructed
   * object will be a clone of this object, but will
   * also own all of its data.
   *
   * @param right Object to be copied.
   */
  StoichSubstance& StoichSubstance::
  operator=(const StoichSubstance &right) {
    if (&right != this) {
      ThermoPhase::operator=(right);
      m_kk      = right.m_kk;
      m_tmin    = right.m_tmin;
      m_tmax    = right.m_tmax;
      m_press   = right.m_press;
      m_p0      = right.m_p0;
      m_tlast   = right.m_tlast;
      m_h0_RT   = right.m_h0_RT;
      m_cp0_R   = right.m_cp0_R;
      m_s0_R    = right.m_s0_R;
    }
    return *this;
  }

  // Duplicator from the %ThermoPhase parent class
  /*
   * Given a pointer to a %ThermoPhase object, this function will
   * duplicate the %ThermoPhase object and all underlying structures.
   * This is basically a wrapper around the copy constructor.
   *
   * @return returns a pointer to a %ThermoPhase
   */
  ThermoPhase *StoichSubstance::duplMyselfAsThermoPhase() const {
    ThermoPhase *igp = new StoichSubstance(*this);
    return (ThermoPhase *) igp;
  }

  // Destructor
  StoichSubstance::~StoichSubstance() {
  }

    void StoichSubstance::initThermo() {
        m_kk = nSpecies();
        if (m_kk > 1) {
            throw CanteraError("initThermo",
                "stoichiometric substances may only contain one species.");
        } 
        doublereal tmin = m_spthermo->minTemp();
        doublereal tmax = m_spthermo->maxTemp();
        if (tmin > 0.0) m_tmin = tmin;
        if (tmax > 0.0) m_tmax = tmax;
        m_p0 = refPressure();

        int leng = m_kk;
        m_h0_RT.resize(leng);
        m_cp0_R.resize(leng);
        m_s0_R.resize(leng);
    }


    void StoichSubstance::_updateThermo() const {
        doublereal tnow = temperature();
        if (m_tlast != tnow) {
            m_spthermo->update(tnow, &m_cp0_R[0], &m_h0_RT[0], 
                &m_s0_R[0]);
            m_tlast = tnow;
        }
    }

    void StoichSubstance::
    getUnitsStandardConc(double *uA, int k, int sizeUA) const {
	for (int i = 0; i < sizeUA; i++) {
	  uA[i] = 0.0;
	}
    }

    void StoichSubstance::setParameters(int n, double * c) {
        double rho = c[0];
        setDensity(rho);
    }

    void StoichSubstance::getParameters(int &n, double * const c) const {
        double rho = density();
        c[0] = rho;
    }

    void StoichSubstance::setParametersFromXML(const XML_Node& eosdata) {
        eosdata._require("model","StoichSubstance");
        doublereal rho = getFloat(eosdata, "density", "-");
        setDensity(rho);
    }

}




