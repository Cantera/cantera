/**
 *  @file SolidTransport.cpp
 *   Definition file for the class SolidTransport, which handles transport
 *   of ions within solid phases
 *  (see \ref tranprops and \link Cantera::SolidTransport SolidTransport \endlink).
 */
/* $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2008 California Institute of Technology


// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ThermoPhase.h"
#include "SolidTransportData.h"
#include "SolidTransport.h"

#include "utilities.h"
#include <iostream>

using namespace std;

namespace Cantera {

  //====================================================================================================================
  SolidTransport::SolidTransport() :
    Transport() ,
    m_nmobile(0),
    m_Adiff(0),
    m_Ndiff(0),
    m_Ediff(0),
    m_sp(0),
    m_Alam(-1.0),
    m_Nlam(0),
    m_Elam(0)
  {
  }
  //====================================================================================================================
  SolidTransport::~SolidTransport() 
  {
  } 
  //====================================================================================================================
  SolidTransport::SolidTransport(const SolidTransport &right) :
    Transport(),
    m_nmobile(0),
    m_Adiff(0),
    m_Ndiff(0),
    m_Ediff(0),
    m_sp(0),
    m_Alam(-1.0),
    m_Nlam(0),
    m_Elam(0)
  {
    /*
     * Use the assignment operator to do the brunt
     * of the work for the copy construtor.
     */
    *this = right;
  }
  //====================================================================================================================
  SolidTransport& SolidTransport::operator=(const SolidTransport& b) 
  {
    if (&b != this) {
      return *this; 
    }
    Transport::operator=(b);

    m_nmobile =  b.m_nmobile; 
    m_Adiff = b.m_Adiff;
    m_Ndiff = b.m_Ndiff;
    m_Ediff = b.m_Ediff;
    m_sp = b.m_sp;
    m_Alam = b.m_Alam;
    m_Nlam = b.m_Nlam;
    m_Elam = b.m_Elam;
   
    return *this;
  }
  //====================================================================================================================
  Transport *SolidTransport::duplMyselfAsTransport() const 
  {
    SolidTransport* tr = new SolidTransport(*this); 
    return (dynamic_cast<Transport *>(tr));
  }


  //====================================================================================================================
  // Initialize the transport object
  /*
   * Here we change all of the internal dimensions to be sufficient.
   * We get the object ready to do property evaluations.
   * A lot of the input required to do property evaluations is 
   * contained in the SolidTransportData class that is 
   * filled in TransportFactory. 
   *
   * @param tr  Transport parameters for the phase
   */
  bool SolidTransport::initSolid(SolidTransportData& tr) {

    m_thermo = tr.thermo;
    tr.thermo = 0;
    //m_nsp   = m_thermo->nSpecies();
    //m_tmin  = m_thermo->minTemp();
    //m_tmax  = m_thermo->maxTemp();

    // make a local copy of the molecular weights
    //m_mw.resize(m_nsp, 0.0);
    //copy(m_thermo->molecularWeights().begin(), 
    //     m_thermo->molecularWeights().end(), m_mw.begin());

    m_ionConductivity =  tr.ionConductivity;
    tr.ionConductivity = 0;
    m_electConductivity =  tr.electConductivity;
    tr.electConductivity = 0;
    m_thermalConductivity = tr.thermalConductivity;
    tr.thermalConductivity = 0;
    m_speciesDiffusivity = tr.speciesDiffusivity;
    tr.speciesDiffusivity = 0;

    return true;
  }




  //====================================================================================================================
  void SolidTransport::setParameters(const int n, const int k, const doublereal * const p) {
    switch (n) {
            
    case 0:
      // set the Arrhenius parameters for the diffusion coefficient
      // of species k.
      m_sp.push_back(k);
      m_Adiff.push_back(p[0]);
      m_Ndiff.push_back(p[1]);
      m_Ediff.push_back(p[2]);
      m_nmobile = m_sp.size();
      break;

    case 1:
      // set the thermal conductivity Arrhenius parameters.
      m_Alam = p[0];
      m_Nlam = p[2];            
      m_Elam = p[2];            
      break;

    default:
      ;
    }

    m_work.resize(m_thermo->nSpecies());
  }


  /******************  ionConductivity ******************************/

  // Returns the ionic conductivity of the phase
  /*
   *  The thermo phase needs to be updated (temperature) prior to calling this.
   *  The ionConductivity calculation is handled by subclasses of 
   *  LTPspecies as specified in the input file. 
   *   
   */ 
  doublereal SolidTransport::ionConductivity() {        
    // LTPspecies method
    return m_ionConductivity->getSpeciesTransProp();
  }



  /******************  electron Conductivity ******************************/

  // Returns the electron conductivity of the phase
  /*
   *  The thermo phase needs to be updated (temperature) prior to calling this.
   *  The ionConductivity calculation is handled by subclasses of 
   *  LTPspecies as specified in the input file. 
   *
   * There is also a legacy multicomponent diffusion approach to electrical conductivity.
   *   
   */ 
  doublereal SolidTransport::electricalConductivity() {        
    if ( m_nmobile == 0 ) {
      // LTPspecies method
      return m_electConductivity->getSpeciesTransProp();
    } else {
      getMobilities(&m_work[0]);
      int nsp = m_thermo->nSpecies();
      doublereal sum = 0.0;
      for (int k = 0; k < nsp; k++) {
	sum += m_thermo->charge(k) * m_thermo->moleFraction(k) * m_work[k];
      }
      return sum * m_thermo->molarDensity();
    } 
  }
  


  /******************  thermalConductivity ******************************/

  // Returns the thermal conductivity of the phase
  /*
   *  The thermo phase needs to be updated (temperature) prior to calling this.
   *  The thermalConductivity calculation is handled by subclasses of 
   *  LTPspecies as specified in the input file. 
   *   
   *  There is also a legacy method to evaluate 
   * \f[
   * \lambda = A T^n \exp(-E/RT)
   * \f]
   */ 
  doublereal SolidTransport::thermalConductivity() {        
    if ( m_Alam > 0.0 ) {
      //legacy test case?
      doublereal t = m_thermo->temperature();
      return m_Alam * pow(t, m_Nlam) * exp(-m_Elam/t);
    } else {
      // LTPspecies method
      return m_thermalConductivity->getSpeciesTransProp();
    } 
  }



  /******************  defectDiffusivity ******************************/

  // Returns the diffusivity of the phase
  /*
   *  The thermo phase needs to be updated (temperature) prior to calling this.
   *  The defectDiffusivity calculation is handled by subclasses of 
   *  LTPspecies as specified in the input file. 
   *   
   */ 
  doublereal SolidTransport::defectDiffusivity() {        
    // LTPspecies method
    return m_speciesDiffusivity->getSpeciesTransProp();
  }


  //====================================================================================================================
  /*
   * Compute the mobilities of the species from the diffusion coefficients, 
   * using the Einstein relation.
   */
  void SolidTransport::getMobilities(doublereal* const mobil) {
    int k;
    getMixDiffCoeffs(mobil);
    doublereal t = m_thermo->temperature();
    int nsp = m_thermo->nSpecies();
    doublereal c1 = ElectronCharge / (Boltzmann * t);
    for (k = 0; k < nsp; k++) {
      mobil[k] *= c1;
    }
  } 
  //====================================================================================================================
  /*
   * The diffusion coefficients are computed from 
   *
   * \f[
   * D_k = A_k T^{n_k} \exp(-E_k/RT).
   * \f]
   *
   * The diffusion coefficients are only non-zero for species for
   * which parameters have been specified using method
   * setParameters.
   */
  void SolidTransport::getMixDiffCoeffs(doublereal* const d) {
    doublereal temp = m_thermo->temperature();
    int nsp = m_thermo->nSpecies();
    int k;
    for (k = 0; k < nsp; k++) d[k] = 0.0;
    for (k = 0; k < m_nmobile; k++) {
      d[m_sp[k]] = 
	m_Adiff[k] * pow(temp, m_Ndiff[k]) * exp(-m_Ediff[k]/temp);
    }
  }
}
