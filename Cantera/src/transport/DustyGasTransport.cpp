/**
 *
 *  @file DustyGasTransport.cpp
 *  Implementation file for class DustyGasTransport
 *
 *  @ingroup transportProps
 *
 */

/*
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2003 California Institute of Technology
 *  See file License.txt for licensing information
 *
 */


// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ThermoPhase.h"
#include "DustyGasTransport.h"

using namespace std;

/** 
 * Mole fractions below MIN_X will be set to MIN_X when computing
 * transport properties.
 */
#define MIN_X 1.e-20


namespace Cantera {

  //====================================================================================================================
  DustyGasTransport::DustyGasTransport(thermo_t* thermo) :
      Transport(thermo),
      m_nsp(0),
      m_mw(0),
      m_dk(0),
      m_temp(-1.0),     
      m_multidiff(0,0),
      m_spwork(0),
      m_spwork2(0),
      m_gradP(0.0),
      m_knudsen_ok(false),
      m_bulk_ok(false),
      m_gradConc_set(false),
      m_gradP_set(false),
      m_porosity(0.0),
      m_tortuosity(1.0),
      m_pore_radius(0.0),
      m_diam(0.0),
      m_perm(-1.0),
      m_gastran(0)
  {
  }
 //====================================================================================================================
  DustyGasTransport::DustyGasTransport(const DustyGasTransport &right) :
      Transport(),
      m_nsp(0),
      m_mw(0),
      m_dk(0),
      m_temp(-1.0),     
      m_multidiff(0,0),
      m_spwork(0),
      m_spwork2(0),
      m_gradP(0.0),
      m_knudsen_ok(false),
      m_bulk_ok(false),
      m_gradConc_set(false),
      m_gradP_set(false),
      m_porosity(0.0),
      m_tortuosity(1.0),
      m_pore_radius(0.0),
      m_diam(0.0),
      m_perm(-1.0),
      m_gastran(0)
  {
    *this = right;
  }
  //====================================================================================================================
  // Assignment operator
  /*
   *  This is NOT a virtual function.
   *
   * @param right    Reference to %DustyGasTransport object to be copied 
   *                 into the current one.
   */
  DustyGasTransport& DustyGasTransport::operator=(const  DustyGasTransport& right)
  {
    if (&right == this) {
      return *this;
    }
    Transport::operator=(right);

    m_nsp = right.m_nsp;
    m_mw = right.m_mw;
    m_d = right.m_d;
    m_x = right.m_x;
    m_dk = right.m_dk;
    m_temp = m_temp;
    m_multidiff = right.m_multidiff;
    m_spwork = right.m_spwork;
    m_spwork2 = right.m_spwork2;
    m_gradP = right.m_gradP;
    m_knudsen_ok = right.m_knudsen_ok;
    m_bulk_ok= right.m_bulk_ok;
    m_gradConc_set = right.m_gradConc_set;
    m_gradP_set = right.m_gradP_set;
    m_porosity = right.m_porosity;  
    m_tortuosity = right.m_tortuosity;   
    m_pore_radius = right.m_pore_radius;  
    m_diam = right.m_diam;
    m_perm = right.m_perm;

    // Warning -> This is a shallow pointer copy. gastran may not point to the correct object
    //            after this copy. The routine initialize() must be called 
    m_gastran = right.m_gastran;

    return *this;
  }
  //====================================================================================================================
  DustyGasTransport::~DustyGasTransport() {
  }
  //====================================================================================================================
  // Duplication routine for objects which inherit from %Transport
  /*
   *  This virtual routine can be used to duplicate %Transport objects
   *  inherited from %Transport even if the application only has
   *  a pointer to %Transport to work with.
   *
   *  These routines are basically wrappers around the derived copy
   *  constructor.
   */
  Transport *DustyGasTransport::duplMyselfAsTransport() const {
    DustyGasTransport* tr = new DustyGasTransport(*this);
    return (dynamic_cast<Transport *>(tr));
  }
  //====================================================================================================================
  void DustyGasTransport::setParameters(const int type, const int k, const doublereal* const p) {
    switch(type) {
    case 0:
      setPorosity(p[0]); break;
    case 1:
      setTortuosity(p[0]); break;
    case 2:
      setMeanPoreRadius(p[0]); break;
    case 3:
      setMeanParticleDiameter(p[0]); break;
    case 4:
      setPermeability(p[0]); break;
    default:
      throw CanteraError("DustyGasTransport::init",
			 "unknown parameter");
    }
  }
  //====================================================================================================================
  //  Initialization routine called by TransportFactory
  /*
   *  The DustyGas model is a subordinate model to the gas phase transport model. Here we 
   *  set the gas phase models.
   *
   *  This is a protected routine, so that initialiation of the Model must occur within Cantera's setup
   *
   *   @param  phase           Pointer to the underlying ThermoPhase model for the gas phase
   *   @param  gastr           Pointer to the underlying Transport model for transport in the gas phse.
   */
  void DustyGasTransport::initialize(ThermoPhase* phase, Transport* gastr) {

    // constant mixture attributes
    m_thermo = phase;
    m_nsp   = m_thermo->nSpecies();
    m_gastran = gastr;

    // make a local copy of the molecular weights
    m_mw.resize(m_nsp);
    copy(m_thermo->molecularWeights().begin(),  m_thermo->molecularWeights().end(), m_mw.begin());

    m_multidiff.resize(m_nsp, m_nsp);
    m_d.resize(m_nsp, m_nsp);
    m_dk.resize(m_nsp, 0.0);

    m_x.resize(m_nsp, 0.0);
    m_thermo->getMoleFractions(DATA_PTR(m_x));

    // set flags all false
    m_knudsen_ok = false;
    m_bulk_ok = false;
    m_gradConc_set = false;
    m_gradP_set = false;
  
    m_spwork.resize(m_nsp);
    m_spwork2.resize(m_nsp);
  }
  //====================================================================================================================

  void DustyGasTransport::updateBinaryDiffCoeffs() {
    if (m_bulk_ok) return;
    int n,m;

    // get the gaseous binary diffusion coefficients
    m_gastran->getBinaryDiffCoeffs(m_nsp, m_d.ptrColumn(0));
    doublereal por2tort = m_porosity / m_tortuosity;
    for (n = 0; n < m_nsp; n++) {
      for (m = 0; m < m_nsp; m++) {
	m_d(n,m) *= por2tort;
      }
    }
    m_bulk_ok = true;
  }
  //====================================================================================================================
  void DustyGasTransport::updateKnudsenDiffCoeffs() {
    if (m_knudsen_ok) return;
    doublereal K_g = m_pore_radius * m_porosity / m_tortuosity;
    const doublereal TwoThirds = 2.0/3.0;
    for (int k = 0; k < m_nsp; k++) {
      m_dk[k] = TwoThirds * K_g * sqrt((8.0 * GasConstant * m_temp)/
				       (Pi * m_mw[k]));
    }
    m_knudsen_ok = true;
  }

  //====================================================================================================================      
  void DustyGasTransport::eval_H_matrix() {
    updateBinaryDiffCoeffs();
    updateKnudsenDiffCoeffs();
    int k,l,j;
    doublereal sum;
    for (k = 0; k < m_nsp; k++) {

      // evaluate off-diagonal terms
      for (l = 0; l < m_nsp; l++) {
	m_multidiff(k,l) = -m_x[k]/m_d(k,l);
      }

      // evaluate diagonal term
      sum = 0.0;
      for (j = 0; j < m_nsp; j++) {
	if (j != k) {
	  sum += m_x[j]/m_d(k,j);
	}
      }
      m_multidiff(k,k) = 1.0/m_dk[k] + sum;
    }
  }
  //====================================================================================================================
  void DustyGasTransport::getMolarFluxes(const doublereal* const state1,
					 const doublereal * const state2, 
					 const doublereal delta,
					 doublereal * const fluxes) {

    int k;
    doublereal conc1, conc2;
    doublereal* cbar = DATA_PTR(m_spwork);
    doublereal* gradc = DATA_PTR(m_spwork2);
    doublereal t1 = state1[0];
    doublereal t2 = state2[0];
    doublereal rho1 = state1[1];
    doublereal rho2 = state2[1];
    const doublereal* y1 = state1 + 2;
    const doublereal* y2 = state2 + 2;
    doublereal c1sum = 0.0, c2sum = 0.0;
    for (k = 0; k < m_nsp; k++) {
      conc1 = rho1*y1[k]/m_mw[k];
      conc2 = rho2*y2[k]/m_mw[k];
      cbar[k] = 0.5*(conc1 + conc2);
      gradc[k] = (conc2 - conc1)/delta;
      c1sum += conc1;
      c2sum += conc2;
    }
    doublereal p1 = c1sum * GasConstant * state1[0];
    doublereal p2 = c2sum * GasConstant * state2[0];
    doublereal pbar = 0.5*(p1 + p2);
    doublereal gradp = (p2 - p1)/delta;
    doublereal tbar = 0.5*(t1 + t2);

    m_thermo->setState_TPX(tbar, pbar, cbar);

    updateMultiDiffCoeffs();

    multiply(m_multidiff, gradc, fluxes);
    divide_each(cbar, cbar + m_nsp, m_dk.begin());

    // if no permeability has been specified, use result for 
    // close-packed spheres
    double b = 0.0;
    if (m_perm < 0.0) {
      double p = m_porosity;
      double d = m_diam;
      double t = m_tortuosity;
      b = p*p*p*d*d/(72.0*t*(1.0-p)*(1.0-p));
    }
    else {
      b = m_perm;
    }
    b *= gradp / m_gastran->viscosity();
    scale(cbar, cbar + m_nsp, cbar, b);
    increment(m_multidiff, cbar, fluxes);
    scale(fluxes, fluxes + m_nsp, fluxes, -1.0);
  }
  //====================================================================================================================

  void DustyGasTransport::updateMultiDiffCoeffs() {
    // see if temperature has changed
    updateTransport_T();

    // update the mole fractions
    updateTransport_C();

    eval_H_matrix();

    // invert H
    int ierr = invert(m_multidiff);

    if (ierr != 0) {
      throw CanteraError("DustyGasTransport::updateMultiDiffCoeffs",
			 "invert returned ierr = "+int2str(ierr));
    }
  }
  //====================================================================================================================
  // Return the Multicomponent diffusion coefficients. Units: [m^2/s]. 
  /*
   * Returns the array of multicomponent diffusion coefficients.
   *
   *  @param ld  The dimension of the inner loop of d (usually equal to m_nsp)
   *  @param d  flat vector of diffusion coefficients, fortran ordering.
   *            d[ld*j+i] is the D_ij diffusion coefficient (the diffusion
   *            coefficient for species i due to species j).
   */
  void DustyGasTransport::getMultiDiffCoeffs(const int ld, doublereal* const d) {
    int i,j;
    updateMultiDiffCoeffs();
    for (i = 0; i < m_nsp; i++) {
      for (j = 0; j < m_nsp; j++) {
	d[ld*j + i] = m_multidiff(i,j);
      }
    }
  }
  //====================================================================================================================    
  // Update temperature-dependent quantities within the object
  /*
   *  The object keeps a value m_temp, which is the temperature at which quantities were last evaluated
   *  at. If the temperature is changed, update Booleans are set false, triggering recomputation.
   */
  void DustyGasTransport::updateTransport_T() 
  {
    if (m_temp == m_thermo->temperature()) return;
    m_temp = m_thermo->temperature();
    m_knudsen_ok = false;
    m_bulk_ok = false;
  }                 
  //====================================================================================================================
  void DustyGasTransport::updateTransport_C()  
  {
    m_thermo->getMoleFractions(DATA_PTR(m_x));

    // add an offset to avoid a pure species condition
    // (check - this may be unnecessary)
    for (int k = 0; k < m_nsp; k++) {
      m_x[k] = fmaxx(MIN_X, m_x[k]);
    }
    // diffusion coeffs depend on Pressure
    m_bulk_ok = false;
  } 
  //====================================================================================================================
}
