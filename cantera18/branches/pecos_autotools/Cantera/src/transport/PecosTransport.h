/**
 *  @file PecosTransport.h
 *   Header file defining class PecosTransport
 */

/* $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_PECOSTRAN_H
#define CT_PECOSTRAN_H


// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

// STL includes
#include <vector>
#include <string>
#include <map>
#include <numeric>
#include <algorithm>

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

using namespace std;

// Cantera includes
#include "TransportBase.h"
#include "DenseMatrix.h"

namespace Cantera {


  class GasTransportParams;

  /**
   *
   * Class PecosTransport implements mixture-averaged transport
   * properties for ideal gas mixtures. 
   *
   */
  class PecosTransport : public Transport {

  public:

    virtual ~PecosTransport() {}

    virtual int model() const { return cPecosTransport; }

    //! Viscosity of the mixture
    /*!
     *
     */
    virtual doublereal viscosity();

    virtual void getSpeciesViscosities(doublereal* visc)
    { update_T();  updateViscosity_T(); copy(m_visc.begin(), m_visc.end(), visc); }

    //! Return the thermal diffusion coefficients
    /*!
     * For this approximation, these are all zero.
     */
    virtual void getThermalDiffCoeffs(doublereal* const dt);

    /*! returns the mixture thermal conductivity
     *
     * This is computed using the lumped model,
     * \f[
     *    k = k^{tr} + k^{ve} 
     * \f]
     * where, 
     * \f[
     *    k^{tr}= 5/2 \mu_s C_{v,s}^{trans} + \mu_s C_{v,s}^{rot}
     * \f]
     * and,
     * \f[
     *    k^{ve}= \mu_s C_{v,s}^{vib} + \mu_s C_{v,s}^{elec}
     * \f]
     * 
     */
    virtual doublereal thermalConductivity();

    virtual void getBinaryDiffCoeffs(const int ld, doublereal* const d);

    
    //! Mixture-averaged diffusion coefficients [m^2/s]. 
    /*!
    * For the single species case or the pure fluid case
    * the routine returns the self-diffusion coefficient.
    * This is need to avoid a Nan result in the formula
    * below.
    */
    virtual void getMixDiffCoeffs(doublereal* const d);
    virtual void getMobilities(doublereal* const mobil);
    virtual void update_T();
    virtual void update_C();

    //! Get the species diffusive mass fluxes wrt to 
    //! the mass averaged velocity, 
    //! given the gradients in mole fraction and temperature
    /*!
     *  Units for the returned fluxes are kg m-2 s-1.
     * 
     *  @param ndim Number of dimensions in the flux expressions
     *  @param grad_T Gradient of the temperature
     *                 (length = ndim)
     * @param ldx  Leading dimension of the grad_X array 
     *              (usually equal to m_nsp but not always)
     * @param grad_X Gradients of the mole fraction
     *             Flat vector with the m_nsp in the inner loop.
     *             length = ldx * ndim
     * @param ldf  Leading dimension of the fluxes array 
     *              (usually equal to m_nsp but not always)
     * @param fluxes  Output of the diffusive mass fluxes
     *             Flat vector with the m_nsp in the inner loop.
     *             length = ldx * ndim
     */
    virtual void getSpeciesFluxes(int ndim, 
				  const doublereal* grad_T,
				  int ldx,
				  const doublereal* grad_X, 
				  int ldf, doublereal* fluxes);

    //! Initialize the transport object
    /*!
     *
     * Here we change all of the internal dimensions to be sufficient.
     * We get the object ready to do property evaluations.
     *
     * @param tr  Transport parameters for all of the species
     *            in the phase.
     * 
     */
    virtual bool initGas( GasTransportParams& tr );


    /**         
     *
     * Reads the transport table specified (currently defaults to internal file)
     * 
     * Reads the user-specified transport table, appending new species                            
     * data and/or replacing default species data.   
     *
     */
    void read_blottner_transport_table ();

    friend class TransportFactory;

    /**
     * Return a structure containing all of the pertinent parameters
     * about a species that was used to construct the Transport
     * properties in this object.
     *
     * @param k Species number to obtain the properties from.
     */
    struct GasTransportData getGasTransportData(int);

  protected:

    /// default constructor
    PecosTransport();

  private:

    //! Calculate the pressure from the ideal gas law
    doublereal pressure_ig() const {
      return (m_thermo->molarDensity() * GasConstant *
	      m_thermo->temperature());
    }

    // mixture attributes
    int m_nsp;
    doublereal m_tmin, m_tmax;
    vector_fp  m_mw;

    // polynomial fits
    vector<vector_fp>            m_visccoeffs;
    vector<vector_fp>            m_condcoeffs;
    vector<vector_fp>            m_diffcoeffs;
    vector_fp                    m_polytempvec;

    // blottner fits
    //int species = 20;
    double a[20], b[20], c[20];

    // property values
    DenseMatrix                  m_bdiff;
    vector_fp                    m_visc;
    vector_fp                    m_sqvisc;
    vector_fp                    m_cond;

    array_fp                    m_molefracs;

    vector<vector<int> > m_poly;
    vector<vector_fp >   m_astar_poly;
    vector<vector_fp >   m_bstar_poly;
    vector<vector_fp >   m_cstar_poly;
    vector<vector_fp >   m_om22_poly;
    DenseMatrix          m_astar;
    DenseMatrix          m_bstar;
    DenseMatrix          m_cstar;
    DenseMatrix          m_om22;

    DenseMatrix m_phi;            // viscosity weighting functions
    DenseMatrix m_wratjk, m_wratkj1;

    vector_fp   m_zrot;
    vector_fp   m_crot;
    vector_fp   m_cinternal;
    vector_fp   m_eps;
    vector_fp   m_alpha;
    vector_fp   m_dipoleDiag;

    doublereal m_temp, m_logt, m_kbt, m_t14, m_t32;
    doublereal m_sqrt_kbt, m_sqrt_t;

    vector_fp  m_sqrt_eps_k;
    DenseMatrix m_log_eps_k;
    vector_fp  m_frot_298;
    vector_fp  m_rotrelax;

    doublereal m_lambda;
    doublereal m_viscmix;

    // work space
    vector_fp  m_spwork;

    void updateThermal_T();
    void updateViscosity_T();
    void updateCond_T();
    void updateSpeciesViscosities();
    void updateDiff_T();
    void correctBinDiffCoeffs();
    bool m_viscmix_ok;
    bool m_viscwt_ok;
    bool m_spvisc_ok;
    bool m_diffmix_ok;
    bool m_bindiff_ok;
    bool m_abc_ok;
    bool m_spcond_ok;
    bool m_condmix_ok;

    int m_mode;

    DenseMatrix m_epsilon;
    DenseMatrix m_diam;
    DenseMatrix incl;
    bool m_debug;

    // specific heats
    vector_fp            cv_rot;
    vector_fp            cp_R;
    vector_fp            cv_int;

  };
}
#endif
