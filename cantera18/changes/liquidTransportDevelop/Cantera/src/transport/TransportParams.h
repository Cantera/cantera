#ifndef CT_TRANSPORTPARAMS_H
#define CT_TRANSPORTPARAMS_H

#include <vector>

#include "ct_defs.h"
#include "TransportBase.h"
#include "xml.h"
#include "XML_Writer.h"
#include "DenseMatrix.h"

namespace Cantera {

  /** 
   * Base class to hold transport model parameters.  
   * Used by TransportFactory.
   */
  class TransportParams {

  public:

    //! Default Constructor
    TransportParams() : 
      nsp_(0),
      thermo(0),
      mw(0),
      velocityBasis_(VB_MASSAVG),
      tmax(1000000.),
      tmin(10.),
      mode_(0),
      xml(0),
      log_level(-1)
    {
    }

    //! Destructor
    virtual ~TransportParams() 
    {
#ifdef DEBUG_MODE
      delete xml;
#endif
    }

    //! Local storage of the number of species
    int nsp_;

    //!  Pointer to the ThermoPhase object
    thermo_t* thermo;

    //! Local storage of the molecular weights of the species
    /*!
     *  Length is nsp_ and units are kg kmol-1.
     */
    vector_fp        mw;

    //! A basis for the average velocity can be specified. 
    /*!
     *  Valid bases include "mole", "mass", and "species" names.
     */
    VelocityBasis velocityBasis_;

    //! Maximum temperatures for parameter fits
    doublereal tmax;

    //! Minimum temperatures for parameter fits
    doublereal tmin;

    //!  Mode parameter
    int mode_;

    //! Pointer to the xml tree describing the implementation of transport for this object
    XML_Writer* xml;

    //! Log level
    int log_level;

  };


  /**
   * Holds transport model parameters relevant to transport in ideal 
   * gases with a kinetic theory of gases derived transport model. 
   * Used by TransportFactory.
   */
  class GasTransportParams : public TransportParams {

  public:

    GasTransportParams() :
      TransportParams(),
      visccoeffs(0),
      condcoeffs(0),
      diffcoeffs(0),
      polytempvec(0),
      poly(0),
      omega22_poly(0),
      astar_poly(0),
      bstar_poly(0),
      cstar_poly(0),
      zrot(0),
      crot(0),
      polar(0),
      alpha(0),
      fitlist(0),
      eps(0),
      sigma(0),
      reducedMass(0, 0),
      diam(0, 0),
      epsilon(0, 0),
      dipole(0, 0),
      delta(0, 0)
    {
    }

    virtual ~GasTransportParams() {}

    // polynomial fits

    //! temperature-fit of the viscosity 
    std::vector<vector_fp> visccoeffs; 

    //temperature-fit heat conduction 
    std::vector<vector_fp>            condcoeffs; 

    //temperature-fit diffusivity 
    std::vector<vector_fp>            diffcoeffs; 
    vector_fp                         polytempvec;

    std::vector<std::vector<int> > poly;
    std::vector<vector_fp >   omega22_poly;
    std::vector<vector_fp >   astar_poly;
    std::vector<vector_fp >   bstar_poly;
    std::vector<vector_fp >   cstar_poly;

    //! Rotational relaxation number for the species in the current phase
    /*!
     * length is the number of species in the phase
     * units are dimensionless
     */
    vector_fp   zrot;

    //! Dimensionless rotational heat capacity of the species in the current phase
    /*!
     *  These values are 0, 1 and 1.5 for single-molecule, linear, and nonlinear species respectively
     *  length is the number of species in the pahse
     *  units are dimensionless  (Cr / R)
     */
    vector_fp   crot;

    std::vector<bool> polar;
    vector_fp    alpha;
    vector_fp    fitlist;

    //! Lennard-Jones well-depth of the species in the current phase
    /*!
     * length is the number of species in the phase
     * Units are Joules (Note this is not Joules/kmol)
     */
    vector_fp    eps;

    //! Lennard-Jones diameter of the species in the current phase
    /*!
     * length is the number of species in the phase
     * units are in meters.
     */
    vector_fp    sigma;

    DenseMatrix  reducedMass;    
    DenseMatrix  diam;           
    DenseMatrix  epsilon;        
    DenseMatrix  dipole;         
    DenseMatrix  delta;          

  };

}

#endif //CT_TRANSPORTPARAMS_H
