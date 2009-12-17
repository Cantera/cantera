/**
 *  @file SolidTransport.h
 *   Header file for defining the class SolidTransport, which handles transport
 *   of ions within solid phases
 *  (see \ref tranprops and \link Cantera::SolidTransport SolidTransport \endlink).
 */

/* 
 * $Revision$
 * $Date$
 */

// Copyright 2003  California Institute of Technology


#ifndef CT_SOLIDTRAN_H
#define CT_SOLIDTRAN_H


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


// Cantera includes
#include "TransportBase.h"
#include "DenseMatrix.h"

namespace Cantera {


  
  //! Class SolidTransport implements transport properties for solids.
  /*!
   *
   *
   *
   */
  class SolidTransport : public Transport {

  public:
    
    //! Default constructor
    SolidTransport();

    //! Copy Constructor
    /*!
     *  @param right  Object to be copied
     */
    SolidTransport(const SolidTransport &right);

    //! Destructor
    virtual ~SolidTransport();

    //! Assignment operator
    /*!
     *  This is NOT a virtual function.
     *
     * @param right    Reference to Transport object to be copied into the
     *                 current one.
     */
    SolidTransport&  operator=(const SolidTransport& right);

    //! Duplication routine for objects which inherit from
    //! %Transport
    /*!
     *  This virtual routine can be used to duplicate %Transport objects
     *  inherited from %Transport even if the application only has
     *  a pointer to %Transport to work with.
     *
     *  These routines are basically wrappers around the derived copy
     *  constructor.
     */
    virtual Transport *duplMyselfAsTransport() const;


    virtual int model() const { return cSolidTransport; }

    virtual doublereal thermalConductivity();
    virtual void getMixDiffCoeffs(doublereal* const d);

    //!  Compute the electrical mobilities of the species from the diffusion coefficients, 
    //!  using the Einstein relation.
    /*!
     *   Frequently, but not always, the mobility is calculated from the
     *   diffusion coefficient using the Einstein relation
     *
     *     \f[
     *          \mu^e_k = \frac{F D_k}{R T}
     *     \f]
     *
     *  units (m^2/V/s).
     *  @param mobil   Returns the mobilities of
     *                 the species in array \c mobil_e. The array must be
     *                 dimensioned at least as large as the number of species.
     */
    virtual void getMobilities(doublereal* const mobil);

    virtual void setParameters(const int n, const int k, const doublereal* const p);

    friend class TransportFactory;

    /**
     * The electrical conductivity (Siemens/m).
     */
    virtual doublereal electricalConductivity();
  

  private:

    int m_nmobile;    // number of mobile species
    vector_fp m_Adiff;
    vector_fp m_Ndiff;
    vector_fp m_Ediff;
    vector_int m_sp;
    doublereal m_Alam;
    doublereal m_Nlam;
    doublereal m_Elam;
    vector_fp m_work;
  };
}
#endif






