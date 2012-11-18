/**
 * @file TortuosityBase.h
 * Virtual base class to compute the increase in diffusive path length associated with
 * tortuous path diffusion through, for example, porous media.
 */

/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef CT_TORTUOSITYBASE_H
#define CT_TORTUOSITYBASE_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

//! Base case to handle tortuosity corrections for diffusive transport
//! in porous media
/*!
 * Class to compute the increase in diffusive path length associated with
 * tortuous path diffusion through, for example, porous media.
 * This base class implementation relates tortuosity to volume fraction
 * through a power-law relationship that goes back to Bruggeman.  The
 * exponent is referred to as the Bruggeman exponent.
 *
 * Note that the total diffusional flux is generally written as
 *
 * \f[
 *   \frac{ \phi C_T D_i \nabla X_i }{ \tau^2 }
 * \f]
 *
 * where \f$ \phi \f$ is the volume fraction of the transported phase,
 * \f$ \tau \f$ is referred to as the tortuosity.  (Other variables are
 * \f$ C_T \f$, the total concentration, \f$ D_i \f$, the diffusion
 * coefficient, and \f$ X_i \f$, the mole fraction with Fickian
 * transport assumed.)
 *
 * The tortuosity comes into play in conjunction the the
 */
class TortuosityBase
{

public:
    //! Default constructor uses Bruggeman exponent of 1.5
    TortuosityBase();

    //! Copy Constructor
    /*!
     * @param right  Object to be copied
     */
    TortuosityBase(const TortuosityBase& right);

    //! Default destructor for TortuosityBase
    virtual ~TortuosityBase();

    //! Assignment operator
    /*!
     * @param right Object to be copied
     */
    TortuosityBase& operator=(const TortuosityBase& right);

    //! Duplication operator
    /*!
     *  @return  Returns a pointer to a duplicate of the current object given a
     *           base class pointer
     */
    virtual TortuosityBase* duplMyselfAsTortuosityBase() const;

    //! The tortuosity factor models the effective increase in the
    //! diffusive transport length.
    /*!
     * This method returns \f$ 1/\tau^2 \f$ in the description of the  flux
     *
     *    \f$  C_T D_i \nabla X_i / \tau^2 \f$.
     *
     *
     */
    virtual doublereal tortuosityFactor(doublereal porosity);

    //! The McMillan number is the ratio of the flux-like
    //! variable to the value it would have without porous flow.
    /**
     * The McMillan number combines the effect of tortuosity
     * and volume fraction of the transported phase.  The net flux
     * observed is then the product of the McMillan number and the
     * non-porous transport rate.  For a conductivity in a non-porous
     * media, \f$ \kappa_0 \f$, the conductivity in the porous media
     * would be \f$ \kappa = (\rm McMillan) \kappa_0 \f$.
     */
    virtual doublereal McMillanFactor(doublereal porosity);

protected:

};



}

#endif

