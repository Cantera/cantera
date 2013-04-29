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

#ifndef CT_TORTUOSITYBRUGGEMAN_H
#define CT_TORTUOSITYBRUGGEMAN_H

#include "TortuosityBase.h"

namespace Cantera
{

//!  Maxwell model for tortuosity
/*!
 *
 * This class implements transport coefficient corrections
 * appropriate for porous media with a dispersed phase.
 * This model goes back to Maxwell.  The formula for the
 * conductivity is expressed in terms of the volume fraction
 * of the continuous phase, \f$ \phi \f$, and the relative
 * conductivities of the dispersed and continuous phases,
 * \f$ r = \kappa_d / \kappa_0 \f$.  For dilute particle
 * suspensions the effective conductivity is
 *
 * \f[
 *    \kappa / \kappa_0 = 1 + 3 ( 1 - \phi ) ( r - 1 ) / ( r + 2 )
 *                        + O(\phi^2)
 * \f]
 *
 * The class is derived from the TortuosityBase class.
 *
 */
class TortuosityMaxwell : public TortuosityBase
{

public:
    //! Default constructor uses Maxwelln exponent of 1.5
    /*!
     *  @param setPower       Exponent in the Maxwell factor. The default is 1.5
     */
    TortuosityMaxwell(double relativeConductivites = 0.0);

    //! Copy Constructor
    /*!
     * @param right  Object to be copied
     */
    TortuosityMaxwell(const TortuosityMaxwell& right);

    //! Assignment operator
    /*!
     * @param right Object to be copied
     */
    TortuosityMaxwell& operator=(const TortuosityMaxwell& right);

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

    //! Relative  conductivities of the dispersed and continuous phases,
    /*!
     *
     *   \f[
     *         \code{relativeConductivites_} = \kappa_d / \kappa_0
     *   \f]
     */
    doublereal relativeConductivities_;

};



}

#endif

