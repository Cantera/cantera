/**
 *  @file LiquidTransportParams.h
 *  Header file defining class LiquidTransportParams
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_LIQUIDTRANSPORTPARAMS_H
#define CT_LIQUIDTRANSPORTPARAMS_H

#include "TransportParams.h"
#include "LiquidTranInteraction.h"

namespace Cantera
{

//! Class LiquidTransportParams holds transport model parameters relevant to
//! transport in mixtures.
/*!
 * @deprecated To be removed after Cantera 2.4
 *
 * This class is used by TransportFactory to initialize transport objects.
 */
class LiquidTransportParams : public TransportParams
{
public:
    LiquidTransportParams();
    ~LiquidTransportParams();
    LiquidTransportParams(const LiquidTransportParams& right) = delete;
    LiquidTransportParams& operator=(const LiquidTransportParams& right) = delete;

    //! Species transport parameters
    std::vector<LiquidTransportData> LTData;

    //! Object that specifies the viscosity interaction for the mixture
    LiquidTranInteraction* viscosity;

    //! Object that specifies the ionic Conductivity of the mixture
    LiquidTranInteraction* ionConductivity;

    //! Vector of pointer to the LiquidTranInteraction object which handles the
    //! calculation of the mobility ratios for the phase
    /*!
     * The mobility ratio is defined via the following quantity where i and j
     * are species indices.
     *
     *   mobRat(i,j) = mu_i / mu_j
     *
     * It is returned in fortran-ordering format. i.e. it is returned as
     * mobRat[k], where
     *
     *        k = j * nsp + i
     *
     * Length = nsp * nsp
     */
    std::vector<LiquidTranInteraction*> mobilityRatio;

    //! Vector of pointer to the LiquidTranInteraction object which handles the
    //! calculation of each species' self diffusion coefficient for the phase
    std::vector<LiquidTranInteraction*> selfDiffusion;

    //! Pointer to the LiquidTranInteraction object which handles the
    //! calculation of the mixture thermal conductivity for the phase
    LiquidTranInteraction* thermalCond;

    //! Pointer to the LiquidTranInteraction object which handles the
    //! calculation of the species diffusivity for the phase
    LiquidTranInteraction* speciesDiffusivity;

    //! Pointer to the LiquidTranInteraction object which handles the
    //! calculation of the electrical conductivity for the phase
    LiquidTranInteraction* electCond;

    //! Pointer to the LiquidTranInteraction object which handles the
    //! calculation of the hydrodynamic radius for the phase
    /*!
     *  @note  I don't understand at the moment how one can define a
     *         hydrodynamic radius for the phase
     */
    LiquidTranInteraction* hydroRadius;

    //! Model for species interaction effects for viscosity
    //! Takes enum LiquidTranMixingModel
    LiquidTranMixingModel model_viscosity;

    //! Model for species interaction effects for ionic conductivity
    //! Takes enum LiquidTranMixingModel
    LiquidTranMixingModel model_ionConductivity;

    //! Model for species interaction effects for mobility ratio
    //! Takes enum LiquidTranMixingModel
    std::vector<LiquidTranMixingModel*> model_mobilityRatio;

    //! Model for species interaction effects for mobility ratio
    //! Takes enum LiquidTranMixingModel
    std::vector<LiquidTranMixingModel*> model_selfDiffusion;

    //! Interaction associated with linear weighting of
    //! thermal conductivity.
    /**
     * This is used for either LTI_MODEL_MASSFRACS or LTI_MODEL_MOLEFRACS. The
     * overall formula for the mixture viscosity is
     *
     * \f[ \eta_{mix} = \sum_i X_i \eta_i + \sum_i \sum_j X_i X_j A_{i,j} \f].
     */
    DenseMatrix thermalCond_Aij;

    //! Model for species interaction effects for mass diffusivity
    //! Takes enum LiquidTranMixingModel
    LiquidTranMixingModel model_speciesDiffusivity;

    //! Interaction associated with linear weighting of
    //! thermal conductivity.
    /**
     * This is used for either LTI_MODEL_PAIRWISE_INTERACTION or
     * LTI_MODEL_STEFANMAXWELL_PPN. These provide species interaction
     * coefficients associated with the Stefan-Maxwell formulation.
     */
    DenseMatrix diff_Dij;

    //! Model for species interaction effects for hydrodynamic radius
    //! Takes enum LiquidTranMixingModel
    LiquidTranMixingModel model_hydroradius;

    //! Interaction associated with hydrodynamic radius.
    /**
     * Not yet implemented
     */
    DenseMatrix radius_Aij;

    //! Default composition dependence of the transport properties
    /*!
     *   Permissible types of composition dependencies
     *    0 - Solvent values (i.e., species 0) contributes only
     *    1 - linear combination of mole fractions;
     */
    LiquidTranMixingModel compositionDepTypeDefault_;
};

}

#endif
