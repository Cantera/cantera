/**
 *  @file LiquidTransportParams.h
 *  Header file defining class LiquidTransportParams
 */

#ifndef CT_LIQUIDTRANSPORTPARAMS_H
#define CT_LIQUIDTRANSPORTPARAMS_H

#include "cantera/base/ct_defs.h"
#include "TransportBase.h"
#include "TransportParams.h"
#include "LiquidTransportData.h"
#include "LiquidTranInteraction.h"
#include "cantera/base/xml.h"
#include "cantera/base/XML_Writer.h"

#include <vector>

namespace Cantera
{
//! Composition dependence type for liquid mixture transport properties
/*!
 *  Types of temperature dependencies:
 *  -   0  - Mixture calculations with this property are not allowed
 *  -   1  - Use solvent (species 0) properties
 *  -   2  - Properties weighted linearly by mole fractions
 *  -   3  - Properties weighted linearly by mass fractions
 *  -   4  - Properties weighted logarithmically by mole fractions (interaction energy weighting)
 *  -   5  - Interactions given pairwise between each possible species (i.e. D_ij)
 *
 *   \verbatim
 *    <transport model="Liquid">
 *       <viscosity>
 *          <compositionDependence model="logMoleFractions">
 *             <interaction>
 *                <speciesA> LiCl(L) </speciesA>
 *                <speciesB> KCl(L)  </speciesB>
 *                <Eij units="J/kmol"> -1.0 </Eij>
 *                <Sij units="J/kmol/K"> 1.0E-1 </Sij>
 *     -or-       <Sij>
 *                  <floatArray units="J/kmol/K"> 1.0E-1, 0.001 0.01 </floatArray>
 *                </Sij>
 *     -same form for Hij,Aij,Bij-
 *             </interaction>
 *          </compositionDependence>
 *       </viscosity>
 *       <speciesDiffusivity>
 *          <compositionDependence model="pairwiseInteraction">
 *             <interaction>
 *                <speciesA> Li+ </speciesA>
 *                <speciesB> K+  </speciesB>
 *                <Dij units="m2/s"> 1.5 </Dij>
 *             </interaction>
 *             <interaction>
 *                <speciesA> K+  </speciesA>
 *                <speciesB> Cl- </speciesB>
 *                <Dij units="m2/s"> 1.0 </Dij>
 *             </interaction>
 *             <interaction>
 *                <speciesA> Li+  </speciesA>
 *                <speciesB> Cl-  </speciesB>
 *                <Dij units="m2/s"> 1.2 </Dij>
 *             </interaction>
 *          </compositionDependence>
 *       </speciesDiffusivity>
 *       <thermalConductivity>
 *          <compositionDependence model="massFractions"/>
 *       </thermalConductivity>
 *       <hydrodynamicRadius>
 *          <compositionDependence model="none"/>
 *       </hydrodynamicRadius>
 *    </transport>
 *   \endverbatim
 *
 */




//! Class LiquidTransportParams holds transport model parameters
//!  relevant to transport in mixtures.
/*!
 * This class is used by TransportFactory to initialize transport objects.
 */
class LiquidTransportParams : public TransportParams
{

public:

    //! Constructor
    LiquidTransportParams();

    //! Destructor
    ~LiquidTransportParams();

    //! Copy constructor
    /*!
     * @param right  Object to be copied
     */
    LiquidTransportParams(const LiquidTransportParams& right);

    //! Assignment operator
    /*!
     * @param right  Object to be copied
     */
    LiquidTransportParams& operator=(const LiquidTransportParams& right);

    //! Species transport parameters
    std::vector<Cantera::LiquidTransportData> LTData;

    //! Object that specifies the viscosity interaction for the mixture
    LiquidTranInteraction* viscosity;

    //! Object that specifes the ionic Conductivity of the mixture
    LiquidTranInteraction* ionConductivity;

    //! Vector of pointer to the LiquidTranInteraction object which handles the calculation of
    //! each species' mobility ratios for the phase
    /*!
     *   mobRat(i,j) = mu_i / mu_j
     *
     *    It is returned in fortran-ordering format. ie. it is returned as mobRat[k], where
     *
     *        k = j * nsp + i
     */
    std::vector<LiquidTranInteraction*> mobilityRatio;

    //! Vector of pointer to the LiquidTranInteraction object which handles the calculation of
    //! each species' self diffusion coefficient for the phase
    std::vector<LiquidTranInteraction*> selfDiffusion;

    //! Pointer to the  LiquidTranInteraction object which handles the calculation of the
    //! mixture thermal conductivity for the phase
    LiquidTranInteraction* thermalCond;

    //! Pointer to the  LiquidTranInteraction object which handles the calculation of the
    //! species diffusivity for the phase
    LiquidTranInteraction* speciesDiffusivity;

    //! Pointer to the  LiquidTranInteraction object which handles the calculation of the
    //! electrical conductivity for the phase
    LiquidTranInteraction* electCond;

    //! Pointer to the  LiquidTranInteraction object which handles the calculation of the hydrodynamic
    //! radius for the phase
    /*!
     *  @note  I don't understand at the moment how one can define a hydrodynamic
     *         radius for the phase
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
     * This is used for either LTI_MODEL_MASSFRACS
     * or LTI_MODEL_MOLEFRACS.
     * The overall formula for the mixture viscosity is
     *
     * \f[ \eta_{mix} = \sum_i X_i \eta_i
     *  + \sum_i \sum_j X_i X_j A_{i,j} \f].
     */
    DenseMatrix  thermalCond_Aij;

    //! Model for species interaction effects for mass diffusivity
    //! Takes enum LiquidTranMixingModel
    LiquidTranMixingModel model_speciesDiffusivity;

    //! Interaction associated with linear weighting of
    //! thermal conductivity.
    /**
     * This is used for either LTI_MODEL_PAIRWISE_INTERACTION or LTI_MODEL_STEFANMAXWELL_PPN.
     * These provide species interaction coefficients associated with
     * the Stefan-Maxwell formulation.
     */
    DenseMatrix  diff_Dij;

    //! Model for species interaction effects for hydrodynamic radius
    //! Takes enum LiquidTranMixingModel
    LiquidTranMixingModel model_hydroradius;

    //! Interaction associated with hydrodynamic radius.
    /**
     * Not yet implemented
     */
    DenseMatrix  radius_Aij;
};




}

#endif
