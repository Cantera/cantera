#ifndef CT_LIQUIDTRANSPORTPARAMS_H
#define CT_LIQUIDTRANSPORTPARAMS_H

#include <vector>

#include "ct_defs.h"
#include "TransportBase.h"
#include "TransportParams.h"
#include "LiquidTransportData.h"
#include "xml.h"
#include "XML_Writer.h"

namespace Cantera {


    //! Composition dependence type for liquid mixture transport properties
    /*!
     *  Types of temperature dependencies:
     *     0  - Use solvent (species 0) properties
     *     1  - Properties weighted linearly by mole fractions
     *     2  - Properties weighted linearly by mass fractions
     *     3  - Properties weighted logarithmically by mole fractions (interaction energy weighting)
     *     4  - Interactions given pairwise between each possible species (i.e. D_ij)
     */
  enum LiquidTranMixingModel {
    LTR_MIXMODEL_NOTSET=-1,
    LTR_MIXMODEL_SOLVENT, 
    LTR_MIXMODEL_MOLEFRACS,
    LTR_MIXMODEL_MASSFRACS,
    LTR_MIXMODEL_LOG_MOLEFRACS,
    LTR_PAIRWISE_INTERACTIONS
  };
  

    /**
     * Holds transport model parameters relevant to transport in 
     * liquids for which activated jump processes limit transport
     * (giving Arrhenius type transport properties). 
     * Used by TransportFactory.
     */
    class LiquidTransportParams :public TransportParams {

    public:

        LiquidTransportParams() {}
        ~LiquidTransportParams() {}

	//! Energies of molecular interaction associated with viscosity.
	/** 
	 * These multiply the mixture viscosity by
	 *  \f[ \exp( \sum_{i} \sum_{j} X_i X_j ( S_{i,j} + E_{i,j} / T ) ) \f].
	 *
	 * The overall formula for the logarithm of the mixture viscosity is 
	 *
	 * \f[ \ln \eta_{mix} = \sum_i X_i \ln \eta_i 
	 *  + \sum_i \sum_j X_i X_j ( S_{i,j} + E_{i,j} / T ) \f].
	 */
	DenseMatrix  visc_Eij; 

	//! Entropies of molecular interaction associated with viscosity.
	DenseMatrix  visc_Sij; 

        std::vector<Cantera::LiquidTransportData> LTData;

    };
}

#endif
