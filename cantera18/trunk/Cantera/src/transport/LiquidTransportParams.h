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


	//section for liquid transport properties

	//Arrhenius parameters for transport coefficients:

	//!Arrhenius pre-exponential parameter for viscosity.
	vector_fp  visc_A; 
	//!Temperature exponent for viscosity.
	vector_fp  visc_n; 
	//!Arrhenius activation temperature for viscosity.
	vector_fp  visc_Tact; 

	//!Arrhenius pre-exponential parameter for thermal conductivity.
	vector_fp  thermCond_A; 
	//!Temperature exponent for thermal conductivity.
	vector_fp  thermCond_n; 
	//!Arrhenius activation temperature for thermal conductivity.
	vector_fp  thermCond_Tact; 

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

	//Hydrodynamic radius of transported molecule
	vector_fp               hydroRadius;

        //! Coefficients for the limiting conductivity of ions 
        //! in solution: A_k
        /*!
         *  This is used in the following formula for the
         *  limiting conductivity of the kth ion.
         *
         *   ln (lambda^o_k nu_solv) = A_k + B_k / T
         *
         * nu_solv is the pure component solvent viscosity
         *
         *  Note the limiting conductivities of ions will also
         *  be used to input the diffusion coefficients. 
         */
        vector_fp A_k_cond;

        //! Coefficients for the limiting conductivity of ions
        //! in solution: B_k
        vector_fp B_k_cond;


        std::vector<Cantera::LiquidTransportData> LTData;

    };
}

#endif
