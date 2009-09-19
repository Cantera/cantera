#ifndef CT_LIQUIDTRANSPORTPARAMS_H
#define CT_LIQUIDTRANSPORTPARAMS_H

#include <vector>

#include "ct_defs.h"
#include "TransportBase.h"
#include "TransportParams.h"
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

	//Arrhenius parameters for transport coefficients
	//	std::vector<vector_fp>  viscParams; 
	vector_fp  visc_A; 
	vector_fp  visc_n; 
	vector_fp  visc_Tact; 
	vector_fp  thermCond_A; 
	vector_fp  thermCond_n; 
	vector_fp  thermCond_Tact; 
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


    };
}

#endif
