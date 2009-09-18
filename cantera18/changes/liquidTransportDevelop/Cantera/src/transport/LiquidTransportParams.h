#ifndef CT_LIQUIDTRANSPORTPARAMS_H
#define CT_LIQUIDTRANSPORTPARAMS_H

#include <vector>

#include "ct_defs.h"
#include "TransportBase.h"
#include "xml.h"
#include "XML_Writer.h"

namespace Cantera {

    /**
     *
     * Holds transport data. Used by TransportFactory.
     *
     */
    class LiquidTransportParams {

    public:

        LiquidTransportParams() : thermo(0), xml(0) {}
        virtual ~LiquidTransportParams();
        int nsp;

        //        phase_t* mix;
        thermo_t* thermo;
        vector_fp        mw;

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


        // polynomial fits
        std::vector<vector_fp>  viscCoeffsVector_;
        std::vector<vector_fp>  condcoeffs;
        std::vector<vector_fp>  diffcoeffs ;


        std::vector<bool> polar;
        //vector_fp    alpha;
        vector_fp    fitlist;
        vector_fp    eps;
        vector_fp    sigma;
        DenseMatrix  reducedMass;  
        DenseMatrix  diam;           
        DenseMatrix  epsilon;        
        DenseMatrix  dipole;         
        DenseMatrix  delta;          
        doublereal tmax, tmin;
        int mode;
        XML_Writer* xml;
        int log_level;
    };
}

#endif
