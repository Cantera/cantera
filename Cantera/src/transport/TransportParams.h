#ifndef CT_TRANSPORTPARAMS_H
#define CT_TRANSPORTPARAMS_H

#include <vector>

#include "ct_defs.h"
#include "TransportBase.h"
#include "xml.h"
#include "XML_Writer.h"

namespace Cantera {

    /** 
     * Base class to hold transport model parameters.  
     * Used by TransportFactory.
     */
    class TransportParams {

    public:

        TransportParams() : thermo(0), xml(0) {}
        virtual ~TransportParams();
        size_t nsp_;

        //        phase_t* mix;
        thermo_t* thermo;
        vector_fp        mw;

        // polynomial fits
	//temperature-fit viscosity 
        std::vector<vector_fp>            visccoeffs; 
	//temperature-fit heat conduction 
        std::vector<vector_fp>            condcoeffs; 
	//temperature-fit diffusivity 
	std::vector<vector_fp>            diffcoeffs; 
        vector_fp                         polytempvec;

	//minimum and maximum temperatures for parameter fits
        doublereal tmax, tmin;
        int mode_;
        XML_Writer* xml;
        int log_level;

    };


    /**
     * Holds transport model parameters relevant to transport in ideal 
     * gases with a kinetic theory of gases derived transport model. 
     * Used by TransportFactory.
     */
    class GasTransportParams : public TransportParams {

    public:

        GasTransportParams() {}
        ~GasTransportParams() {}


        std::vector<std::vector<int> > poly;
        std::vector<vector_fp >   omega22_poly;
        std::vector<vector_fp >   astar_poly;
        std::vector<vector_fp >   bstar_poly;
        std::vector<vector_fp >   cstar_poly;

        vector_fp   zrot;
        vector_fp   crot;

        std::vector<bool> polar;
        vector_fp    alpha;
        vector_fp    fitlist;
        vector_fp    eps;
        vector_fp    sigma;
        DenseMatrix  reducedMass;    
        DenseMatrix  diam;           
        DenseMatrix  epsilon;        
        DenseMatrix  dipole;         
        DenseMatrix  delta;          

    };

}

#endif //CT_TRANSPORTPARAMS_H
