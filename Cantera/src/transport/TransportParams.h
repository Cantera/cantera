#ifndef CT_TRANSPORTPARAMS_H
#define CT_TRANSPORTPARAMS_H

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
    class TransportParams {

    public:

        TransportParams() : thermo(0), xml(0) {}
        virtual ~TransportParams();
        int nsp;

        //        phase_t* mix;
        thermo_t* thermo;
        vector_fp        mw;

        // polynomial fits
        std::vector<vector_fp>            visccoeffs;
        std::vector<vector_fp>            condcoeffs;
        std::vector<vector_fp>            diffcoeffs;
        vector_fp                    polytempvec;
        
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
        doublereal tmax, tmin;
        int mode;
        XML_Writer* xml;
        int log_level;
    };
}

#endif
