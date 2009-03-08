#ifndef CT_TRANSPORTPARAMS_H
#define CT_TRANSPORTPARAMS_H

#include <vector>
using namespace std;

#include "../ct_defs.h"
#include "TransportBase.h"
#include "../xml.h"
#include "../XML_Writer.h"

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
        vector<vector_fp>            visccoeffs;
        vector<vector_fp>            condcoeffs;
        vector<vector_fp>            diffcoeffs;
        vector_fp                    polytempvec;
        
        vector<vector<int> > poly;
        vector<vector_fp >   omega22_poly;
        vector<vector_fp >   astar_poly;
        vector<vector_fp >   bstar_poly;
        vector<vector_fp >   cstar_poly;

        vector_fp   zrot;
        vector_fp   crot;

        vector<bool> polar;
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
