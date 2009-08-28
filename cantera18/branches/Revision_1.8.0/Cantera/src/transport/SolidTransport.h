/**
 *
 *  @file SolidTransport.h
 *   Header file defining class SolidTransport
 */

/* $Author: hkmoffa $
 * $Revision: 1.7 $
 * $Date: 2009/03/27 18:24:39 $
 */

// Copyright 2003  California Institute of Technology


#ifndef CT_SOLIDTRAN_H
#define CT_SOLIDTRAN_H


// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

// STL includes
#include <vector>
#include <string>
#include <map>
#include <numeric>
#include <algorithm>


// Cantera includes
#include "TransportBase.h"
#include "DenseMatrix.h"

namespace Cantera {


    /**
     * Class SolidTransport implements transport
     * properties for solids.
     */
    class SolidTransport : public Transport {

    public:
virtual ~SolidTransport() {}

        virtual int model() { return cSolidTransport; }

        virtual doublereal thermalConductivity();
        virtual void getMixDiffCoeffs(doublereal* const d);
        virtual void getMobilities(doublereal* const mobil);
        virtual void setParameters(const int n, const int k, const doublereal* const p);

        friend class TransportFactory;

    protected:

        /// default constructor
        SolidTransport();

    private:

        int m_nmobile;    // number of mobile species
        vector_fp m_Adiff;
        vector_fp m_Ndiff;
        vector_fp m_Ediff;
        vector_int m_sp;
        doublereal m_Alam;
        doublereal m_Nlam;
        doublereal m_Elam;
    };
}
#endif






