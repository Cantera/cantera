/**
 *  @file FtnTransport.h
 *
 *  Customizable Fortran transport manager. This manager calls
 *  external Fortran functions to evaluate the transport
 *  properties. This is designed to be used to build custom transport
 *  managers in Fortran.
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_FTNTRANSPORT_H
#define CT_FTNTRANSPORT_H

#include "ct_defs.h"
#include "DenseMatrix.h"


/**
 * Change these definitions to change the names of the Fortran
 * procedures. If you want to define more than one custom Fortran
 * transport manager, copy this file, rename the class, and change
 * these definitions.
 *
 * The Fortran procedure names must follow the conventions of the Fortran
 * compiler. On most unix systems, this means the names must be lowercase,
 * and must include a trailing underscore.
 *
 */
#define __VISC__           visc_
#define __BULKVISC__       bvisc_
#define __TCON__           tcon_
#define __SPVISC__         spvisc_
#define __SPCOND__         spcond_
#define __SPFLUXES__       spfluxes_
#define __TDIFF__          tdiff_
#define __MULTIDIFF__      multidiff_
#define __MIXDIFF__        mixdiff_
#define __UPT__            updatet_
#define __UPC__            updatec_
#define __INIT__           init_

extern "C" {
    doublereal __VISC__(doublereal* t, doublereal* p, doublereal* x);
    doublereal __BULKVISC__(doublereal* t, doublereal* p, doublereal* x);
    doublereal __TCON__(doublereal* t, doublereal* p, doublereal* x);

    void __SPVISC__(doublereal* t, doublereal* p, doublereal* x, doublereal* visc);
    void __SPCOND__(doublereal* t, doublereal* p, doublereal* x, doublereal* cond);
    void __SPFLUXES__(doublereal* t, doublereal* p, doublereal* x,
        integer* ndim, doublereal* gradt, integer* ldx, doublereal* gradx,
        integer* ldf, doublereal* fluxes);
    void __TDIFF__(doublereal* t, doublereal* p, doublereal* x, doublereal* dt);
    void __MULTIDIFF__(doublereal* t, doublereal* p, doublereal* x, 
        integer* ld, doublereal* d);
    void __MIXDIFF__(doublereal* t, doublereal* p, doublereal* x, doublereal* d);
    void __BINDIFF__(doublereal* t, doublereal* p, doublereal* x, doublereal* d);
    void __UPT__(doublereal* t);
    void __UPC__(doublereal* p, doublereal* x);
    void __INIT__();
}

namespace Cantera {

    /**
     * A class that calls external Fortran functions to evaluate 
     * transport properties. Not currently used - may need updating.
     */
    class FtnTransport : public Transport {

    public:

        FtnTransport(int model) { m_model = model; }

        virtual int model() { return m_model; }

        virtual doublereal viscosity() { 
            return __VISC__(&m_temp, &m_pres, m_x.begin()); 
        }

        virtual void getSpeciesViscosities(doublereal* visc) { 
            __SPVISC__(&m_temp, &m_pres, m_x.begin(), visc); 
        } 

        virtual void getSpeciesConductivities(doublereal* cond) { 
            __SPCOND__(&m_temp, &m_pres, m_x.begin(), cond); 
        } 

        virtual doublereal bulkViscosity()  
            { return __BULKVISC__(&m_temp, &m_pres, m_x.begin()); }

        virtual doublereal thermalConductivity()
            { return __TCON__(&m_temp, &m_pres, m_x.begin()); }

        virtual void getSpeciesFluxes(doublereal p, int ndim, 
        doublereal* grad_T, int ldx, doublereal* grad_X,
            int ldf, doublereal* fluxes) {
            doublereal pp = p;
            integer ldxx = ldx, ndimm = ndim, ldff = ldf;
            __SPFLUXES__(&m_temp, &pp, &m_x, &ndimm, grad_T, &ldxx, grad_X,
                &ldff, fluxes);
        }

        virtual void getThermalDiffCoeffs(doublereal* dt) 
            { __TDIFF__(&m_temp, &m_pres, m_x.begin(), dt); }

        virtual void getBinaryDiffCoeffs(doublereal p, int ld, doublereal* d) 
            { m_pres = p;
              integer ldd = ld;
              __BINDIFF__(&m_temp, &m_pres, m_x.begin(), &ldd, d);
            }

        virtual void getMultiDiffCoeffs(doublereal p, int ld, doublereal* d)
            { m_pres = p;
              integer ldd = ld;
              __MULTIDIFF__(&m_temp, &m_pres, m_x.begin(), &ldd, d);
            }

        virtual void getMixDiffCoeffs(doublereal p, doublereal* d) 
            { m_pres = p;
              integer ldd = ld;
              __MIXDIFF__(&m_temp, &m_pres, m_x.begin(), d);
            }            

        virtual void update_T() 
            { m_temp = m_mix->temperature();
              __UPT__(m_temp);
            }

        virtual void update_C()
            { m_mix->getMoleFractions(m_x.begin());
              m_pres = m_mix->pressure();
              __UPC__(m_pres, m_x);
            }

        virtual bool init(TransportParams& tr) { 
            m_mix = tr.mix;
            __INIT__();
        }

    private:

        doublereal m_temp;
        doublereal m_pres;
        vector_fp  m_x;
        int m_model;

    };

}
#endif






