/**
 *
 */

#ifndef CT_FLOW_BOUNDS
#define CT_FLOW_BOUNDS

#include <vector>
#include "ct_defs.h"
//#include "surfKinetics.h"

using namespace Cantera;

namespace FlowBdry {

    class Boundary {

    public:
        Boundary(int nsp = 0) 
            : rtau(1.e5), m_lr(1), m_nsp(nsp), m_nv(nsp+4),
              m_bv(0), m_bp(0), m_temp(0.0), m_mdot(0.0), m_V(0.0) { 
            m_y.resize(nsp, 0.1);
        }
        virtual ~Boundary(){}

        doublereal rtau;
        doublereal T() { return m_temp; }
        doublereal mdot() { return m_mdot; }
        doublereal V() { return m_V; }
        doublereal* Y() { return m_y.begin(); }

        int orient() { return m_lr; }
        int faceLeft() { m_lr = -1; return m_lr; }
        int faceRight() { m_lr = 1; return m_lr; }

        void setNSpecies(int n) { m_nsp = n; }
        void set_mdot(doublereal mdot) { m_mdot = mdot; }
        void set_V(doublereal V) { m_V = V; }
        void set_T(doublereal T) { m_temp = T; }
        void set_Y(const doublereal* y) {
            copy(y, y+m_nsp, m_y.begin());
        }

        virtual void evalInt(int j, doublereal* x, doublereal* r) {
            int i;
            for (i = 0; i < m_bv; i++) {
                r[i] = x[i];
            }
        }
        virtual void eval(doublereal* x0, doublereal density, 
            doublereal* data, doublereal* r) {
            throw CanteraError("Boundary::eval",
                " ERROR!! Base class eval called!");
        }

        virtual int nIntVars() { return m_bv; }
        virtual int nIntPoints() { return m_bp; }

    protected:

        int offset(int n) { return m_lr*n*m_nv; }

        int m_lr, m_nsp, m_nv;
        int m_bv, m_bp;
        doublereal m_temp, m_mdot, m_V;
        vector_fp m_y;
    };


    /**
     * zero axial velocity, zero gradients for everything else.
     */
    class SymmPlane : public Boundary {

    public:
        SymmPlane(int nsp) : Boundary(nsp) {}
        virtual ~SymmPlane() {}

        virtual void eval(doublereal* x0, doublereal density,
            doublereal* data, doublereal* r) {
            for (int n = 0; n < m_nv; n++) {
                r[n] = x0[n] - x0[offset(1) + n];
            }
            r[0] = x0[0];
        }

    protected:
        int m_iloc;
    };



    /**
     * Zero gradients for all components, except V, which is zero.
     */
    class Outlet : public Boundary {

    public:

        Outlet(int nsp) : Boundary(nsp) {}
        virtual ~Outlet() {}

        virtual void eval(doublereal* x0, doublereal density,
            doublereal* data, doublereal* r) {
            for (int n = 0; n < m_nv; n++) {
                r[n] = -rtau*(x0[n] - x0[offset(1) + n]);
            }
            r[1] = -rtau*x0[1];
        }
    };


#ifdef INCL_SURF

    class Surface : public Boundary {

    public:

        Surface(int nsp, SurfKinetics* kin=0) 
            : Boundary(nsp), m_kin(kin), m_dt(1.e3) {
            if (m_kin) 
                m_sdot.resize(m_kin->nTotal());
            else
                m_sdot.resize(nsp);
        }
        virtual ~Surface(){}


        virtual void eval(doublereal* x0, doublereal density, 
            doublereal* data, doublereal* r) {

            doublereal sum = 0.0, mdot = 0.0;
            //            updateSurfaceRates(x0[2], density, x0+4);

            for (int k = 0; k < m_nsp; k++) {
                //                r[4+k] = m_sdot[k] - m_lr*data[k];
                //mdot += m_sdot[k];
                //sum += x0[4+k];
                r[4+k] = - m_lr*data[k];
                sum += x0[4+k];
            }
            r[4] = 1.0 - sum;

            // match the Stefan velocity at the surface                
            r[0] = density*x0[0] - m_lr*mdot;
            r[1] = x0[1];              // no slip 
            r[2] = x0[2] - m_temp;     // specified T
            r[3] = x0[3];
        }

        /**
         * Update the mass species production rates, given
         * the species mass fractions at the surface.
         * @param y
         */
        virtual void updateSurfaceRates(doublereal t, 
            doublereal rho, doublereal* y) {
            int k;
            if (m_kin) {
                for (k = 0; k < m_nsp; k++) {
                    m_sdot[k] = fmaxx(Tiny, y[k]);
                }
                m_kin->bulkPhase(0)->setState_TR(t, rho);
                m_kin->bulkPhase(0)->setMassFractions_NoNorm(m_sdot.begin());
                //m_kin->integrate(m_dt);
                m_kin->getNetProductionRates(m_sdot.begin());
            }
            else
                for (k = 0; k < m_nsp; k++) m_sdot[k] = 0.0;

        }
            
            

    protected:

        const doublereal* m_wt;
        vector_fp   m_sdot;
        SurfKinetics* m_kin;
        doublereal m_dt;
    };

#endif

    class Inlet : public Boundary {

    public:
        Inlet(int nsp, doublereal* wt=0) 
            : Boundary(nsp) {}
        virtual ~Inlet(){}

        virtual void eval(doublereal* x0, 
            doublereal density, doublereal* flux, doublereal* r) {
            int k;
            for (k = 0; k < m_nsp; k++) {
                r[4+k] = rtau*(m_y[k] - x0[k+4] - flux[k]/m_mdot);
            }
            r[0] = -rtau*(density*x0[0] - m_lr*m_mdot);
            r[1] = -rtau*(x0[1] - V());
            r[2] = -rtau*(x0[2] - T());     // specified T
            r[3] = -rtau*x0[3];
        }
    };

}



#endif
