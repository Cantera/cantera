/**
 *  @file FlowDevice.h
 *
 *  $Author$
 *  $Date$
 *  $Revision$
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_WALL_H
#define CT_WALL_H

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "../ct_defs.h"
#include "../Func1.h"

namespace Cantera {

    class ReactorBase;  // forward reference
    class Kinetics;
    class Func1;
    class SurfPhase;

    const int Rigid_Type = 1;
    const int Flexible_Type = 2;

    class Wall {

    public:

        /// Constructor
        Wall();

        /// Destructor 
        virtual ~Wall() {}

        /**
         * Rate of volume change (kg/s). Positive value increases
         * volume of reactor on left, and decreases volume on right.
         */
        virtual doublereal vdot(doublereal t);
        virtual doublereal Q(doublereal t);

        /// Area in m^2.
        doublereal area() { return m_area; }

        /// Set the area [m^2].
        void setArea(doublereal a) { m_area = a; }

        void setThermalResistance(doublereal Rth) { m_rrth = 1.0/Rth; }

        /// Set the overall heat transfer coefficient [W/m^2/K].
        void setHeatTransferCoeff(doublereal U) { m_rrth = U; }

        void setEmissivity(doublereal epsilon) { m_emiss = epsilon; }

        //  /** Set the rate of volume change to a specified function.*/
        //        void setExpansionRate(Func1* f=0) {if (f) m_vf = f;}

        /** Set the piston velocity to a specified function. */
        void setVelocity(Func1* f=0) {if (f) m_vf = f;}

        /** 
         * Set the expansion rate coefficient.
         */
        void setExpansionRateCoeff(doublereal k) {m_k = k;}

        /**
         * Specify the heat flux q(t).
         */
        void setHeatFlux(Func1* q) { m_qf = q;}

        bool install(ReactorBase& in, ReactorBase& out);
        virtual bool ready() { return (m_left != 0 && m_right != 0); }

        int type() { return 0; }

        /// Return a reference to the left reactor.
        ReactorBase& left() const { return *m_left; }

        /// Return a reference to the right-hand reactor.
        const ReactorBase& right() { return *m_right; }

        /// set parameters
        virtual void setParameters(int n, doublereal* coeffs) {
            m_coeffs.resize(n);
            copy(coeffs, coeffs + n, m_coeffs.begin());
        }

        void setKinetics(Kinetics* left = 0,
            Kinetics* right = 0);

        SurfPhase* surface(int leftright) {
            return m_surf[leftright];
        }

        Kinetics* kinetics(int leftright) {
            return m_chem[leftright];
        }

        void setCoverages(int leftright, const doublereal* cov);

        void getCoverages(int leftright, doublereal* cov);

        void syncCoverages(int leftright);


    protected:

        vector_fp m_coeffs;

        ReactorBase* m_left;
        ReactorBase* m_right;
        Kinetics * m_chem[2];
        SurfPhase* m_surf[2];
        int m_nsp[2];
        doublereal m_area, m_k, m_rrth;
        doublereal m_emiss;
        Func1 *m_vf;
        Func1 *m_qf;
        vector_fp m_leftcov, m_rightcov;

    private:

    };


//     class Piston : public Wall {
//     public: 
//         Piston() 
//             : m_omega(2.0*3.1415926*freq), Wall() {
//             //m_vdisp = stroke * Pi * bore * bore / 4.0;
//             //m_vclear = tdc * bore;
//             //m_ra = 1.0/crankradius;
//         }
//         ~Piston() {}
//         virtual doublereal vdot(double t) {
//             doublereal theta = m_omega * t;
//             doublereal sinth = sin(theta);
//             return 0.0;
//             //return m_vclear + 0.5*m_vdist*(1.0 + m_ra - m_omega*sin(theta)
//             //    - sqrt(m_ra * m_ra - sinth*sinth));
//         }
//     protected:
//         doublereal m_tdc, m_bdc, m_stroke, m_bore, m_rodlen,
//             m_radius, m_omega;
//     };


    class Piston : public Wall {
    public: 
        Piston(doublereal freq, 
            doublereal tdc, doublereal bdc,
            doublereal stroke, doublereal bore,
            doublereal rodlen, doublereal crankradius) 
            : Wall(), m_omega(2.0*3.1415926*freq) {
            //m_vdisp = stroke * Pi * bore * bore / 4.0;
            //m_vclear = tdc * bore;
            //m_ra = 1.0/crankradius;
        }
        virtual ~Piston() {}
        virtual doublereal vdot(double t) {
            //            doublereal theta = m_omega * t;
            //doublereal sinth = sin(theta);
            return 0.0;
            //return m_vclear + 0.5*m_vdist*(1.0 + m_ra - m_omega*sin(theta)
            //    - sqrt(m_ra * m_ra - sinth*sinth));
        }
    protected:
        doublereal m_tdc, m_bdc, m_stroke, m_bore, m_rodlen,
            m_radius, m_omega;
    };
        
}

#endif
