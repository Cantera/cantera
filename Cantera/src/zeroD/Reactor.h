/**
 *  @file Reactor.h
 *
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_REACTOR_H
#define CT_REACTOR_H

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ReactorBase.h"
#include "../FuncEval.h"
#include "../CVode.h"
#include "../Kinetics.h"

namespace Cantera {

    /**
     * Class Reactor is a general-purpose class for stirred
     * reactors. The reactor may have an arbitrary number of inlets
     * and outlets, each of which may be connected to a "flow device"
     * such as a mass flow controller, a pressure regulator,
     * etc. Additional reactors may be connected to the other end of
     * the flow device, allowing construction of arbitrary reactor
     * networks.
     *
     * The reactor class integrates the same governing equations no
     * mattter what type of reactor is simulated. The differences
     * among reactor types are completely specified by the attached
     * flow devices and the time-dependent user-specified boundary
     * conditions. 
     *
     * If an instance of class Reactor is used directly, it will
     * simulate an adiabatic, constant volume reactor with gas-phase
     * chemistry but no surface chemistry. Other reactor types may be
     * simulated by deriving a class from Reactor and overloading
     * method getParams.  This method allows specifying the following
     * in terms of the instantaneous reactor state:
     *
     *  - rate of change of the total volume (m^3/s) 
     *  - surface heat loss rate (W) 
     *  - species surface production rates (kmol/s)
     * 
     * class Reactor inherits from both ReactorBase and
     * FuncEval. ReactorBase provides the basic reactor-like methods
     * that FlowDevice instances can access to determine their mass
     * flow rate. Class FuncEval is the class used to define a system
     * of ODE's to be integrated.
     */

    class Reactor : public ReactorBase, public FuncEval {

    public:

        /**
         * Default constructor.
         */
        Reactor();


        /**
         * Destructor. Deletes the integrator.
         */
        virtual ~Reactor(){ delete m_integ; }
        
        virtual int type() const { return ReactorType; }

        /** 
         * Advance the state of the reactor in time. On the first
         * call, internal method 'initialize' is called, and the maximum
         * integrator step size is set. By default, this is set to
         * 'time'. To specify a different maximum step size, precede the
         * call to advance with a call to setMaxStep. Note that this
         * cannot be reset after advance has been called.
         * 
         * @param time Final time (s).
         */
        virtual void advance(doublereal time) {
            if (!m_init) {
                setMaxStep(time);
                initialize();
            }
            m_integ->integrate(time);
            m_time = time;
            updateState(m_integ->solution());
            m_mix->saveState(m_state);
        }

        virtual double step(doublereal time) {
            if (!m_init) {
                setMaxStep(time);
                initialize();
            }
            m_time = m_integ->step(time);
            updateState(m_integ->solution());
            m_mix->saveState(m_state);
            return m_time;
        }

        /**
         * Insert something into the reactor. The 'something' must
         * belong to a class that is a subclass of both ThermoPhase
         * and Kinetics.
         */
        template<class G>
        void insert(G& contents) {
            setThermoMgr(contents);
            setKineticsMgr(contents);
        }

        void setKineticsMgr(Kinetics& kin) {
            m_kin = &kin;
        }

        /**
         * Set the maximum step size for integration.
         */
        void setMaxStep(doublereal maxstep) {
            m_maxstep = maxstep;
        }

        /**
         * Set the reactor surface area [m$^2$]. Can be changed at any time.
         */
        void setArea(doublereal area) {
            m_area = area;
        }

        /**
         * Set the external temperature \f$ T_0 \f$ 
         * used for heat loss calculations.
         * The heat loss rate is calculated from
         * \f[
         * \dot Q_{out} = h A (T - T_0) + \epsilon A (T^4 - T_{0,R}^4).
         * \f]
         * @see setArea, setEmissivity, setExtRadTemp
         */
        void setExtTemp(doublereal ts) {
            m_ext_temp = ts;
            if (!m_trad_set) m_ext_temp4 = ts*ts*ts*ts;
        }

        /**
         * Set the external temperature for radiation. By default, this
         * is the same as the temperature set by setExtTemp. But if 
         * setExtRadTemp is called, then subsequent of calls to 
         * setExtTemp do not modify the value set here.
         */
        void setExtRadTemp(doublereal tr) {        
            m_ext_temp4 = tr*tr*tr*tr;
        }

        void setHeatTransferCoeff(doublereal h) {
            m_h = h;
        }

        void setVDotCoeff(doublereal k) {
            m_kv = k;
        }

        void setEmissivity(doublereal emis) {
            m_emis = emis;
        }

        void setExtPressure(doublereal p0) {
            m_p0 = p0;
        }

        void disableChemistry() { m_chem = false; }
        void enableChemistry() { m_chem = true; }

        /// Set the energy equation on or off.
        void setEnergy(int eflag = 1) { 
            if (eflag > 0) m_energy = true;
            else m_energy = false;
        } 

        //-----------------------------------------------------

        /** @name References to internal objects */
        //@{

        /// Return a reference to the integrator.
        Integrator& integrator() { return *m_integ; }

        //@}


        //-----------------------------------------------------

        // overloaded methods of class FuncEval
        virtual int neq() { return m_nsp + 2; }
	virtual void eval(doublereal t, doublereal* y, doublereal* ydot);
        virtual void getInitialConditions(doublereal t0, size_t leny, 
            doublereal* y);
        virtual void initialize(doublereal t0 = 0.0);



        //-----------------------------------------------------

        /**
         * @name Methods to specify simulation options.
         * These virtual methods may be overloaded in 
         * derived classes to implement models for heat gain/loss,
         * surface chemistry, and compression/expansion.
         */
        //@{

        /**
         * Initialize the boundary conditions, if necessary. This
         * method does nothing, but may be overloaded in derived classes if
         * initialization is needed.
         */
        //        virtual void initBC() {}

        /**
         * Evaluate the reactor boundary conditions. This procedure is
         * called during integration to evaluate the rate of volume
         * change \f$ dV/dt \f$ [m^3/s], the heat loss rate [W], and
         * the species production rates due to surface chemistry.
         *
         * It may be overloaded in derived classes to implement other
         * boundary conditions. If not overloaded, this routine 
         * implements the following boundary conditions.
         *
         * The rate of volume change is
         * \f[
         * dV/dt = K ( P - P_{ext})
         * \f]
         * where K is set in procedure setVDotCoeff.
         *
         * 
         * The heat loss rate is calculated from
         * \f[
         * \dot Q_{out} = h A (T - T_0) + \epsilon A (T^4 - T_{0,R}^4).
         * \f]
         * @see setArea, setEmissivity, setExtRadTemp
         */

//         virtual void evalBC(doublereal& vdot, 
//             doublereal& heatLossRate, doublereal* sdot) {
//             doublereal t = m_mix->temperature();

//             //m_p0 = m_env->pressure();
//             vdot = m_kv * (m_thermo->pressure()/m_p0 - 1.0)*m_vol0;
//             heatLossRate = m_area * (
//                 m_h * (t - m_ext_temp)
//                 + m_emis * StefanBoltz * (t*t*t*t - m_ext_temp4)
//                 );
//         }

        //@}

        //-----------------------------------------------------

        /**
         * Set the mixture to a state consistent with solution
         * vector y.
         */

    protected:

        virtual void updateState(doublereal* y);
        
        Kinetics*   m_kin;
        //        ReactorBase*     m_env;
        //        Thermo*     m_thermo;

        Integrator* m_integ;         // pointer to integrator
        doublereal m_temp_atol;      // tolerance on T
        doublereal m_maxstep;        // max step size
        doublereal m_vdot, m_Q;
        doublereal m_emis, m_h, m_area;
        doublereal m_ext_temp, m_ext_temp4;
        doublereal m_kv, m_p0;
        vector_fp m_atol;
        doublereal m_rtol;
        vector_fp m_sdot;            // surface production rates
        bool m_trad_set;
        bool m_chem;
        bool m_energy;

    private:
    };
}

#endif

