/**
 *  @file Reactor.cpp
 *
 *  A zero-dimensional reactor
 */
 
// Copyright 2001  California Institute of Technology


#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "Reactor.h"
#include "../CVode.h"
#include "FlowDevice.h"
#include "Wall.h"
#include "../InterfaceKinetics.h"
#include "../SurfPhase.h"

namespace Cantera {

    doublereal quadInterp(doublereal x0, doublereal* x, doublereal* y);

    Reactor::Reactor() : ReactorBase(), 
                         FuncEval(),
                         m_kin(0),
                         m_integ(0),
                         m_temp_atol(1.e-11), 
                         m_maxstep(0.0),
                         m_vdot(0.0), 
                         m_Q(0.0), 
                         m_rtol(1.e-9),
                         m_chem(true),
                         m_energy(true)
    {
        m_integ = new CVodeInt;

        // use backward differencing, with a full Jacobian computed
        // numerically, and use a Newton linear iterator
        m_integ->setMethod(BDF_Method);
        m_integ->setProblemType(DENSE + NOJAC);
        m_integ->setIterator(Newton_Iter);        
    }


    // overloaded method of FuncEval. Called by the integrator to
    // get the initial conditions.
    void Reactor::getInitialConditions(double t0, size_t leny, double* y) 
    {
        m_init = true;
        if (m_mix == 0) {
            cout << "Error: reactor is empty." << endl;
            return;
        }  
        m_time = t0;

        // total mass
        doublereal mass = m_mix->density() * m_vol;
        
        // set components y + 2 ... y + K + 1 to the  
        // mass M_k of each species
        m_mix->getMassFractions(y+2);
        scale(y + 2, y + m_nsp + 2, y + 2, mass);
            
        // set the first component to the total internal 
        // energy
        y[0] = m_thermo->intEnergy_mass() * mass;
        
        // set the second component to the total volume
        y[1] = m_vol;

        // set the remaining components to the surface species
        // coverages on the walls
        int loc = m_nsp + 2;
        SurfPhase* surf;
        for (int m = 0; m < m_nwalls; m++) {
            surf = m_wall[m]->surface(m_lr[m]);
            if (surf) {
                surf->getCoverages(y+loc);
                loc += surf->nSpecies();
            }
        }
    }


    /**
     *  Must be called before calling method 'advance'
     */
    void Reactor::initialize(doublereal t0) {
        m_mix->restoreState(m_state);
        m_sdot.resize(m_nsp, 0.0);
        m_nv = m_nsp + 2;
        for (int w = 0; w < m_nwalls; w++)
            if (m_wall[w]->surface(m_lr[w]))
                m_nv += m_wall[w]->surface(m_lr[w])->nSpecies();
        m_atol.resize(neq());
        fill(m_atol.begin(), m_atol.end(), 1.e-15);
        m_integ->setTolerances(m_rtol, neq(), m_atol.begin());
        m_integ->setMaxStep(m_maxstep);
        m_integ->initialize(t0, *this);

        m_enthalpy = m_thermo->enthalpy_mass();
        m_pressure = m_thermo->pressure();
        m_intEnergy = m_thermo->intEnergy_mass();

        int nt = 0, maxnt = 0;
        for (int m = 0; m < m_nwalls; m++) {
            if (m_wall[m]->kinetics(m_lr[m])) {
                nt = m_wall[m]->kinetics(m_lr[m])->nTotalSpecies();
                if (nt > maxnt) maxnt = nt;
                if (m_wall[m]->kinetics(m_lr[m])) {
                    if (&m_kin->thermo(0) != 
                        &m_wall[m]->kinetics(m_lr[m])->thermo(0)) {
                        throw CanteraError("Reactor::initialize",
                            "First phase of all kinetics managers must be"
                            " the gas.");
                    }
                }
            }   
        }
        m_work.resize(maxnt);

        m_init = true;
    }
    
    void Reactor::updateState(doublereal* y) {

        phase_t& mix = *m_mix;  // define for readability

        // The components of y are the total internal energy,
        // the total volume, and the mass of each species.

        // Set the mass fractions and  density of the mixture.

        doublereal u   = y[0];
        m_vol          = y[1];
        doublereal* mss = y + 2;
        doublereal mass = accumulate(y+2, y+2+m_nsp, 0.0);
        m_mix->setMassFractions(mss);
        m_mix->setDensity(mass/m_vol);

        doublereal temp = temperature();
        mix.setTemperature(temp);

        if (m_energy) {
            doublereal u_mass = u/mass;       // specific int. energy
            doublereal delta;

            do {
                delta = -(m_thermo->intEnergy_mass() 
                    - u_mass)/m_thermo->cv_mass();
                temp += delta;
                mix.setTemperature(temp);
            }
            while (fabs(delta) > m_temp_atol);
        }
        mix.setTemperature(temp);
        m_state[0] = temp;

        int loc = m_nsp + 2;
        SurfPhase* surf;
        for (int m = 0; m < m_nwalls; m++) {
            surf = m_wall[m]->surface(m_lr[m]);
            if (surf) {
                surf->setTemperature(temp);
                surf->setCoverages(y+loc);
                loc += surf->nSpecies();
            }
        }

        // save parameters needed by other connected reactors
        m_enthalpy = m_thermo->enthalpy_mass();
        m_pressure = m_thermo->pressure();
        m_intEnergy = m_thermo->intEnergy_mass();

        m_mix->saveState(m_state);
    }


    void Reactor::eval(doublereal time, doublereal* y, doublereal* ydot) 
    {
        updateState(y);          // synchronize the reactor state with y
        evalEqs(time, y, ydot);
    }

    /**
     * Called by the integrator to evaluate ydot given y at time 'time'.
     */
    void Reactor::evalEqs(doublereal time, doublereal* y, doublereal* ydot) 
    {
        int i, k, nk;
        m_time = time;
        m_mix->restoreState(m_state);

        //        updateState(y);          // synchronize the reactor state with y

        m_vdot = 0.0;
        m_Q    = 0.0;

        // compute wall terms
        doublereal vdot, rs0, sum, wallarea;
        Kinetics* kin;
        SurfPhase* surf;
        int lr, ns, loc = m_nsp+2, surfloc;
        fill(m_sdot.begin(), m_sdot.end(), 0.0);
        for (i = 0; i < m_nwalls; i++) {
            lr = 1 - 2*m_lr[i];
            vdot = lr*m_wall[i]->vdot(time);
            m_vdot += vdot;
            m_Q += lr*m_wall[i]->Q(time);
            kin = m_wall[i]->kinetics(m_lr[i]);
            surf = m_wall[i]->surface(m_lr[i]);
            if (surf && kin) {
                rs0 = 1.0/surf->siteDensity();
                nk = surf->nSpecies();
                sum = 0.0;
                kin->getNetProductionRates(m_work.begin());
                ns = kin->surfacePhaseIndex();
                surfloc = kin->kineticsSpeciesIndex(0,ns);
                for (k = 1; k < nk; k++) {
                    ydot[loc + k] = m_work[surfloc+k]*rs0*surf->size(k);
                    sum -= ydot[loc + k];
                }
                ydot[loc] = sum;
                loc += nk;

                wallarea = m_wall[i]->area();
                for (k = 0; k < m_nsp; k++) {
                    m_sdot[k] += m_work[k]*wallarea;
                }
            }
        } 

        // volume equation
        ydot[1] = m_vdot;

        /* species equations
         *  Equation is:
         *  \dot M_k = \hat W_k \dot\omega_k + \dot m_{in} Y_{k,in}
         *             - \dot m_{out} Y_{k} + A \dot s_k.
         */
        const doublereal* mw = m_mix->molecularWeights().begin();

        int n;
        if (m_chem) {
            m_kin->getNetProductionRates(ydot+2);   // "omega dot"
        }
        else {
            fill(ydot + 2, ydot + 2 + m_nsp, 0.0);
        }
        for (n = 0; n < m_nsp; n++) {
            ydot[n+2] *= m_vol;     //           moles/s/m^3 -> moles/s
            ydot[n+2] += m_sdot[n]; 
            ydot[n+2] *= mw[n];
        }


        /**
         *  Energy equation.
         * \dot U = -P\dot V + A \dot q + \dot m_{in} h_{in}
         * - \dot m_{out} h.
         */
        if (m_energy) {
            ydot[0] = - m_thermo->pressure() * m_vdot - m_Q;
        }
        else {
            ydot[0] = 0.0;    
        }

        // add terms for open system
        if (m_open) {

            const doublereal* mf = m_mix->massFractions();
            doublereal enthalpy = m_thermo->enthalpy_mass();

            // outlets 

            int n;
            doublereal mdot_out;
            for (i = 0; i < m_nOutlets; i++) {
                mdot_out = m_outlet[i]->massFlowRate();
                for (n = 0; n < m_nsp; n++) {
                    ydot[2+n] -= mdot_out * mf[n];
                }
                if (m_energy)
                    ydot[0] -= mdot_out * enthalpy;
            }


            // inlets

            doublereal mdot_in;
            for (i = 0; i < m_nInlets; i++) {
                mdot_in = m_inlet[i]->massFlowRate();
                for (n = 0; n < m_nsp; n++) {
                    ydot[2+n] += m_inlet[i]->massFlowRate(n);
                }
                if (m_energy)
                    ydot[0] += mdot_in * m_inlet[i]->enthalpy_mass();
            }
        }
    }
}
