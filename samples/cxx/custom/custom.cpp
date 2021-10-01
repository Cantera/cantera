/*!
 * @file custom.cpp
 *
 * Custom Reactor
 *
 * Solve a closed-system constant pressure ignition problem where the governing
 * equations are custom-implemented, using Cantera's interface to CVODES to
 * integrate the equations.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/base/Solution.h"
#include "cantera/numerics/Integrator.h"
#include <fstream>

using namespace Cantera;

class ReactorODEs : public FuncEval {
public:
    /**
     * Constructor
     * @param[in] sol Solution object specifying initial system state.
     */
    ReactorODEs(shared_ptr<Solution> sol) {
        /* ---------------------- INITIALIZE MEMBER VARS ---------------------- */

        // pointer to the system's ThermoPhase object. updated by the solver during
        // simulation to provide iteration-specific thermodynamic properties.
        m_gas = sol->thermo();

        // pointer to the kinetics manager. provides iteration-specific species
        // production rates based on the current state of the ThermoPhase.
        m_kinetics = sol->kinetics();

        // the system's constant pressure, taken from the provided initial state.
        m_pressure = m_gas->pressure();

        // number of chemical species in the system.
        m_nSpecies = m_gas->nSpecies();

        // resize the vector_fp storage containers for species partial molar enthalpies
        // and net production rates. internal values are updated and used by the solver
        // per iteration.
        m_hbar.resize(m_nSpecies);
        m_wdot.resize(m_nSpecies);

        // number of equations in the ODE system. a conservation equation for each
        // species, plus a single energy conservation equation for the system.
        m_nEqs = m_nSpecies + 1;
    }

    /**
     * Evaluate the ODE right-hand-side function, ydot = f(t,y).
     *   - overridden from FuncEval, called by the integrator during simulation.
     * @param[in] t time.
     * @param[in] y solution vector, length neq()
     * @param[out] ydot rate of change of solution vector, length neq()
     * @param[in] p sensitivity parameter vector, length nparams()
     *   - note: sensitivity analysis isn't implemented in this example
     */
    void eval(double t, double* y, double* ydot, double* p) {
        // the solution vector *y* is [T, Y_1, Y_2, ... Y_K], where T is the
        // system temperature, and Y_k is the mass fraction of species k.
        // similarly, the time derivative of the solution vector, *ydot*, is
        // [dT/dt, Y_1/dt, Y_2/dt, ... Y_K/dt].
        // the following variables are defined for clear and convenient access
        // to these vectors:
        double temperature = y[0];
        double *massFracs = &y[1];
        double *dTdt = &ydot[0];
        double *dYdt = &ydot[1];
        
        /* ------------------------- UPDATE GAS STATE ------------------------- */
        // the state of the ThermoPhase is updated to reflect the current solution
        // vector, which was calculated by the integrator.
        m_gas->setMassFractions_NoNorm(massFracs);
        m_gas->setState_TP(temperature, m_pressure);

        /* ----------------------- GET REQ'D PROPERTIES ----------------------- */
        double rho = m_gas->density();
        double cp = m_gas->cp_mass();
        m_gas->getPartialMolarEnthalpies(&m_hbar[0]);
        m_kinetics->getNetProductionRates(&m_wdot[0]);
        
        /* -------------------------- ENERGY EQUATION ------------------------- */
        // the rate of change of the system temperature is found using the energy
        // equation for a closed-system constant pressure ideal gas:
        //     m*cp*dT/dt = - sum[h(k) * dm(k)/dt]
        // or equivalently:
        //     dT/dt = - sum[hbar(k) * dw(k)/dt] / (rho * cp)
        double hdot_vol = 0;
        for (size_t k = 0; k < m_nSpecies; k++) {
            hdot_vol += m_hbar[k] * m_wdot[k];
        }
        *dTdt = - hdot_vol / (rho * cp);

        /* --------------------- SPECIES CONSERVATION EQS --------------------- */
        // the rate of change of each species' mass fraction is found using the closed-system
        // species conservation equation, applied once for each species:
        //     m*dY(k)/dt = dm(k)/dt
        // or equivalently:
        //     dY(k)/dt = dw(k)/dt * MW(k) / rho
        for (size_t k = 0; k < m_nSpecies; k++) {
            dYdt[k] = m_wdot[k] * m_gas->molecularWeight(k) / rho;
        }
    }
    
    /**
     * Number of equations in the ODE system.
     *   - overridden from FuncEval, called by the integrator during initialization.
     */
    size_t neq() {
        return m_nEqs;
    }

    /**
     * Provide the current values of the state vector, *y*.
     *   - overridden from FuncEval, called by the integrator during initialization.
     * @param[out] y solution vector, length neq()
     */
    void getState(double* y) {
        // the solution vector *y* is [T, Y_1, Y_2, ... Y_K], where T is the
        // system temperature, and Y_k is the mass fraction of species k.
        y[0] = m_gas->temperature();
        m_gas->getMassFractions(&y[1]);
    }

private:
    // private member variables, to be used internally.
    shared_ptr<ThermoPhase> m_gas;
    shared_ptr<Kinetics> m_kinetics;
    vector_fp m_hbar;
    vector_fp m_wdot;
    double m_pressure;
    size_t m_nSpecies;
    size_t m_nEqs;
};

int main() {
    /* -------------------- CREATE GAS & SPECIFY STATE -------------------- */
    auto sol = newSolution("h2o2.yaml", "ohmech", "None");
    auto gas = sol->thermo();
    gas->setState_TPX(1001, OneAtm, "H2:2, O2:1, N2:4");

    /* ---------------------- CREATE CSV OUTPUT FILE ---------------------- */
    // simulation results will be outputted to a .csv file as complete state vectors
    // for each simulation time point.
    // create the csv file, overwriting any existing ones with the same name.
    std::ofstream outputFile("custom_cxx.csv");
    
    // for convenience, a title for each of the state vector's components is written to
    // the first line of the csv file.
    outputFile << "time (s), temp (K)";
    for (size_t k = 0; k < gas->nSpecies(); k++) {
        outputFile << ", Y_" << gas->speciesName(k);
    }
    outputFile << std::endl;

    /* --------------------- CREATE ODE RHS EVALUATOR --------------------- */
    ReactorODEs odes = ReactorODEs(sol);

    /* ---------------------- SPECIFY TIME INTERVAL ----------------------- */
    // the simulation is run over the time interval specified below. tnow is initialized
    // with the simulation's start time, and is updated on each timestep to reflect
    // the new absolute time of the system.
    double tnow = 0.0;
    double tfinal = 1e-3;
    
    /* ------------------- CREATE & INIT ODE INTEGRATOR ------------------- */
    // create an ODE integrator object, which will be used to solve the system of ODES defined
    // in the ReactorODEs class. a C++ interface to the C-implemented SUNDIALS CVODES integrator
    // (CVodesIntegrator) is built into Cantera, and will be used to solve this example.
    //  - the default settings for CVodesIntegrator are used:
    //     solution method: BDF_Method
    //     problem type: DENSE + NOJAC
    //     relative tolerance: 1.0e-9
    //     absolute tolerance: 1.0e-15
    //     max step size: +inf
    std::unique_ptr<Integrator> integrator(newIntegrator("CVODE"));

    // initialize the integrator, specifying the start time and the RHS evaluator object.
    // internally, the integrator will apply settings, allocate needed memory, and populate
    // this memory with the appropriate initial values for the system.
    integrator->initialize(tnow, odes);
    
    /* ----------------------- SIMULATION TIME LOOP ----------------------- */
    while (tnow < tfinal) {
        // advance the simulation to the current absolute time, tnow, using the integrator's
        // ODE system time-integration capability. a pointer to the resulting system state vector
        // is fetched in order to access the solution components.
        integrator->integrate(tnow);
        double *solution = integrator->solution();

        // output the current absolute time and solution state vector to the csv file. you can view
        // these results by opening the "custom_cxx.csv" file that appears in this program's directory
        // after compiling and running.
        outputFile << tnow;
        for (size_t i = 0; i < odes.neq(); i++) {
            outputFile << ", " << solution[i];
        }
        outputFile << std::endl;

        // increment the simulation's absolute time, tnow, then return to the start of the loop to
        // advance the simulation to this new time point.
        tnow += 1e-5;
    }
}
