

/**
 *
 *  @file Thermo.h
 */

/*
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2001 California Institute of Technology
 *
 */

#ifndef CT_THERMO_H
#define CT_THERMO_H

#include "ThermoPhase.h"

// #include "ct_defs.h"
// #include "mix_defs.h"
// #include "Phase.h"

// #include "ctml.h"
// using namespace ctml;


// namespace Cantera {

//     /**
//      * Exception thrown if a method of class Thermo is called.  The
//      * methods of Thermo should be overloaded to implement a
//      * particular thermo manager. But a given manager may not overload
//      * every method. If an unimplemented virtual method is called, the
//      * base class method will throw an exception.
//      */
//     class ThermoNotImplemented : public CanteraError {
//     public:
//         ThermoNotImplemented(string method) : CanteraError("Thermo",
//             "**** Method "+method+" not implemented. ****\n") {}
//     };


//     /**
//      * Base class for thermodynamic property managers.
//      */
//     class Thermo  {

//     public:

//         Thermo(phase_t* phase=0, SpeciesThermo* sptherm = 0) {
//             if (phase == 0) m_s = new Phase();
//             m_s = phase; 
//             m_spthermo = sptherm;
//             m_xml = new XML_Node("thermo");
//             m_index = -1;
//         }

//         virtual ~Thermo() {}

//         int index() { return m_index; }
//         void setIndex(int m) { m_index = m; }

//         XML_Node& xml() { return *m_xml; }

//         /**
//          * Initialize.  @param s Object defining the composition and
//          * species properties of the phase that thermodynamic
//          * properties will be computed for.
//          */
//         virtual void initThermo(Phase& s) { 
//             m_s = &s;
//         }


//         /** Return a reference to the phase object. */
//         phase_t& phase() { return *m_s; }


//         /** Return a read-only reference to the phase object. */
//         const phase_t& phase() const { return *m_s; }


//         /** 
//          * Equation of state type. The base class returns
//          * zero. Subclasses should define this to return a non-zero
//          * value.
//          */
//         virtual int eosType() const { return 0; }

//         /**
//          * @name Virtual Methods
//          * 
//          * The methods in this section should be overloaded by subclasses. 
//          * The base class methods throw an exception. 
//          */

//         //@{ 

//         /**
//          * Molar enthalpy. Units: J/kmol. 
//          */
//         virtual doublereal enthalpy_mole() const {
//             return err("enthalpy_mole");
//         }

//         /**
//          * Molar internal energy. Units: J/kmol. 
//          */
//         virtual doublereal intEnergy_mole() const {
//             return err("intEnergy_mole");
//         }

//         /**
//          * Molar entropy. Units: J/kmol/K. 
//          */
//         virtual doublereal entropy_mole() const {
//             return err("entropy_mole");
//         }

//         /**
//          * Molar Gibbs function. Units: J/kmol. 
//          */
//         virtual doublereal gibbs_mole() const {
//             return err("gibbs_mole");
//         }

//         /**
//          * Molar heat capacity at constant pressure. Units: J/kmol/K. 
//          */
//         virtual doublereal cp_mole() const {
//             return err("cp_mole");
//         }

//         /**
//          * Molar heat capacity at constant volume. Units: J/kmol/K. 
//          */
//         virtual doublereal cv_mole() const {
//             return err("cv_mole");
//         }

//         /**
//          * Pressure. Units: Pa. 
//          */
//         virtual doublereal pressure() const {
//             return err("pressure");
//         }

//         /**
//          * Set the pressure. Units: Pa. 
//          */
//         virtual void setPressure(doublereal p) {
//             err("setPressure");
//         }

        
//         /**
//          * Get the species chemical potentials. Units: J/kmol.
//          */
//         virtual void getChemPotentials(doublereal* mu) const {
//             err("getChemPotentials_RT");
//         }

//         virtual void getChemPotentials_RT(doublereal* mu) const {
//             err("getChemPotentials_RT");
//         }

//         /**
//          * Get the species partial molar enthalpies. Units: J/kmol.
//          */
//         virtual void getPartialMolarEnthalpies(doublereal* hbar) const {
//             err("getPartialMolarEnthalpies");
//         }

//         /**
//          * Get the species partial molar entropies. Units: J/kmol.
//          */
//         virtual void getPartialMolarEntropies(doublereal* sbar) const {
//             err("getPartialMolarEntropies");
//         }

//         /**
//          * Get the species partial molar enthalpies. Units: J/kmol.
//          */
//         virtual void getPartialMolarVolumes(doublereal* vbar) const {
//             err("getPartialMolarVolumes");
//         }

//         /**
//          * Get the nondimensional Gibbs functions for the pure species
//          * at the current T and P.
//          */
//         virtual void getEnthalpy_RT(doublereal* hrt) const {
//             err("getEnthalpy_RT");
//         }


//         /**
//          * Get the nondimensional Gibbs functions for the pure species
//          * at the current T and P.
//          */
//         virtual void getEntropy_R(doublereal* sr) const {
//             err("getEntropy_R");
//         }

//         /**
//          * Get the nondimensional Gibbs functions for the pure species
//          * at the current T and P.
//          */
//         virtual void getGibbs_RT(doublereal* grt) const {
//             err("getGibbs_RT");
//         }

//         virtual void getPureGibbs(doublereal* gpure) const {
//             err("getPureGibbs");
//         }

//         /**
//          * Get the nondimensional Gibbs functions for the pure species
//          * at the current T and P.
//          */
//         virtual void getCp_R(doublereal* cpr) const {
//             err("getCp_RT");
//         }

//         //@}

//         virtual doublereal refPressure() const {
//             err("refPressure");
//             return 0.0;
//         }

//         virtual doublereal minTemp(int k = -1) {
//             err("minTemp");
//             return 0.0;
//         }

//         virtual doublereal maxTemp(int k = -1) {
//             err("maxTemp");
//             return 0.0;
//         }

//         /**
//          * Specific enthalpy. Units: J/kg. 
//          */
//         doublereal enthalpy_mass() const {
//             return enthalpy_mole()/m_s->meanMolecularWeight();
//         }

//         /**
//          * Specific internal energy. Units: J/kg. 
//          */
//         doublereal intEnergy_mass() const {
//             return intEnergy_mole()/m_s->meanMolecularWeight();
//         }

//         /**
//          * Specific entropy. Units: J/kg/K. 
//          */
//         doublereal entropy_mass() const {
//             return entropy_mole()/m_s->meanMolecularWeight();
//         }

//         /**
//          * Specific Gibbs function. Units: J/kg. 
//          */
//         doublereal gibbs_mass() const {
//             return gibbs_mole()/m_s->meanMolecularWeight();
//         }

//         /**
//          * Specific heat at constant pressure. Units: J/kg/K. 
//          */
//         doublereal cp_mass() const {
//             return cp_mole()/m_s->meanMolecularWeight();
//         }

//         /**
//          * Specific heat at constant volume. Units: J/kg/K. 
//          */
//         doublereal cv_mass() const {
//             return cv_mole()/m_s->meanMolecularWeight();
//         }

//         doublereal _temp() const {
//             return m_s->temperature();
//         }

//         doublereal _dens() const {
//             return m_s->density();
//         }

//         doublereal _RT() const {
//             return m_s->temperature() * GasConstant;
//         }

//         /** Set the temperature (K), pressure (Pa), and mole fractions.  */
//         void setState_TPX(doublereal t, doublereal p, const doublereal* x) {
//             m_s->setMoleFractions(x); m_s->setTemperature(t); setPressure(p);
//         }

//         /** Set the temperature (K), pressure (Pa), and mole fractions.  */
//         void setState_TPX(doublereal t, doublereal p, compositionMap& x) {
//             m_s->setMoleFractionsByName(x); m_s->setTemperature(t); setPressure(p);
//         }

//         /** Set the temperature (K), pressure (Pa), and mole fractions.  */
//         void setState_TPX(doublereal t, doublereal p, const string& x) {
//             compositionMap xx;
//             parseCompString(x, xx);
//             m_s->setMoleFractionsByName(xx); m_s->setTemperature(t); setPressure(p);
//         }        

//         /** Set the temperature (K), pressure (Pa), and mass fractions. */
//         void setState_TPY(doublereal t, doublereal p, const doublereal* y) {
//             m_s->setMassFractions(y); m_s->setTemperature(t); setPressure(p);
//         }

//         /** Set the temperature (K), pressure (Pa), and mass fractions. */
//         void setState_TPY(doublereal t, doublereal p, compositionMap& y) {
//             m_s->setMassFractionsByName(y); m_s->setTemperature(t); setPressure(p);
//         }
        
//         /** Set the temperature (K), pressure (Pa), and mass fractions.  */
//         void setState_TPY(doublereal t, doublereal p, const string& y) {
//             compositionMap yy;
//             parseCompString(y, yy);
//             m_s->setMassFractionsByName(yy); m_s->setTemperature(t); setPressure(p);
//         }

//         /** Set the temperature (K) and pressure (Pa) */
//         void setState_TP(doublereal t, doublereal p) {
//             m_s->setTemperature(t); setPressure(p);
//         }

//         /** Set the pressure (Pa) and mole fractions.  */
//         void setState_PX(doublereal p, doublereal* x) {
//             m_s->setMoleFractions(x); setPressure(p);
//         }

//         /** Set the pressure (Pa) and mass fractions.  */
//         void setState_PY(doublereal p, doublereal* y) {
//             m_s->setMassFractions(y); setPressure(p);
//         }

//         void setState_HP(doublereal h, doublereal p, doublereal tol = 1.e-8) {
//             doublereal dt;
//             setPressure(p);
//             for (int n = 0; n < 20; n++) {
//                 dt = (h - enthalpy_mass())/cp_mass();
//                 if (dt > 100.0) dt = 100.0;
//                 else if (dt < -100.0) dt = -100.0; 
//                 setState_TP(_temp() + dt, p);
//                 if (fabs(dt) < tol) {
//                     return;
//                 }
//             }
//             throw CanteraError("setState_HP","no convergence. dt = " + fp2str(dt));
//         }

//         void setState_UV(doublereal u, doublereal v, doublereal tol = 1.e-8) {
//             doublereal dt;
//             m_s->setDensity(1.0/v);
//             for (int n = 0; n < 20; n++) {
//                 dt = (u - intEnergy_mass())/cv_mass();
//                 if (dt > 100.0) dt = 100.0;
//                 else if (dt < -100.0) dt = -100.0; 
//                 m_s->setTemperature(_temp() + dt);
//                 if (fabs(dt) < tol) {
//                     return;
//                 }
//             }
//             throw CanteraError("setState_UV","no convergence. dt = " + fp2str(dt));
//         }

//         void setState_SP(doublereal s, doublereal p, doublereal tol = 1.e-8) {
//             doublereal dt;
//             setPressure(p);
//             for (int n = 0; n < 20; n++) {
//                 dt = (s - entropy_mass())*_temp()/cp_mass();
//                 if (dt > 100.0) dt = 100.0;
//                 else if (dt < -100.0) dt = -100.0; 
//                 setState_TP(_temp() + dt, p);
//                 if (fabs(dt) < tol) {
//                     return;
//                 }
//             }
//             throw CanteraError("setState_SP","no convergence. dt = " + fp2str(dt));
//         }

//         void setState_SV(doublereal s, doublereal v, doublereal tol = 1.e-8) {
//             doublereal dt;
//             m_s->setDensity(1.0/v);
//             for (int n = 0; n < 20; n++) {
//                 dt = (s - entropy_mass())*_temp()/cv_mass();
//                 if (dt > 100.0) dt = 100.0;
//                 else if (dt < -100.0) dt = -100.0; 
//                 m_s->setTemperature(_temp() + dt);
//                 if (fabs(dt) < tol) {
//                     return;
//                 }
//             }
//             throw CanteraError("setState_SV","no convergence. dt = " + fp2str(dt));
//         }

//         virtual void setToEquilState(const doublereal* lambda_RT) {
//             err("setToEquilState");
//         }

//         /// Install a standard-state species thermodynamic property
//         /// manager
//         void setSpeciesThermo(SpeciesThermo* spthermo) 
//             { m_spthermo = spthermo; }

//         SpeciesThermo& speciesThermo() { return *m_spthermo; }

//         virtual void setParameters(int n, doublereal* c) {}

//     protected:

//         Phase* m_s;
//         XML_Node* m_xml;
//         SpeciesThermo* m_spthermo;
//         int m_index;

//     private:

//         doublereal err(string msg) const {
//             throw ThermoNotImplemented(msg);
//             return 0;
//         }
//     };

//     typedef Thermo thermo_t;
// }
        
#endif





