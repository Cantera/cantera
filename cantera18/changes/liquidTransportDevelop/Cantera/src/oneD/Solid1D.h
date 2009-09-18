/**
 * @file Solid1D.h
 *
 */

/*
 * $Author: dggoodwin $
 * $Revision: 1.3 $
 * $Date: 2006/11/27 21:43:34 $
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_SOLID1D_H
#define CT_SOLID1D_H

#include "../transport/TransportBase.h"
#include "Domain1D.h"
#include "../Array.h"
#include "../sort.h"
#include "../ThermoPhase.h"
#include "../Kinetics.h"
#include "../funcs.h"


namespace Cantera {

    class MultiJac;


    //-----------------------------------------------------------
    //  Class Solid1D
    //-----------------------------------------------------------


    /**
     *  A class for one-dimensional reacting solids with current
     *  transport.  This class implements the one-dimensional
     *  similarity solution for a chemically-reacting, axisymmetric,
     *  stagnation-point flow.
     */
    class Solid1D : public Domain1D {

    public:


        //------------------------------------------
        //   constants
        //------------------------------------------

        /**
         * Offsets of solution components in the solution array.
         */
        const unsigned int c_phi_loc;   // electric potential
        const unsigned int c_T_loc;     // temperature
        const unsigned int c_C_loc;     // concentrations


        //--------------------------------
        // construction and destruction
        //--------------------------------

        // Constructor.
        Solid1D(ThermoPhase* ph = 0, int nsp = 1, int points = 1);

        /// Destructor.
        virtual ~Solid1D(){}


        /**
         * @name Problem Specification
         */
        //@{

        virtual void setupGrid(int n, const doublereal* z);

        thermo_t& phase() { return *m_thermo; }
        kinetics_t& kinetics() { return *m_kin; }

        /**
         * Set the thermo manager.
         */
        void setThermo(thermo_t& th) { 
            m_thermo = &th;
        }

        /// set the kinetics manager
        void setKinetics(kinetics_t& kin) { m_kin = &kin; }

        /// set the transport manager
        void setTransport(Transport& trans);

        virtual void setState(int point, const doublereal* state) {
            setTemperature(point, state[c_T_loc]);
            setElectricPotential(point, state[c_phi_loc]);
            int k;
            for (k = 0; k < m_nsp; k++) {
                setConcentration(point, k, state[c_C_loc+k]);
            }
        }


        virtual void _getInitialSoln(doublereal* x) {
            int k, j;
            for (j = 0; j < m_points; j++) {
                x[index(c_T_loc,j)] = T_fixed(j);
                x[index(c_phi_loc,j)] = phi_fixed(j);
                for (k = 0; k < m_nsp; k++) {
                    x[index(c_C_loc+k,j)] = C_fixed(k,j);
                }
            }
        }   

        virtual void _finalize(const doublereal* x) {
            int k, j;
            doublereal zz, tt;
            int nz = m_zfix.size();
            bool e = m_do_energy[0];
            for (j = 0; j < m_points; j++) {
                if (e || nz == 0) 
                    setTemperature(j, T(x, j));
                else {
                    zz = (z(j) - z(0))/(z(m_points - 1) - z(0));
                    tt = linearInterp(zz, m_zfix, m_tfix);
                    setTemperature(j, tt);
                }   
                setElectricPotential(j, phi(x,j));
                for (k = 0; k < m_nsp; k++) {
                    setConcentration(j, k, C(x, k, j));
                }
            }
            if (e) solveEnergyEqn();
        }


        void setFixedTempProfile(vector_fp& zfixed, vector_fp& tfixed) {
            m_zfix = zfixed;
            m_tfix = tfixed;
        }

        /**
         * Set the temperature fixed point at grid point j, and
         * disable the energy equation so that the solution will be
         * held to this value.
         */
        void setTemperature(int j, doublereal t) {
            m_fixedtemp[j] = t;
            m_do_energy[j] = false;
        }

        /**
         * Set the electric potential fixed point at grid point j, and
         * disable Gauss's equation so that the solution will be
         * held to this value.
         */
        void setElectricPotential(int j, doublereal phi) {
            m_fixedphi[j] = phi;
            m_do_gauss[j] = false;
        }

        /**
         * Set the mass fraction fixed point for species k at grid
         * point j, and disable the species equation so that the
         * solution will be held to this value.
         */
        void setConcentration(int j, int k, doublereal c) {
            m_fixedc(k,j) = c;
            m_do_species[k] = true; // false;
        }

        /**
         * The fixed temperature value at point j.
         */
        doublereal T_fixed(int j) const {return m_fixedtemp[j];}

        /**
         * The fixed potential value at point j.
         */
        doublereal phi_fixed(int j) const {return m_fixedphi[j];}

        /**
         * The fixed mass fraction value of species k at point j.
         */
        doublereal C_fixed(int k, int j) const {return m_fixedc(k,j);}

        virtual std::string componentName(int n) const;

        void setDielectricConstant(doublereal e) { m_eps = e; }
        doublereal dielectricConstant() { return m_eps; }

        /**
         * Write a Tecplot zone corresponding to the current solution.
         * May be called multiple times to generate animation.
         */
        void outputTEC(ostream &s, const doublereal* x, 
            std::string title, int zone);

        virtual void showSolution(const doublereal* x);

        virtual void save(XML_Node& o, doublereal* sol);

        virtual void restore(XML_Node& dom, doublereal* soln);

        // overloaded in subclasses
        virtual std::string solidType() { return "<none>"; }

        void solveEnergyEqn(int j=-1) {
            if (j < 0)
                for (int i = 0; i < m_points; i++)
                    m_do_energy[i] = true;
            else 
                m_do_energy[j] = true;
            m_refiner->setActive(c_T_loc, true);
            needJacUpdate();
        }

        void fixTemperature(int j=-1) {
            if (j < 0)
                for (int i = 0; i < m_points; i++) {
                    m_do_energy[i] = false;
                }
            else m_do_energy[j] = false;
            m_refiner->setActive(c_T_loc, false);
            needJacUpdate();
        }

        void solveGaussEqn(int j=-1) {
            if (j < 0)
                for (int i = 0; i < m_points; i++)
                    m_do_gauss[i] = true;
            else 
                m_do_gauss[j] = true;
            m_refiner->setActive(c_phi_loc, true);
            needJacUpdate();
        }

        void fixElectricPotential(int j=-1) {
            if (j < 0)
                for (int i = 0; i < m_points; i++) {
                    m_do_gauss[i] = false;
                }
            else m_do_gauss[j] = false;
            m_refiner->setActive(c_phi_loc, false);
            needJacUpdate();
        }

        bool doSpecies(int k) { return m_do_species[k]; }
        bool doEnergy(int j) { return m_do_energy[j]; }
        bool doGauss(int j) { return m_do_gauss[j]; }

        void solveSpecies(int k=-1) {
            if (k == -1) {
                for (int i = 0; i < m_nsp; i++) 
                    m_do_species[i] = true;
            }
            else m_do_species[k] = true;
            needJacUpdate();
        }

        void fixSpecies(int k=-1) {
            if (k == -1) {
                for (int i = 0; i < m_nsp; i++) 
                    m_do_species[i] = false;
            }
            else m_do_species[k] = false;
            needJacUpdate();
        }

        void resize(int points);

        void setJac(MultiJac* jac);
        void setThermoState(const doublereal* x,int j);
        void setStateAtMidpoint(const doublereal* x,int j);


    protected:

        doublereal component(const doublereal* x, int i, int j) const {
            doublereal xx = x[index(i,j)];
            return xx;
        }

        doublereal wdot(int k, int j) const {return m_wdot(k,j);}

        /// write the net production rates at point j into array m_wdot
        void getWdot(doublereal* x,int j) { 
            setThermoState(x,j);
            m_kin->getNetProductionRates(&m_wdot(0,j));
        }

        /**
         * update the thermodynamic properties from point
         * j0 to point j1 (inclusive), based on solution x.
         */
        void updateThermo(const doublereal* x, int j0, int j1) {
            int j;
            for (j = j0; j <= j1; j++) {
                setThermoState(x,j);
                m_cp[j]  = m_thermo->cp_mass();
            }
        }


        //--------------------------------
        //      solution components
        //--------------------------------


        doublereal T(const doublereal* x,int j) const {
            return x[index(c_T_loc, j)];
        }

        doublereal& T(doublereal* x,int j) {return x[index(c_T_loc, j)];}

        doublereal T_prev(int j) const {return prevSoln(c_T_loc, j);} 

        doublereal C(const doublereal* x,int k, int j) const {
            return x[index(c_C_loc + k, j)];
        }

        doublereal& C(doublereal* x,int k, int j) {
            return x[index(c_C_loc + k, j)];
        }

        doublereal C_prev(int k, int j) const {
            return prevSoln(c_C_loc + k, j);
        }

        doublereal flux(int k, int j) const {
            return m_flux(k, j);
        }

        doublereal phi(const doublereal* x, int j) {
            return x[index(c_phi_loc, j)];
        }

        doublereal divHeatFlux(const doublereal* x, int j) const {
            doublereal c1 = m_tcon[j-1]*(T(x,j) - T(x,j-1));
            doublereal c2 = m_tcon[j]*(T(x,j+1) - T(x,j));
            return -2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
        }

        doublereal divDisplCurr(const doublereal* x, int j) const {
            doublereal c1 = (phi(x,j) - phi(x,j-1));
            doublereal c2 = (phi(x,j+1) - phi(x,j));
            return -2.0*m_eps*epsilon_0*
                (c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
        }

        void updateDiffFluxes(const doublereal* x, int j0, int j1);

        //---------------------------------------------------------
        //
        //             member data
        //
        //---------------------------------------------------------

        doublereal m_eps;           // relative dielectric constant

        // grid parameters
        vector_fp m_dz;

        // mixture thermo properties
        vector_fp m_cdens;

        // transport properties
        vector_fp m_tcon;
        vector_fp m_diff;
        Array2D m_flux;

        // production rates
        Array2D m_wdot;

        int m_nsp;

        thermo_t*       m_thermo;
        kinetics_t*     m_kin;
        Transport*      m_trans;

        MultiJac*       m_jac;

        bool m_ok;

        // flags
        std::vector<bool> m_do_energy;
        std::vector<bool> m_do_species;
        std::vector<bool> m_do_gauss;

        // fixed T and Y values
        Array2D   m_fixedy;
        Array2D   m_fixedphi;
        vector_fp m_fixedtemp;
        vector_fp m_zfix;
        vector_fp m_tfix;

    private:

        vector_fp m_cbar;
    };
}

#endif
