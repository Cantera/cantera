/**
 * @file StFlow.h
 *
 */

/*
 * $Author: dggoodwin $
 * $Revision: 1.18 $
 * $Date: 2007/05/22 20:16:47 $
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_STFLOW_H
#define CT_STFLOW_H

#include "TransportBase.h"
#include "Domain1D.h"
#include "Array.h"
#include "IdealGasPhase.h"
#include "Kinetics.h"
#include "funcs.h"
//#include "../flowBoundaries.h"


namespace Cantera {

    typedef IdealGasPhase igthermo_t;

    class MultiJac;


    //------------------------------------------
    //   constants
    //------------------------------------------

    // Offsets of solution components in the solution array.
    const unsigned int c_offset_U = 0;    // axial velocity
    const unsigned int c_offset_V = 1;    // strain rate
    const unsigned int c_offset_T = 2;    // temperature
    const unsigned int c_offset_L = 3;    // (1/r)dP/dr
    const unsigned int c_offset_Y = 4;    // mass fractions

    // Transport option flags
    const int c_Mixav_Transport = 0;
    const int c_Multi_Transport = 1;
    const int c_Soret = 2;



    //-----------------------------------------------------------
    //  Class StFlow
    //-----------------------------------------------------------


    /**
     *  This class represents 1D flow domains that satisfy the
     *  one-dimensional similarity solution for chemically-reacting,
     *  axisymmetric, flows.
     */
    class StFlow : public Domain1D {

    public:

        //--------------------------------
        // construction and destruction
        //--------------------------------

        /// Constructor. Create a new flow domain.
        /// @param gas Object representing the gas phase. This object
        /// will be used to evaluate all thermodynamic, kinetic, and transport
        /// properties.
        /// @param nsp Number of species.
        StFlow(igthermo_t* ph = 0, int nsp = 1, int points = 1);

        /// Destructor.
        virtual ~StFlow(){}

        /**
         * @name Problem Specification
         */
        //@{

        virtual void setupGrid(int n, const doublereal* z);

        thermo_t& phase() { return *m_thermo; }
        kinetics_t& kinetics() { return *m_kin; }

        virtual void init(){
		}

		/**
         * Set the thermo manager. Note that the flow equations assume
         * the ideal gas equation.
         */
        void setThermo(igthermo_t& th) { m_thermo = &th; }

        /// Set the kinetics manager. The kinetics manager must
        void setKinetics(kinetics_t& kin) { m_kin = &kin; }

        /// set the transport manager
        void setTransport(Transport& trans, bool withSoret = false);
        void enableSoret(bool withSoret);
        bool withSoret() const { return m_do_soret; }

        /// Set the pressure. Since the flow equations are for the limit of
        /// small Mach number, the pressure is very nearly constant
        /// throughout the flow.
        void setPressure(doublereal p) { m_press = p; }


        /// @todo remove? may be unused
        virtual void setState(int point, const doublereal* state,
                              doublereal *x) {
            setTemperature(point, state[2]);
            int k;
            for (k = 0; k < m_nsp; k++) {
                setMassFraction(point, k, state[4+k]);
            }
        }


        /// Write the initial solution estimate into
        /// array x.
        virtual void _getInitialSoln(doublereal* x) {
            int k, j;
            for (j = 0; j < m_points; j++) {
                x[index(2,j)] = T_fixed(j);
                for (k = 0; k < m_nsp; k++) {
                    x[index(4+k,j)] = Y_fixed(k,j);
                }
            }
        }

        virtual void _finalize(const doublereal* x);


        /// Sometimes it is desired to carry out the simulation
        /// using a specified temperature profile, rather than
        /// computing it by solving the energy equation. This
        /// method specifies this profile.
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
         * Set the mass fraction fixed point for species k at grid
         * point j, and disable the species equation so that the
         * solution will be held to this value.
         * note: in practice, the species are hardly ever held fixed.
         */
        void setMassFraction(int j, int k, doublereal y) {
            m_fixedy(k,j) = y;
            m_do_species[k] = true; // false;
        }


         /// The fixed temperature value at point j.
        doublereal T_fixed(int j) const {return m_fixedtemp[j];}


        /// The fixed mass fraction value of species k at point j.
        doublereal Y_fixed(int k, int j) const {return m_fixedy(k,j);}


        virtual std::string componentName(int n) const;

        //added by Karl Meredith
        int componentIndex(std::string name) const;


        virtual void showSolution(const doublereal* x);

        virtual void save(XML_Node& o, doublereal* sol);

        virtual void restore(const XML_Node& dom, doublereal* soln);

        // overloaded in subclasses
        virtual std::string flowType() { return "<none>"; }

        void solveEnergyEqn(int j=-1) {
            if (j < 0)
                for (int i = 0; i < m_points; i++)
                    m_do_energy[i] = true;
            else
                m_do_energy[j] = true;
            m_refiner->setActive(0, true);
            m_refiner->setActive(1, true);
            m_refiner->setActive(2, true);
            needJacUpdate();
        }

        void fixTemperature(int j=-1) {
            if (j < 0)
                for (int i = 0; i < m_points; i++) {
                    m_do_energy[i] = false;
                }
            else m_do_energy[j] = false;
            m_refiner->setActive(0, false);
            m_refiner->setActive(1, false);
            m_refiner->setActive(2, false);
            needJacUpdate();
        }

        bool doSpecies(int k) { return m_do_species[k]; }
        bool doEnergy(int j) { return m_do_energy[j]; }

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

        void integrateChem(doublereal* x,doublereal dt);

        void resize(int components, int points);

        virtual void setFixedPoint(int j0, doublereal t0){}


        void setJac(MultiJac* jac);
        void setGas(const doublereal* x,int j);
        void setGasAtMidpoint(const doublereal* x,int j);

        //Karl Meredith
        //        doublereal density_unprotected(int j) const {
        //    return m_rho[j];
        // }
        doublereal density(int j) const {
            return m_rho[j];
        }

        virtual bool fixed_mdot() { return true; }
        void setViscosityFlag(bool dovisc) { m_dovisc = dovisc; }

    protected:

        doublereal component(const doublereal* x, int i, int j) const {
            doublereal xx = x[index(i,j)];
            return xx;
        }

        doublereal conc(const doublereal* x,int k,int j) const {
            return Y(x,k,j)*density(j)/m_wt[k];
        }

        doublereal cbar(const doublereal* x,int k, int j) const {
            return std::sqrt(8.0*GasConstant * T(x,j) / (Pi * m_wt[k]));
        }

        doublereal wdot(int k, int j) const {return m_wdot(k,j);}

        /// write the net production rates at point j into array m_wdot
        void getWdot(doublereal* x,int j) {
            setGas(x,j);
            m_kin->getNetProductionRates(&m_wdot(0,j));
        }

        /**
         * update the thermodynamic properties from point
         * j0 to point j1 (inclusive), based on solution x.
         */
        void updateThermo(const doublereal* x, int j0, int j1) {
            int j;
            for (j = j0; j <= j1; j++) {
                setGas(x,j);
                m_rho[j] = m_thermo->density();
                m_wtm[j] = m_thermo->meanMolecularWeight();
                m_cp[j]  = m_thermo->cp_mass();
            }
        }


        //--------------------------------
        // central-differenced derivatives
        //--------------------------------

        doublereal cdif2(const doublereal* x, int n, int j,
            const doublereal* f) const {
            doublereal c1 = (f[j] + f[j-1])*(x[index(n,j)] - x[index(n,j-1)]);
            doublereal c2 = (f[j+1] + f[j])*(x[index(n,j+1)] - x[index(n,j)]);
            return (c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
        }


        //--------------------------------
        //      solution components
        //--------------------------------


        doublereal T(const doublereal* x,int j) const {
            return x[index(c_offset_T, j)];
        }
        doublereal& T(doublereal* x,int j) {return x[index(c_offset_T, j)];}
        doublereal T_prev(int j) const {return prevSoln(c_offset_T, j);}

        doublereal rho_u(const doublereal* x,int j) const {
            return m_rho[j]*x[index(c_offset_U, j)];}

        doublereal u(const doublereal* x,int j) const {
            return x[index(c_offset_U, j)];}

        doublereal V(const doublereal* x,int j) const {
            return x[index(c_offset_V, j)];}
        doublereal V_prev(int j) const {
            return prevSoln(c_offset_V, j);}

        doublereal lambda(const doublereal* x,int j) const {
            return x[index(c_offset_L, j)];
        }

        doublereal Y(const doublereal* x,int k, int j) const {
            return x[index(c_offset_Y + k, j)];
        }

        doublereal& Y(doublereal* x,int k, int j) {
            return x[index(c_offset_Y + k, j)];
        }

        doublereal Y_prev(int k, int j) const {
            return prevSoln(c_offset_Y + k, j);
        }

        doublereal X(const doublereal* x,int k, int j) const {
            return m_wtm[j]*Y(x,k,j)/m_wt[k];
        }

        doublereal flux(int k, int j) const {
            return m_flux(k, j);
        }


        // convective spatial derivatives. These use upwind
        // differencing, assuming u(z) is negative

        doublereal dVdz(const doublereal* x,int j) const {
            int jloc = (u(x,j) > 0.0 ? j : j + 1);
            return (V(x,jloc) - V(x,jloc-1))/m_dz[jloc-1];
        }

        doublereal dYdz(const doublereal* x,int k, int j) const {
            int jloc = (u(x,j) > 0.0 ? j : j + 1);
            return (Y(x,k,jloc) - Y(x,k,jloc-1))/m_dz[jloc-1];
        }

        doublereal dTdz(const doublereal* x,int j) const {
            int jloc = (u(x,j) > 0.0 ? j : j + 1);
            return (T(x,jloc) - T(x,jloc-1))/m_dz[jloc-1];
        }

        doublereal shear(const doublereal* x,int j) const {
            doublereal c1 = m_visc[j-1]*(V(x,j) - V(x,j-1));
            doublereal c2 = m_visc[j]*(V(x,j+1) - V(x,j));
            return 2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
        }

        doublereal divHeatFlux(const doublereal* x, int j) const {
            doublereal c1 = m_tcon[j-1]*(T(x,j) - T(x,j-1));
            doublereal c2 = m_tcon[j]*(T(x,j+1) - T(x,j));
            return -2.0*(c2/(z(j+1) - z(j)) - c1/(z(j) - z(j-1)))/(z(j+1) - z(j-1));
        }

        int mindex(int k, int j, int m) {
            return m*m_nsp*m_nsp + m_nsp*j + k;
        }

        void updateDiffFluxes(const doublereal* x, int j0, int j1);


        //---------------------------------------------------------
        //
        //             member data
        //
        //---------------------------------------------------------

        // inlet
        doublereal m_inlet_u;
        doublereal m_inlet_V;
        doublereal m_inlet_T;
        doublereal m_rho_inlet;
        vector_fp m_yin;

        // surface
        doublereal m_surface_T;

        doublereal m_press;        // pressure

        // grid parameters
        vector_fp m_dz;
        //vector_fp m_z;

        // mixture thermo properties
        vector_fp m_rho;
        vector_fp m_wtm;

        // species thermo properties
        vector_fp m_wt;
        vector_fp m_cp;
        vector_fp m_enth;

        // transport properties
        vector_fp m_visc;
        vector_fp m_tcon;
        vector_fp m_diff;
        vector_fp m_multidiff;
        Array2D m_dthermal;
        Array2D m_flux;

        // production rates
        Array2D m_wdot;
        vector_fp m_surfdot;

        int m_nsp;

        igthermo_t*     m_thermo;
        kinetics_t*     m_kin;
        Transport*      m_trans;

        MultiJac*       m_jac;

        bool m_ok;

        // flags
        std::vector<bool> m_do_energy;
        bool m_do_soret;
        std::vector<bool> m_do_species;
        int m_transport_option;

        // solution estimate
        //vector_fp m_zest;
        //Array2D   m_yest;

        // fixed T and Y values
        Array2D   m_fixedy;
        vector_fp m_fixedtemp;
        vector_fp m_zfix;
        vector_fp m_tfix;

        doublereal m_efctr;
        bool m_dovisc;
        void updateTransport(doublereal* x,int j0, int j1);

    private:
        vector_fp m_ybar;

    };


    /**
     * A class for axisymmetric stagnation flows.
     */
    class AxiStagnFlow : public StFlow {
    public:
        AxiStagnFlow(igthermo_t* ph = 0, int nsp = 1, int points = 1) :
            StFlow(ph, nsp, points) { m_dovisc = true; }
        virtual ~AxiStagnFlow() {}
        virtual void eval(int j, doublereal* x, doublereal* r,
            integer* mask, doublereal rdt);
        virtual std::string flowType() { return "Axisymmetric Stagnation"; }
    };

    /**
     * A class for freely-propagating premixed flames.
     */
    class FreeFlame : public StFlow {
    public:
        FreeFlame(igthermo_t* ph = 0, int nsp = 1, int points = 1) :
            StFlow(ph, nsp, points) { 
            m_dovisc = false;
            setID("flame");
        }
        virtual ~FreeFlame() {}
        virtual void eval(int j, doublereal* x, doublereal* r,
            integer* mask, doublereal rdt);
        virtual std::string flowType() { return "Free Flame"; }
        virtual bool fixed_mdot() { return false; }
    };


    /*
    class OneDFlow : public StFlow {
    public:
        OneDFlow(igthermo_t* ph = 0, int nsp = 1, int points = 1) :
            StFlow(ph, nsp, points) {
        }
        virtual ~OneDFlow() {}
        virtual void eval(int j, doublereal* x, doublereal* r,
            integer* mask, doublereal rdt);
        virtual std::string flowType() { return "OneDFlow"; }
        doublereal mdot(doublereal* x, int j) {
            return x[index(c_offset_L,j)];
        }

    private:
        void updateTransport(doublereal* x,int j0, int j1);
    };
    */

    void importSolution(doublereal* oldSoln, igthermo_t& oldmech,
        doublereal* newSoln, igthermo_t& newmech);

}

#endif
