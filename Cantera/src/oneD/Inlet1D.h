/**
 * @file Inlet1D.h
 *
 * Boundary objects for one-dimensional simulations.
 *
 */

/*
 * $Author$
 * $Revision$
 * $Date$
 *
 * Copyright 2002-3  California Institute of Technology
 */


#ifndef CT_BDRY1D_H
#define CT_BDRY1D_H

#include "Domain1D.h"
#include "../SurfPhase.h"
#include "../InterfaceKinetics.h"
#include "StFlow.h"
#include "OneDim.h"
#include "../ctml.h"

namespace Cantera {

    const int LeftInlet = 1;
    const int RightInlet = -1;

    /**
     * The base class for boundaries between one-dimensional spatial
     * domains. The boundary may have its own internal variables, such
     * as surface species coverages.
     *
     * The boundary types are an inlet, an outlet, a symmetry plane,
     * and a surface.  
     *
     * The public methods are all virtual, and the base class
     * implementations throw exceptions.
     */
    class Bdry1D : public Domain1D {
    public:

        Bdry1D();
 
        virtual ~Bdry1D() {}

        /// Initialize.
        virtual void init() { _init(1); }

        /// Set the temperature.
        virtual void setTemperature(doublereal t){m_temp = t;}

        /// Temperature [K].
        virtual doublereal temperature() {return m_temp;}

        /// Set the mole fractions by specifying a string.
        virtual void setMoleFractions(string xin){err("setMoleFractions");}

        /// Set the mole fractions by specifying an array.
        virtual void setMoleFractions(doublereal* xin){err("setMoleFractions");}

        /// Mass fraction of species k.
        virtual doublereal massFraction(int k) {err("massFraction"); return 0.0;}

        /// Set the total mass flow rate.
        virtual void setMdot(doublereal mdot){m_mdot = mdot;}

        /// The total mass flow rate [kg/m2/s].
        virtual doublereal mdot() {return m_mdot;}

        virtual void _getInitialSoln(doublereal* x) {
            cout << "Bdry1D::_getInitialSoln called!  " << m_index << endl;
        }

    protected:

        void _init(int n);

        StFlow *m_flow_left, *m_flow_right;
        int m_ilr, m_left_nv, m_right_nv;
        int m_left_loc, m_right_loc;
        int m_left_points;
        int m_nv, m_left_nsp, m_right_nsp;
        int m_sp_left, m_sp_right;
        int m_start_left, m_start_right;
        ThermoPhase *m_phase_left, *m_phase_right;
        doublereal m_temp, m_mdot;

    private:
        void err(string method) {
            throw CanteraError("Bdry1D::"+method, 
                "attempt to call base class method "+method);
        }
    };


    /**
     * An inlet.
     */
    class Inlet1D : public Bdry1D {

    public:

        /**
         * Constructor. Create a new Inlet1D instance. If invoked
         * without parameters, a left inlet (facing right) is
         * constructed).
         */
        Inlet1D() : m_V0(0.0), m_nsp(0), m_flow(0) {
            m_type = cInletType;
            m_xstr = "";
        }
        virtual ~Inlet1D(){}

        /// set spreading rate
        virtual void setSpreadRate(doublereal V0) {
            m_V0 = V0;
            needJacUpdate();
        }

        /// spreading rate
        virtual double spreadRate() {
            return m_V0;
        }

        virtual void showSolution(ostream& s, const doublereal* x) {
            s << "-------------------  Inlet " << domainIndex() << " ------------------- " << endl;
            s << "  mdot:        " << m_mdot << " kg/m^2/s" << "   " << x[0] << endl;
            s << "  temperature: " << m_temp << " K" << "    " << x[1] << endl;
            if (m_flow) {
                s << "  mass fractions: " << endl;
                for (int k = 0; k < m_flow->phase().nSpecies(); k++) {
                    if (m_yin[k] != 0.0) {
                        s << "      " << m_flow->phase().speciesName(k) 
                          << "  " << m_yin[k] << endl;
                    }
                }
            }
            s << endl;
        }

        virtual void showSolution(const doublereal* x) {
            char buf[80];
            sprintf(buf, "    Mass Flux:   %10.4g kg/m^2/s \n", x[0]);
            writelog(buf);
            sprintf(buf, "    Temperature: %10.4g K \n", x[1]);
            writelog(buf);
            if (m_flow) {
                writelog("    Mass Fractions: \n");
                for (int k = 0; k < m_flow->phase().nSpecies(); k++) {
                    if (m_yin[k] != 0.0) {
                        sprintf(buf, "        %16s  %10.4g \n",
                            m_flow->phase().speciesName(k).c_str(), m_yin[k]);
                        writelog(buf);
                    }
                }
            }
            writelog("\n");
        }

        virtual void _getInitialSoln(doublereal* x) {
            x[0] = m_mdot;
            x[1] = m_temp;
        }

        virtual void _finalize(const doublereal* x) {
            ;//m_mdot = x[0];
            //m_temp = x[1];
        }

        virtual void setMoleFractions(string xin);
        virtual void setMoleFractions(doublereal* xin);
        virtual doublereal massFraction(int k) {return m_yin[k];}
        virtual string componentName(int n) const;
        virtual void init();
        virtual void eval(int jg, doublereal* xg, doublereal* rg, 
            integer* diagg, doublereal rdt);
        virtual void save(XML_Node& o, doublereal* soln);
        virtual void restore(XML_Node& dom, doublereal* soln);    

    protected:

        int m_ilr;
        doublereal m_V0;
        int m_nsp;
        vector_fp m_yin;
        string m_xstr;
        StFlow *m_flow;
    };


    /**
     * A symmetry plane. The axial velocity u = 0, and all other
     * components have zero axial gradients.
     */
    class Symm1D : public Bdry1D {

    public:

        Symm1D() {
            m_type = cSymmType; 
        }
        virtual ~Symm1D(){}

        virtual string componentName(int n) const;

        virtual void init();

        virtual void eval(int jg, doublereal* xg, doublereal* rg, 
            integer* diagg, doublereal rdt);

        virtual void save(XML_Node& o, doublereal* soln);
        virtual void restore(XML_Node& dom, doublereal* soln);    
        virtual void _finalize(const doublereal* x) {
            ; //m_temp = x[0];
        }
    protected:

    };


    /**
     */
    class Outlet1D : public Bdry1D {

    public:

        Outlet1D() {
            m_type = cOutletType; 
        }
        virtual ~Outlet1D(){}

        virtual string componentName(int n) const;

        virtual void init();

        virtual void eval(int jg, doublereal* xg, doublereal* rg, 
            integer* diagg, doublereal rdt);

        virtual void save(XML_Node& o, doublereal* soln);
        virtual void restore(XML_Node& dom, doublereal* soln);    
        virtual void _finalize(const doublereal* x) {
            ; //m_temp = x[0];
        }
    protected:

    };


    /**
     * A non-reacting surface. The axial velocity is zero
     * (impermeable), as is the transverse velocity (no slip). The
     * temperature is specified, and a zero flux condition is imposed
     * for the species.
     */
    class Surf1D : public Bdry1D {

    public:

        Surf1D() {
            m_type = cSurfType; 
        }
        virtual ~Surf1D(){}

        virtual string componentName(int n) const;

        virtual void init();

        virtual void eval(int jg, doublereal* xg, doublereal* rg, 
            integer* diagg, doublereal rdt);

        virtual void save(XML_Node& o, doublereal* soln);
        virtual void restore(XML_Node& dom, doublereal* soln);    

        virtual void _getInitialSoln(doublereal* x) {
            x[0] = m_temp;
        }

        virtual void _finalize(const doublereal* x) {
            ; //m_temp = x[0];
        }

        virtual void showSolution(ostream& s, const doublereal* x) {
            s << "-------------------  Surface " << domainIndex() << " ------------------- " << endl;
            s << "  temperature: " << m_temp << " K" << "    " << x[0] << endl;
        }

        virtual void showSolution(const doublereal* x) {
            char buf[80];
            sprintf(buf, "    Temperature: %10.4g K \n", x[0]);
            writelog(buf);
            writelog("\n");
        }

    protected:

    };

}

// //     };



// //     /////////////////////////////////////////////////////////////
// //     //
// //     // surf1D
// //     //
// //     ////////////////////////////////////////////////////////////
// // #ifdef WITH_CHEMSURF

// //     class ChemSurf1D : public Bdry1D {

// //     public:

// //         ChemSurf1D(InterfaceKinetics* skin = 0) {
// //             m_type = cSurfType;
// //             m_kin = 0;
// //             m_sphase = 0;
// //             m_nsp = 0;
// //             if (skin) setKinetics(skin);
// //         }
// //         virtual ~ChemSurf1D(){}

// //         void setKinetics(InterfaceKinetics* kin);

// //         virtual string componentName(int n) const;
// //         virtual void init();

// //         virtual void eval(int jg, doublereal* xg, doublereal* rg, 
// //             integer* diagg, doublereal rdt);

// //         virtual void save(XML_Node& o, doublereal* soln);
                
// //     protected:

// //         int m_ilr, m_nsp;
//         InterfaceKinetics* m_kin;
//         SurfPhase* m_sphase;
//         vector_fp m_work;
//         const doublereal *m_molwt_right, *m_molwt_left; 
//         int m_start_surf;
//         vector<ThermoPhase*> m_bulk;
//         vector<int> m_nbulk;
//         int m_nsurf;
//         vector_fp m_mult;
//         vector<bool> m_do_surf_species;
//         vector_fp m_fixed_cov;
//     };

#endif
