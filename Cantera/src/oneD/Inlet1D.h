
#ifndef CT_BDRY1D_H
#define CT_BDRY1D_H

#include "Resid1D.h"
//#include "surfacePhase.h"
//#include "surfKinetics.h"
#include "StFlow.h"
#include "OneDim.h"
#include "ctml.h"

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
    class Bdry1D : public Resid1D {
    public:
        Bdry1D() : Resid1D(1, 1, 0.0) {}
        virtual ~Bdry1D() {}

        /// Initialize.
        virtual void init(){err("init");}

        /// Set the temperature.
        virtual void setTemperature(doublereal t){err("setTemperature");}

        /// Temperature [K].
        virtual doublereal temperature() {err("temperature"); return 0.0;}

        /// Set the mole fractions by specifying a string.
        virtual void setMoleFractions(string xin){err("setMoleFractions");}

        /// Set the mole fractions by specifying an array.
        virtual void setMoleFractions(doublereal* xin){err("setMoleFractions");}
        /// Mass fraction of species k.
        virtual doublereal massFraction(int k) {err("massFraction"); return 0.0;}
        /// Set the total mass flow rate.
        virtual void setMdot(doublereal mdot){err("setMdot");}

        /// The total mass flow rate [kg/m2/s].
        virtual doublereal mdot() {err("mdot"); return 0.0;}

    protected:
    private:
        void err(string method) {
            throw CanteraError("Bdry1D::"+method, "attempt to call base class method "+method);
        }
    };


    class Inlet1D : public Bdry1D {

    public:

        Inlet1D(int ilr = 1) {
            m_type = cInletType; 
            m_flow = 0;
            m_ilr = ilr;
        }
        virtual ~Inlet1D(){}

        /// Set the inlet temperature
        virtual void setTemperature(doublereal t) {
            m_temp = t;
            needJacUpdate();
        }

        /// set spreading rate
        virtual void setSpreadRate(doublereal V0) {
            m_V0 = V0;
            needJacUpdate();
        }

        /// Temperature [K].
        doublereal temperature() {return m_temp;}

        virtual void setMoleFractions(string xin) {
            m_xstr = xin;
            if (m_flow) {
                m_flow->phase().setMoleFractionsByName(xin);
                m_flow->phase().getMassFractions(m_yin.begin());
                needJacUpdate();
            }
        }

        virtual void setMoleFractions(doublereal* xin) {
            if (m_flow) {
                m_flow->phase().setMoleFractions(xin);
                m_flow->phase().getMassFractions(m_yin.begin());
                needJacUpdate();
            }
        }

        virtual doublereal massFraction(int k) {return m_yin[k];}

        virtual void setMdot(doublereal mdot) { m_mdot = mdot; }

        virtual string componentName(int n) const { 
            switch (n) {
            case 0: return "mdot"; break;
            case 1: return "temperature"; break;
            default: return "unknown";
            }
        }

        virtual void init() {
            if (m_index < 0) {
                throw CanteraError("Inlet1D", 
                    "install in container before calling init.");
            }
            resize(2,1);

            // set bounds
            const doublereal lower[2] = {-1.0e5, 200.0};
            const doublereal upper[2] = {1.0e5, 1.e5};
            setBounds(2, lower, 2, upper);

            // set tolerances
            vector_fp rtol(2, 1e-4);
            vector_fp atol(2, 1.e-5);
            setTolerances(2, rtol.begin(), 2, atol.begin());

            // if a flow domain is present on the left, then this must
            // be a right inlet
            if (m_index > 0) {
                Resid1D& r = container().domain(m_index-1);
                if (r.domainType() == cFlowType) {
                    m_ilr = RightInlet;
                    m_flow = (StFlow*)&r;
                }
                else 
                    throw CanteraError("Inlet1D::init",
                        "Inlet domains can only be connected to a flow domain.");
            }
            else {
                if (container().nDomains() > 1) {
                    Resid1D& r = container().domain(1);
                    if (r.domainType() == cFlowType) {
                        m_ilr = LeftInlet;
                        m_flow = (StFlow*)&r;
                    }
                    else 
                        throw CanteraError("Inlet1D::init",
                            "An inlet domain can only be connected to a flow domain.");
                }
                else 
                    throw CanteraError("Inlet1D::init",
                        "An inlet domain must be connected to a flow domain.");
            }
            // components = u, V, T, lambda, + mass fractions
            m_nsp = m_flow->nComponents() - 4;
            m_yin.resize(m_nsp, 0.0);
            if (m_xstr != "") 
                setMoleFractions(m_xstr);
            else
                m_yin[0] = 1.0;
        }


        virtual void eval(int jg, doublereal* xg, doublereal* rg, 
            integer* diagg, doublereal rdt) {
            int k;
            if (jg >= 0 && (jg < firstPoint() - 2 || jg > lastPoint() + 2)) return;

            // start of local part of global arrays
            doublereal* x = xg + loc();
            doublereal* r = rg + loc();
            integer* diag = diagg + loc();
            doublereal *xb, *rb;

            // residual equations for the two local variables
            r[0] = m_mdot - x[0];
            r[1] = m_temp - x[1];

            // both are algebraic constraints
            diag[0] = 0;
            diag[1] = 0;

            // if it is a left inlet, then the flow solution vector
            // starts 2 to the right in the global solution vector
            if (m_ilr == LeftInlet) {
                xb = x + 2;
                rb = r + 2;

                // If the energy equation is being solved, then
                // the flow domain set this residual to T(0). 
                // Subtract the inlet temperature.
                if (m_flow->doEnergy(0)) {
                    rb[2] -= x[1];   // T
                }

                // spreading rate. Flow domain sets this to V(0),
                // so for finite spreading rate subtract m_V0.
                rb[1] -= m_V0;

                rb[3] += x[0];       // lambda
                for (k = 0; k < m_nsp; k++) {
                    rb[4+k] += x[0]*m_yin[k];
                }
            }

            // right inlet.
            else {
                int boffset = m_flow->nComponents();
                xb = x - boffset;
                rb = r - boffset;
                rb[1] -= m_V0;
                rb[2] -= x[1]; // T
                xb[0] += x[0]; // u
                for (k = 0; k < m_nsp; k++) 
                    rb[4+k] += x[0]*m_yin[k];
            }                
        }

        virtual void save(XML_Node& o, doublereal* soln) {
            doublereal* s = soln + loc();
            XML_Node& inlt = o.addChild("inlet");
            for (int k = 0; k < 2; k++) {
                ctml::addFloat(inlt, componentName(k), s[k], "", "",0.0, 1.0);
            }
        }
                
    protected:

        int m_ilr;
        doublereal m_mdot, m_temp, m_V0;
        StFlow *m_flow;
        int m_nsp;
        vector_fp m_yin;
        string m_xstr;
    };


    class Symm1D : public Bdry1D {

    public:

        Symm1D(int ilr = 1) {
            m_type = cSymmType; 
            m_flow = 0;
            m_ilr = ilr;
        }
        virtual ~Symm1D(){}

        virtual string componentName(int n) const { 
            switch (n) {
            case 0: return "dummy"; break;
            default: return "<unknown>";
            }
        }

        virtual void init() {
            if (m_index < 0) {
                throw CanteraError("Symm1D", 
                    "install in container before calling init.");
            }
            resize(1,1);

            // set bounds
            const doublereal lower = -1.0e5;
            const doublereal upper =  1.0e5;
            setBounds(1, &lower, 1, &upper);

            // set tolerances
            doublereal rtol = 1e-4;
            doublereal atol = 1.e-5;
            setTolerances(1, &rtol, 1, &atol);

            if (m_index > 0) {
                Resid1D& r = container().domain(m_index-1);
                if (r.domainType() == cFlowType) {
                    m_ilr = -1;
                    m_flow = (StFlow*)&r;
                }
                else 
                    throw CanteraError("Symm1D::init",
                        "Symmetry planes can only be connected to flow domains.");
            }
            else {
                if (container().nDomains() > 1) {
                    Resid1D& r = container().domain(1);
                    if (r.domainType() == cFlowType) {
                        m_ilr = 1;
                        m_flow = (StFlow*)&r;
                    }
                    else 
                        throw CanteraError("Symm1D::init",
                            "Symmetry planes can only be connected to flow domains.");
                }
                else 
                    throw CanteraError("Symm1D::init",
                        "A symmetry plane must be connected to a flow domain.");
            }
            m_nsp = m_flow->nComponents() - 4;
        }


        virtual void eval(int jg, doublereal* xg, doublereal* rg, 
            integer* diagg, doublereal rdt) {
            int k;
            if (jg >= 0 && (jg < firstPoint() - 2 || jg > lastPoint() + 2)) return;

            // start of local part of global arrays
            doublereal* x = xg + loc();
            doublereal* r = rg + loc();
            integer* diag = diagg + loc();
            doublereal *xb, *rb;
            // integer *db = diag + loc();

            r[0] = x[0];
            diag[0] = 0;
            int nc = m_flow->nComponents();

            if (m_ilr == 1) {
                xb = x + 1;
                rb = r + 1;
                rb[0] = xb[0];
                rb[1] = xb[1] - xb[1 + nc];
                if (m_flow->doEnergy(0)) {
                    rb[2] = xb[2] - xb[2 + nc];
                }
                rb[3] = xb[3] - xb[3 + nc];
                for (k = 0; k < m_nsp; k++) {
                    rb[4+k] = xb[4+k] - xb[4+k+nc];
                }
            }
            else {
                xb = x - nc;
                rb = r - nc;
                rb[0] = xb[0];
                rb[1] = xb[1] - xb[1 - nc];
                //               if (m_flow->doEnergy(0)) {
                rb[2] = xb[2] - xb[2 - nc];
                    // }
                rb[3] = xb[3] - xb[3 - nc];
                for (k = 0; k < m_nsp; k++) {
                    rb[4+k] = xb[4+k] - xb[4+k-nc];
                }
            }                
        }

        virtual void save(XML_Node& o, doublereal* soln) {
	    // doublereal* s = soln + loc();
            // XML_Node& symm = o.addChild("symmetry_plane");
            (void) o.addChild("symmetry_plane");
        }

    protected:

        int m_ilr;
        StFlow *m_flow;
        int m_nsp;
    };



    /////////////////////////////////////////////////////////////
    //
    // surf1D
    //
    ////////////////////////////////////////////////////////////

    class Surf1D : public Bdry1D {

    public:

        Surf1D(int ilr = 1) {
            m_type = cSurfType; 
            m_flow = 0;
            m_ilr = ilr;
        }
        virtual ~Surf1D(){}

        /// Set the surface temperature
        virtual void setTemperature(doublereal t) {
            m_temp = t;
            needJacUpdate();
        }

        /// Temperature [K].
        doublereal temperature() {return m_temp;}

        virtual string componentName(int n) const { 
            switch (n) {
            case 0: return "dummy"; break;
            case 1: return "temperature"; break;
            default: return "<unknown>";
            }
        }

        virtual void init() {
            if (m_index < 0) {
                throw CanteraError("Surf1D", 
                    "install in container before calling init.");
            }
            resize(2,1);

            // set bounds
            const doublereal lower[2] = {-1.0e5, 200.0};
            const doublereal upper[2] = {1.0e5, 1.e5};
            setBounds(2, lower, 2, upper);

            // set tolerances
            vector_fp rtol(2, 1e-4);
            vector_fp atol(2, 1.e-5);
            setTolerances(2, rtol.begin(), 2, atol.begin());

            if (m_index > 0) {
                Resid1D& r = container().domain(m_index-1);
                if (r.domainType() == cFlowType) {
                    m_ilr = -1;
                    m_flow = (StFlow*)&r;
                }
                else 
                    throw CanteraError("Surf1D::init",
                        "Surface domains can only be connected to a flow domain.");
            }
            else {
                if (container().nDomains() > 1) {
                    Resid1D& r = container().domain(1);
                    if (r.domainType() == cFlowType) {
                        m_ilr = 1;
                        m_flow = (StFlow*)&r;
                    }
                    else 
                        throw CanteraError("Surf1D::init",
                            "A surface domain can only be connected to a flow domain.");
                }
                else 
                    throw CanteraError("Surf1D::init",
                        "A surface domain must be connected to a flow domain.");
            }
            m_nsp = m_flow->nComponents() - 4;
        }


        virtual void eval(int jg, doublereal* xg, doublereal* rg, 
            integer* diagg, doublereal rdt) {
            int k;
            if (jg >= 0 && (jg < firstPoint() - 2 || jg > lastPoint() + 2)) return;

            // start of local part of global arrays
            doublereal* x = xg + loc();
            doublereal* r = rg + loc();
            integer* diag = diagg + loc();
            doublereal *xb, *rb;
            //integer *db = diag + loc();

            r[0] = x[0];
            r[1] = m_temp - x[1];
            diag[0] = 0;
            diag[1] = 0;
            int nc = m_flow->nComponents();

            if (m_ilr == 1) {
                xb = x + 2;
                rb = r + 2;
                rb[0] = xb[0];
                rb[1] = xb[1];   // T
                //                if (m_flow->doEnergy(0)) {
                rb[2] = xb[2] - x[1];   // T
                    //}
                rb[3] = xb[3] - xb[3 + nc];   // lambda
                for (k = 0; k < m_nsp; k++) {
                    rb[4+k] += xb[4+k] - xb[4+k+nc];
                }
            }
            else {
                xb = x - nc;
                rb = r - nc;
                rb[0] = xb[0];
                rb[1] = xb[1];   // T
                //                if (m_flow->doEnergy(0)) {
                rb[2] = xb[2] - x[1];   // T
                    //}
                rb[3] = xb[3] - xb[3 - nc];   // lambda
                for (k = 0; k < m_nsp; k++) {
                    rb[4+k] += xb[4+k] - xb[4+k-nc];
                }
            }
        }

        virtual void save(XML_Node& o, doublereal* soln) {
            doublereal* s = soln + loc();
            XML_Node& srf = o.addChild("surface");
            for (int k = 1; k < 2; k++) {
                ctml::addFloat(srf, componentName(k), s[k], "", "",0.0, 1.0);
            }
        }
                
    protected:

        int m_ilr;
        doublereal m_temp;
        StFlow *m_flow;
        int m_nsp;
    };


}

#endif
