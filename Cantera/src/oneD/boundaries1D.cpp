/**
 * @file boundaries1D.cpp
 */

/*
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2002-3  California Institute of Technology


// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "Inlet1D.h"

namespace Cantera {

    Bdry1D::Bdry1D() : Domain1D(1, 1, 0.0), 
                       m_flow_left(0), m_flow_right(0),
                       m_ilr(0), m_left_nv(0), m_right_nv(0),
                       m_left_loc(0), m_right_loc(0),
                       m_left_points(0), m_nv(0), 
                       m_left_nsp(0), m_right_nsp(0),
                       m_sp_left(0), m_sp_right(0),
                       m_start_left(0), m_start_right(0),
                       m_phase_left(0), m_phase_right(0), m_temp(0.0), m_mdot(0.0) {
        m_type = cConnectorType;
    }


    void Bdry1D::
    _init(int n) {
        if (m_index < 0) {
            throw CanteraError("Bdry1D", 
                "install in container before calling init.");
        }
        resize(n,1);

        m_left_nsp = 0;
        m_right_nsp = 0;

        // check for left and right flow objects
        if (m_index > 0) {
            Domain1D& r = container().domain(m_index-1);
            if (r.domainType() == cFlowType) {
                m_flow_left = (StFlow*)&r;
                m_left_nv = m_flow_left->nComponents();
                m_left_points = m_flow_left->nPoints();
                m_left_loc = container().start(m_index-1);
                m_left_nsp = m_left_nv - 4;
                m_phase_left = &m_flow_left->phase();
            }
            else 
                throw CanteraError("Bdry1D::init",
                    "Boundary domains can only be "
                    "connected on the left to flow domains, not type "+int2str(r.domainType()) 
                    + " domains.");
        }
        
        if (m_index < container().nDomains() - 1) {
            Domain1D& r = container().domain(m_index+1);
            if (r.domainType() == cFlowType) {
                m_flow_right = (StFlow*)&r;
                m_right_nv = m_flow_right->nComponents();
                m_right_loc = container().start(m_index+1);
                m_right_nsp = m_right_nv - 4;
                m_phase_right = &m_flow_right->phase();
            }
            else 
                throw CanteraError("Bdry1D::init",
                    "Boundary domains can only be "
                    "connected on the right to flow domains, not type "+int2str(r.domainType()) 
                    + " domains.");
        }
    }




    //----------------------------------------------------------
    //
    //   Inlet1D methods
    //
    //----------------------------------------------------------


    void Inlet1D::
    setMoleFractions(string xin) {
        m_xstr = xin;
        if (m_flow) {
            m_flow->phase().setMoleFractionsByName(xin);
            m_flow->phase().getMassFractions(m_yin.begin());
            needJacUpdate();
        }
    }

    void Inlet1D::
    setMoleFractions(doublereal* xin) {
        if (m_flow) {
            m_flow->phase().setMoleFractions(xin);
            m_flow->phase().getMassFractions(m_yin.begin());
            needJacUpdate();
        }
    }
 
    string Inlet1D::
    componentName(int n) const { 
        switch (n) {
        case 0: return "mdot"; break;
        case 1: return "temperature"; break;
        default: return "unknown";
        }
    }

    void Inlet1D::
    init() {

        _init(2);

        // set bounds (mdot, T)
        const doublereal lower[2] = {-1.0e5, 200.0};
        const doublereal upper[2] = {1.0e5, 1.e5};
        setBounds(2, lower, 2, upper);

        // set tolerances
        vector_fp rtol(2, 1e-4);
        vector_fp atol(2, 1.e-5);
        setTolerances(2, rtol.begin(), 2, atol.begin());

        // if a flow domain is present on the left, then this must be
        // a right inlet. Note that an inlet object can only be a
        // terminal object.
        if (m_flow_left) {
            m_ilr = RightInlet;
            m_flow = m_flow_left;
        }
        else if (m_flow_right) {
            m_ilr = LeftInlet;
            m_flow = m_flow_right;
        }
        else {
            throw CanteraError("Inlet1D::init","no flow!");
        }

        // components = u, V, T, lambda, + mass fractions
        m_nsp = m_flow->nComponents() - 4;
        m_yin.resize(m_nsp, 0.0);
        if (m_xstr != "") 
            setMoleFractions(m_xstr);
        else
            m_yin[0] = 1.0;
    }


    void Inlet1D::
    eval(int jg, doublereal* xg, doublereal* rg, 
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

            rb[2] = xb[2] - x[1];   // T

            // spreading rate. Flow domain sets this to V(0),
            // so for finite spreading rate subtract m_V0.
            rb[1] -= m_V0;

            rb[3] += x[0];       // lambda
            for (k = 1; k < m_nsp; k++) {
                if (m_flow->doSpecies(k)) {
                    rb[4+k] += x[0]*m_yin[k];
                }
            }
        }

        // right inlet.
        else {
            int boffset = m_flow->nComponents();
            xb = x - boffset;
            rb = r - boffset;
            rb[1] -= m_V0;
            rb[2] = xb[2] - x[1]; // T
            rb[0] += x[0];        // u
            for (k = 1; k < m_nsp; k++) {
                if (m_flow->doSpecies(k)) {
                    //                    rb[4+k] += x[0]*(-xb[4+k] + m_yin[k]);
                    rb[4+k] += x[0]*(m_yin[k]);
                }
            }
        }                
    }

    void Inlet1D::
    save(XML_Node& o, doublereal* soln) {
        doublereal* s = soln + loc();
        //XML_Node& inlt = o.addChild("inlet");
        XML_Node& inlt = o.addChild("domain");
        inlt.addAttribute("id",id());
        inlt.addAttribute("points",1);
        inlt.addAttribute("type","inlet");
        inlt.addAttribute("components",nComponents());
        for (int k = 0; k < nComponents(); k++) {
            ctml::addFloat(inlt, componentName(k), s[k], "", "",lowerBound(k), upperBound(k));
        }
    }
             
    void Inlet1D::
    restore(XML_Node& dom, doublereal* soln) {
        //map<string, double> x;
        //getFloats(dom, x);
        soln[0] = getFloat(dom, "mdot", "massflowrate"); // x["mdot"];
        soln[1] = getFloat(dom, "temperature", "temperature"); // x["temperature"];
        resize(2,1);
    }




    //--------------------------------------------------   
    //      Empty1D
    //--------------------------------------------------

    string Empty1D::componentName(int n) const { 
        switch (n) {
        case 0: return "dummy"; break;
        default: return "<unknown>";
        }
    }

    void Empty1D::
    init() { //_init(1); 
       // set bounds (T)
        const doublereal lower = -1.0;
        const doublereal upper = 1.0;
        setBounds(1, &lower, 1, &upper);

        // set tolerances
        const doublereal rtol = 1e-4;
        const doublereal atol = 1.e-4;
        setTolerances(1, &rtol, 1, &atol);
    }

    void Empty1D::
    eval(int jg, doublereal* xg, doublereal* rg, 
        integer* diagg, doublereal rdt) {
        if (jg >= 0 && (jg < firstPoint() - 2 || jg > lastPoint() + 2)) return;

        // start of local part of global arrays
        doublereal* x = xg + loc();
        doublereal* r = rg + loc();
        integer* diag = diagg + loc();
//        integer *db;

        r[0] = x[0];
        diag[0] = 0;
    }

    void Empty1D::
    save(XML_Node& o, doublereal* soln) {
        XML_Node& symm = o.addChild("domain");
        symm.addAttribute("id",id());
        symm.addAttribute("points",1);
        symm.addAttribute("type","empty");
        symm.addAttribute("components",nComponents());
    }

    void Empty1D::
    restore(XML_Node& dom, doublereal* soln) {
        resize(1,1);
    }



    //--------------------------------------------------   
    //      Symm1D
    //--------------------------------------------------

    string Symm1D::componentName(int n) const { 
        switch (n) {
        case 0: return "dummy"; break;
        default: return "<unknown>";
        }
    }

    void Symm1D::
    init() { _init(1); 
       // set bounds (T)
        const doublereal lower = -1.0;
        const doublereal upper = 1.0;
        setBounds(1, &lower, 1, &upper);

        // set tolerances
        const doublereal rtol = 1e-4;
        const doublereal atol = 1.e-4;
        setTolerances(1, &rtol, 1, &atol);
    }

    void Symm1D::
    eval(int jg, doublereal* xg, doublereal* rg, 
        integer* diagg, doublereal rdt) {
        if (jg >= 0 && (jg < firstPoint() - 2 || jg > lastPoint() + 2)) return;

        // start of local part of global arrays
        doublereal* x = xg + loc();
        doublereal* r = rg + loc();
        integer* diag = diagg + loc();
        doublereal *xb, *rb;
        integer *db;

        r[0] = x[0];
        diag[0] = 0;
        int nc;

        if (m_flow_right) {
            nc = m_flow_right->nComponents();
            xb = x + 1;
            rb = r + 1;
            db = diag + 1;
            db[1] = 0;
            db[2] = 0;
            rb[1] = xb[1] - xb[1 + nc];      // zero dV/dz
            rb[2] = xb[2] - xb[2 + nc];      // zero dT/dz
        }

        if (m_flow_left) {
            nc = m_flow_left->nComponents();
            xb = x - nc;
            rb = r - nc;
            db = diag - nc;
            db[1] = 0;
            db[2] = 0;
            rb[1] = xb[1] - xb[1 - nc];      // zero dV/dz
            rb[2] = xb[2] - xb[2 - nc];      // zero dT/dz
        }
    }


    void Symm1D::
    save(XML_Node& o, doublereal* soln) {
        XML_Node& symm = o.addChild("domain");
        symm.addAttribute("id",id());
        symm.addAttribute("points",1);
        symm.addAttribute("type","outlet");
        symm.addAttribute("components",nComponents());
    }

    void Symm1D::
    restore(XML_Node& dom, doublereal* soln) {
        resize(1,1);
    }


    //--------------------------------------------------   
    //      Outlet1D
    //--------------------------------------------------

    string Outlet1D::componentName(int n) const { 
        switch (n) {
        case 0: return "dummy"; break;
        default: return "<unknown>";
        }
    }

    void Outlet1D::
    init() { 
        _init(1); 
       // set bounds (T)
        const doublereal lower = -1.0;
        const doublereal upper = 1.0;
        setBounds(1, &lower, 1, &upper);

        // set tolerances
        const doublereal rtol = 1e-4;
        const doublereal atol = 1.e-4;
        setTolerances(1, &rtol, 1, &atol);
    }


    void Outlet1D::
    eval(int jg, doublereal* xg, doublereal* rg, 
        integer* diagg, doublereal rdt) {
        if (jg >= 0 && (jg < firstPoint() - 2 || jg > lastPoint() + 2)) return;

        // start of local part of global arrays
        doublereal* x = xg + loc();
        doublereal* r = rg + loc();
        integer* diag = diagg + loc();
        doublereal *xb, *rb;
        integer *db;

        r[0] = x[0];
        diag[0] = 0;
        int nc, k;

        if (m_flow_right) {
            nc = m_flow_right->nComponents();
            xb = x + 1;
            rb = r + 1;
            db = diag + 1;
            rb[0] = xb[3];      
            rb[2] = xb[2] - xb[2 + nc];
            for (k = 4; k < nc; k++) {
                //if (m_flow_right->doSpecies(k-4)) {
                rb[k] = xb[k] - xb[k + nc];
                //}
            }
        }

        if (m_flow_left) {
            nc = m_flow_left->nComponents();
            xb = x - nc;
            rb = r - nc;
            db = diag - nc;

            // zero Lambda
            rb[0] = xb[3];               // zero Lambda
            rb[2] = xb[2] - xb[2 - nc];  // zero T gradient
            for (k = 5; k < nc; k++) {
                rb[k] = xb[k] - xb[k - nc]; // zero mass fraction gradient
            }
        }
    }


    void Outlet1D::
    save(XML_Node& o, doublereal* soln) {
        XML_Node& outlt = o.addChild("domain");
        outlt.addAttribute("id",id());
        outlt.addAttribute("points",1);
        outlt.addAttribute("type","outlet");
        outlt.addAttribute("components",nComponents());
    }

    void Outlet1D::
    restore(XML_Node& dom, doublereal* soln) {
        resize(1,1);
    }




    //--------------------------------------------------   
    //      OutletRes1D
    //--------------------------------------------------


    void OutletRes1D::
    setMoleFractions(string xres) {
        m_xstr = xres;
        if (m_flow) {
            m_flow->phase().setMoleFractionsByName(xres);
            m_flow->phase().getMassFractions(m_yres.begin());
            needJacUpdate();
        }
    }

    void OutletRes1D::
    setMoleFractions(doublereal* xres) {
        if (m_flow) {
            m_flow->phase().setMoleFractions(xres);
            m_flow->phase().getMassFractions(m_yres.begin());
            needJacUpdate();
        }
    }
 
    string OutletRes1D::componentName(int n) const { 
        switch (n) {
        case 0: return "dummy"; break;
        default: return "<unknown>";
        }
    }

    void OutletRes1D::
    init() { 
        _init(1); 
       // set bounds (dummy)
        const doublereal lower = -1.0;
        const doublereal upper = 1.0;
        setBounds(1, &lower, 1, &upper);

        // set tolerances
        const doublereal rtol = 1e-4;
        const doublereal atol = 1.e-4;
        setTolerances(1, &rtol, 1, &atol);

        if (m_flow_left) {
            m_flow = m_flow_left;
        }
        else if (m_flow_right) {
            m_flow = m_flow_right;
        }
        else {
            throw CanteraError("OutletRes1D::init","no flow!");
        }

        m_nsp = m_flow->nComponents() - 4;
        m_yres.resize(m_nsp, 0.0);
        if (m_xstr != "") 
            setMoleFractions(m_xstr);
        else
            m_yres[0] = 1.0;
    }


    void OutletRes1D::
    eval(int jg, doublereal* xg, doublereal* rg, 
        integer* diagg, doublereal rdt) {

        if (jg >= 0 && (jg < firstPoint() - 2 || jg > lastPoint() + 2)) return;

        // start of local part of global arrays
        doublereal* x = xg + loc();
        doublereal* r = rg + loc();
        integer* diag = diagg + loc();
        doublereal *xb, *rb;
        integer *db;

        // drive dummy component to zero
        r[0] = x[0];
        diag[0] = 0;
        int nc, k;

        if (m_flow_right) {
            nc = m_flow_right->nComponents();
            xb = x + 1;
            rb = r + 1;
            db = diag + 1;

            // this seems wrong...
            // zero Lambda
            rb[0] = xb[3];
      
            // zero gradient for T
            rb[2] = xb[2] - xb[2 + nc];

            // specified mass fractions
            for (k = 4; k < nc; k++) {
                rb[k] = xb[k] - m_yres[k-4];
            }
        }

        if (m_flow_left) {

            nc = m_flow_left->nComponents();
            xb = x - nc;
            rb = r - nc;
            db = diag - nc;

            rb[0] = xb[3];                     // zero Lambda
            rb[2] = xb[2] - m_temp; //xb[2] - xb[2 - nc];        // zero dT/dz
            for (k = 5; k < nc; k++) {
                rb[k] = xb[k] - m_yres[k-4];     // fixed Y
            }
        }
    }


    void OutletRes1D::
    save(XML_Node& o, doublereal* soln) {
        XML_Node& outlt = o.addChild("domain");
        outlt.addAttribute("id",id());
        outlt.addAttribute("points",1);
        outlt.addAttribute("type","outletres");
        outlt.addAttribute("components",nComponents());
    }

    void OutletRes1D::
    restore(XML_Node& dom, doublereal* soln) {
        resize(1,1);
    }


    //-----------------------------------------------------------
    //
    //  Surf1D
    //
    //-----------------------------------------------------------



    string Surf1D::componentName(int n) const { 
        switch (n) {
        case 0: return "temperature"; break;
        default: return "<unknown>";
        }
    }

    void Surf1D::
    init() { 
        _init(1); 
       // set bounds (T)
        const doublereal lower = 200.0;
        const doublereal upper = 1.e5;
        setBounds(1, &lower, 1, &upper);

        // set tolerances
        const doublereal rtol = 1e-4;
        const doublereal atol = 1.e-4;
        setTolerances(1, &rtol, 1, &atol);
    }


    void Surf1D::
    eval(int jg, doublereal* xg, doublereal* rg, 
        integer* diagg, doublereal rdt) {
        if (jg >= 0 && (jg < firstPoint() - 2 || jg > lastPoint() + 2)) return;

        // start of local part of global arrays
        doublereal* x = xg + loc();
        doublereal* r = rg + loc();
        integer* diag = diagg + loc();
        doublereal *xb, *rb;

        r[0] = x[0] - m_temp;
        diag[0] = 0;
        int nc;

        if (m_flow_right) {
            rb = r + 1;
            xb = x + 1;
            rb[2] = xb[2] - x[0];            // specified T
        }

        if (m_flow_left) {
            nc = m_flow_left->nComponents();
            rb = r - nc;
            xb = x - nc;
            rb[2] = xb[2] - x[0];            // specified T
        }
    }

    void Surf1D::
    save(XML_Node& o, doublereal* soln) {
        doublereal* s = soln + loc();
        //XML_Node& inlt = o.addChild("inlet");
        XML_Node& inlt = o.addChild("domain");
        inlt.addAttribute("id",id());
        inlt.addAttribute("points",1);
        inlt.addAttribute("type","surface");
        inlt.addAttribute("components",nComponents());
        for (int k = 0; k < nComponents(); k++) {
            ctml::addFloat(inlt, componentName(k), s[k], "", "",0.0, 1.0);
        }
    }

    void Surf1D::
    restore(XML_Node& dom, doublereal* soln) {
        map<string, double> x;
        getFloats(dom, x);
        soln[0] = x["temperature"];
        resize(1,1);
    }




    //-----------------------------------------------------------
    //
    //  ReactingSurf1D
    //
    //-----------------------------------------------------------



    string ReactingSurf1D::componentName(int n) const { 
        if (n == 0) return "temperature";
        else if (n < m_nsp + 1) 
            return m_sphase->speciesName(n-1);
        else
            return "<unknown>";
    }

    void ReactingSurf1D::
    init() {
        m_nv = m_nsp + 1;
        _init(m_nsp+1); 
        m_fixed_cov.resize(m_nsp, 0.0);
        m_fixed_cov[0] = 1.0;
        int nt = m_kin->nTotalSpecies();
        m_work.resize(nt, 0.0);

       // set bounds 
        vector_fp lower(m_nv), upper(m_nv);
        lower[0] = 200.0;
        upper[0] = 1.e5;
        int n;
        for (n = 0; n < m_nsp; n++) {
            lower[n+1] = -1.0e-5;
            upper[n+1] = 2.0;
        }
        setBounds(m_nv, lower.begin(), m_nv, upper.begin());
        vector_fp rtol(m_nv), atol(m_nv);
        for (n = 0; n < m_nv; n++) {
            rtol[n] = 1.0e-5;
            atol[n] = 1.0e-9;
        }
        atol[0] = 1.0e-4;
        setTolerances(m_nv, rtol.begin(), m_nv, atol.begin());
    }


    void ReactingSurf1D::
    eval(int jg, doublereal* xg, doublereal* rg, 
        integer* diagg, doublereal rdt) {
        if (jg >= 0 && (jg < firstPoint() - 2 || jg > lastPoint() + 2)) return;

        // start of local part of global arrays
        doublereal* x = xg + loc();
        doublereal* r = rg + loc();
        integer* diag = diagg + loc();
        doublereal *xb, *rb;

        // specified surface temp
        r[0] = x[0] - m_temp;

        // set the coverages
        doublereal sum = 0.0;
        int k;
        for (k = 0; k < m_nsp; k++) {
            m_work[k] = x[k+1];
            sum += x[k+1];
        }
        m_sphase->setTemperature(x[0]);
        m_sphase->setCoverages(m_work.begin());
        //m_kin->advanceCoverages(1.0);
        //m_sphase->getCoverages(m_fixed_cov.begin());

        // set the left gas state to the adjacent point

        int leftloc = 0, rightloc = 0;
        int pnt = 0;

        if (m_flow_left) {
            leftloc = m_flow_left->loc();
            pnt = m_flow_left->nPoints() - 1;
            m_flow_left->setGas(xg + leftloc, pnt);
        }

        if (m_flow_right) {
            rightloc = m_flow_right->loc();
            m_flow_right->setGas(xg + rightloc, 0);
        }

        m_kin->getNetProductionRates(m_work.begin());
        doublereal rs0 = 1.0/m_sphase->siteDensity();
            
        //scale(m_work.begin(), m_work.end(), m_work.begin(), m_mult[0]);
        
        //        bool enabled = true;
        int ioffset = m_kin->kineticsSpeciesIndex(0, m_surfindex);

        if (m_enabled) {
            doublereal maxx = -1.0;
            int imx = -1;
            for (k = 0; k < m_nsp; k++) {
                r[k+1] = m_work[k + ioffset] * m_sphase->size(k) * rs0;
                r[k+1] -= rdt*(x[k+1] - prevSoln(k+1,0));
                diag[k+1] = 1;
                if (x[k+1] > maxx) {
                    maxx = x[k+1];
                    imx = k+1;
                }
            }
            r[1] = 1.0 - sum;
            diag[1] = 0;
        }
        else {
            for (k = 0; k < m_nsp; k++) {
                r[k+1] = x[k+1] - m_fixed_cov[k];
                diag[k+1] = 0;
            }
        }
        
        if (m_flow_right) {
            rb = r + 1;
            xb = x + 1;
            rb[2] = xb[2] - x[0];            // specified T
        }
        int nc;
        if (m_flow_left) {
            nc = m_flow_left->nComponents();
            const doublereal* mwleft = m_phase_left->molecularWeights().begin();
            rb =r - nc;
            xb = x - nc;
            rb[2] = xb[2] - x[0];            // specified T
            for (int nl = 1; nl < m_left_nsp; nl++) {
                rb[4+nl] += m_work[nl]*mwleft[nl];
            }
        }

        // gas-phase residuals
//         doublereal rho;
//         if (m_flow_left) {
//             rho = m_phase_left->density();
            //            doublereal rdz = 2.0/
            //                 (m_flow_left->z(m_left_points-1) - 
            //                     m_flow_left->z(m_left_points - 2));


//             for (k = 0; k < m_left_nsp; k++) 
//                 m_work[k + m_start_left] *= m_molwt_left[k];

//             int ileft = loc() - m_left_nv;

//             // if the energy equation is enabled at this point,
//             // set the gas temperature to the surface temperature
//             if (m_flow_left->doEnergy(pnt)) {
//                 rg[ileft + 2] = xg[ileft + 2] - m_sphase->temperature();
//             }

//             for (k = 1; k < m_left_nsp; k++) {
//                 if (enabled && m_flow_left->doSpecies(k)) {
//                     rg[ileft + 4 + k]  += m_work[k + m_start_left];
//                     //+= rdz*m_work[k + m_sp_left]/rho;
                    
//                 }
//             }
//         }

//         if (m_flow_right) {
//             for (k = 0; k < m_right_nsp; k++) 
//                 m_work[k + m_start_right] *= m_molwt_right[k];
            
//             int iright = loc() + m_nsp;
//             rg[iright + 2] -= m_sphase->temperature();
//             //r[iright + 3] = x[iright];
//             for (k = 0; k < m_right_nsp; k++) {
//                 rg[iright + 4 + k] -= m_work[k + m_start_right];
//             }
//         }
    

//         diag[0] = 0;
//         int nc;

//         if (m_flow_right) {
//             rb = r + 1;
//             xb = x + 1;
//             rb[2] = xb[2] - x[0];            // specified T
//         }

//         if (m_flow_left) {
//             nc = m_flow_left->nComponents();
//             rb = r - nc;
//             xb = x - nc;
//             rb[2] = xb[2] - x[0];            // specified T
//         }
    }

    void ReactingSurf1D::
    save(XML_Node& o, doublereal* soln) {
        doublereal* s = soln + loc();
        //XML_Node& inlt = o.addChild("inlet");
        XML_Node& inlt = o.addChild("domain");
        inlt.addAttribute("id",id());
        inlt.addAttribute("points",1);
        inlt.addAttribute("type","surface");
        inlt.addAttribute("components",nComponents());
        for (int k = 0; k < nComponents(); k++) {
            ctml::addFloat(inlt, componentName(k), s[k], "", "",0.0, 1.0);
        }
    }

    void ReactingSurf1D::
    restore(XML_Node& dom, doublereal* soln) {
        map<string, double> x;
        getFloats(dom, x);
        soln[0] = x["temperature"];
        resize(1,1);
    }
}


    /////////////////////////////////////////////////////////////
    //
    // surf1D
    //
    ////////////////////////////////////////////////////////////


// string ChemSurf1D::
//     componentName(int n) const {
//         /// @todo why dummy?
//         switch (n) {
//         case 0: return "dummy"; break;
//         case 1: return "temperature"; break;
//         default: return "<unknown>";
//         }
//     }

//     // Set the kinetics manager for the surface.
//     void ChemSurf1D::
//     setKinetics(InterfaceKinetics* kin) {
//         m_kin = kin;
//         int np = kin->nPhases();
//         m_sphase = 0;
//         for (int n = 0; n < np; n++) {
//             if (kin->phase(n).eosType() == cSurf) {
//                 m_sphase = (SurfPhase*)&m_kin->phase(n);
//                 m_nsurf = n;
//             }
//             else {
//                 m_bulk.push_back(&kin->phase(n));
//                 m_nbulk.push_back(n);
//             }
//         }
//         if (!m_sphase) 
//             throw CanteraError("setKinetics","no surface phase defined");
        
//         m_nsp = m_sphase->nSpecies();
//         resize(m_nsp,1);
//         if (m_bulk.size() == 1) {
//             m_bulk.push_back(0);
//         }
//     }

// void ChemSurf1D::
//     init() {
//         cout << "ChemSurf1D::init" << endl;
//         if (m_index < 0) {
//             throw CanteraError("Surf1D", 
//                 "install in container before calling init.");
//         }
//         resize(m_nsp,1);
//         m_mult.resize(m_nsp, 1.0);
//         m_do_surf_species.resize(m_nsp, true);
//         if (!m_sphase) m_do_surf_species[0] = false;
//         m_fixed_cov.resize(m_nsp, 1.0/m_nsp);

//         // set bounds
//         vector_fp lower(m_nsp, -1.e-3);
//         vector_fp upper(m_nsp, 1.0);
//         setBounds(m_nsp, lower.begin(), m_nsp, upper.begin());

//         // set tolerances
//         vector_fp rtol(m_nsp, 1e-4);
//         vector_fp atol(m_nsp, 1.e-10);
//         setTolerances(m_nsp, rtol.begin(), m_nsp, atol.begin());

//         m_left_nsp = 0;
//         m_right_nsp = 0;

//         // check for left and right flow objects
//         if (m_index > 0) {
//             Domain1D& r = container().domain(m_index-1);
//             if (r.domainType() == cFlowType) {
//                 m_flow_left = (StFlow*)&r;
//                 m_left_nv = m_flow_left->nComponents();
//                 m_left_points = m_flow_left->nPoints();
//                 m_left_loc = container().start(m_index-1);
//                 m_left_nsp = m_left_nv - 4;
//                 m_phase_left = &m_flow_left->phase();
//                 m_molwt_left = m_phase_left->molecularWeights().begin();
//                 if (m_phase_left == m_bulk[0]) 
//                     m_start_left = m_kin->start(m_nbulk[0]);
//                 else if (m_phase_left == m_bulk[1]) 
//                     m_start_left = m_kin->start(m_nbulk[1]);
//                 else 
//                     throw CanteraError("ChemSurf1D::init",
//                         "left gas does not match one in surface mechanism");
//             }
//             else 
//                 throw CanteraError("ChemSurf1D::init",
//                     "Surface domains can only be "
//                     "connected to flow domains.");
//         }
        
//         if (m_index < container().nDomains() - 1) {
//             Domain1D& r = container().domain(m_index+1);
//             if (r.domainType() == cFlowType) {
//                 m_flow_right = (StFlow*)&r;
//                 m_right_nv = m_flow_right->nComponents();
//                 m_right_loc = container().start(m_index+1);
//                 m_right_nsp = m_right_nv - 4;
//                 m_phase_right = &m_flow_right->phase();
//                 m_molwt_right = m_phase_right->molecularWeights().begin();
//                 if (m_phase_right == m_bulk[0]) 
//                     m_start_right = m_kin->start(m_nbulk[0]);
//                 else if (m_phase_right == m_bulk[1]) 
//                     m_start_right = m_kin->start(m_nbulk[1]);
//                 else 
//                     throw CanteraError("ChemSurf1D::init",
//                         "right gas does not match one in surface mechanism");
//             }
//             else 
//                 throw CanteraError("ChemSurf1D::init",
//                     "Surface domains can only be "
//                     "connected to flow domains.");
//         }
//         m_work.resize(m_kin->nTotalSpecies());
//     }


// void ChemSurf1D::eval(int jg, doublereal* xg, doublereal* rg, 
//         integer* diagg, doublereal rdt) {
//         int k;

//         // if computing a Jacobian (jg > 0), and the global point is
//         // outside the points the surface can influence, then skip
//         // evaluating the residual            
//         if (jg >= 0 && (jg < firstPoint() - 2 
//                 || jg > lastPoint() + 2)) return;

//         // start of local part of global arrays
//         doublereal* x = xg + loc();
//         doublereal* r = rg + loc();
//         integer* diag = diagg + loc();

//         // set the coverages
//         doublereal sum = 0.0;
//         for (k = 0; k < m_nsp; k++) {
//             m_work[k] = x[k];
//             sum += x[k];
//         }
//         m_sphase->setCoverages(m_work.begin());

//         // set the left gas state to the adjacent point

//         int leftloc = 0, rightloc = 0;
//         int pnt = 0;

//         if (m_flow_left) {
//             leftloc = m_flow_left->loc();
//             pnt = m_flow_left->nPoints() - 1;
//             m_flow_left->setGas(xg + leftloc, pnt);
//         }

//         if (m_flow_right) {
//             rightloc = m_flow_right->loc();
//             m_flow_right->setGas(xg + rightloc, 0);
//         }

//         m_kin->getNetProductionRates(m_work.begin());
//         doublereal rs0 = 1.0/m_sphase->siteDensity();
            
//         scale(m_work.begin(), m_work.end(), m_work.begin(), m_mult[0]);
        
//         bool enabled = true;
//         int ioffset = m_kin->start(m_nsurf); // m_left_nsp + m_right_nsp;
//         doublereal maxx = -1.0;
//         int imx = -1;
//         for (k = 0; k < m_nsp; k++) {
//             r[k] = m_work[k + ioffset] * m_sphase->size(k) * rs0;
//             r[k] -= rdt*(x[k] - prevSoln(k,0));
//             diag[k] = 1;
//             if (x[k] > maxx) {
//                 maxx = x[k];
//                 imx = k;
//             }
//             if (!m_do_surf_species[k]) {
//                 r[k] = x[k] - m_fixed_cov[k];
//                 diag[k] = 0;
//                 enabled = false;
//             }
//         }
//         if (enabled) {
//             r[imx] = 1.0 - sum;
//             diag[imx] = 0;
//         }
        
//         // gas-phase residuals
//         doublereal rho;
//         if (m_flow_left) {
//             rho = m_phase_left->density();
//             //            doublereal rdz = 2.0/
//             //                 (m_flow_left->z(m_left_points-1) - 
//             //                     m_flow_left->z(m_left_points - 2));

//             for (k = 0; k < m_left_nsp; k++) 
//                 m_work[k + m_start_left] *= m_molwt_left[k];

//             int ileft = loc() - m_left_nv;

//             // if the energy equation is enabled at this point,
//             // set the gas temperature to the surface temperature
//             if (m_flow_left->doEnergy(pnt)) {
//                 rg[ileft + 2] = xg[ileft + 2] - m_sphase->temperature();
//             }

//             for (k = 1; k < m_left_nsp; k++) {
//                 if (enabled && m_flow_left->doSpecies(k)) {
//                     rg[ileft + 4 + k]  += m_work[k + m_start_left];
//                     //+= rdz*m_work[k + m_sp_left]/rho;
                    
//                 }
//             }
//         }

//         if (m_flow_right) {
//             for (k = 0; k < m_right_nsp; k++) 
//                 m_work[k + m_start_right] *= m_molwt_right[k];
            
//             int iright = loc() + m_nsp;
//             rg[iright + 2] -= m_sphase->temperature();
//             //r[iright + 3] = x[iright];
//             for (k = 0; k < m_right_nsp; k++) {
//                 rg[iright + 4 + k] -= m_work[k + m_start_right];
//             }
//         }
//     }

// void ChemSurf1D::
// save(XML_Node& o, doublereal* soln) {
//         doublereal* s = soln + loc();
//         XML_Node& srf = o.addChild("surface");
//         for (int k = 0; k < m_nsp; k++) {
//             ctml::addFloat(srf, componentName(k), s[k], "", "coverage",
//                 0.0, 1.0);
//         }
//     }


//     /////////////////////////////////////////////////////////////
//     //
//     // surf1D
//     //
//     ////////////////////////////////////////////////////////////

//     string ChemSurf1D::
//     componentName(int n) const {
//         /// @todo why dummy?
//         switch (n) {
//         case 0: return "temperature"; break;
//         default: return "<unknown>";
//         }
//     }

//     void ChemSurf1D::
//     init() {
//         if (m_index < 0) {
//             throw CanteraError("Surf1D", 
//                 "install in container before calling init.");
//         }
//         if (m_index > 0) m_left_flow = true;
//         resize(1,1);

//         // set bounds
//         vector_fp lower(1, 200.0);
//         vector_fp upper(1, 10000.0);
//         setBounds(1, lower.begin(), 1, upper.begin());

//         // set tolerances
//         vector_fp rtol(1, 1e-4);
//         vector_fp atol(1, 1.e-5);
//         setTolerances(1, rtol.begin(), 1, atol.begin());
//     }


//     void ChemSurf1D::eval(int jg, doublereal* xg, doublereal* rg, 
//         integer* diagg, doublereal rdt) {
//         int k;

//         // if computing a Jacobian (jg > 0), and the global point is
//         // outside the points the surface can influence, then skip
//         // evaluating the residual            
//         if (jg >= 0 && (jg < firstPoint() - 2 
//                 || jg > lastPoint() + 2)) return;

//         // start of local part of global arrays
//         doublereal* x = xg + loc();
//         doublereal* r = rg + loc();
//         integer* diag = diagg + loc();

//         // set the left gas state to the adjacent point

//         // gas-phase residuals
//         doublereal rho;
//     }

//     void ChemSurf1D::
//     save(XML_Node& o, doublereal* soln) {
//         doublereal* s = soln + loc();
//         XML_Node& srf = o.addChild("surface");
//         for (int k = 0; k < m_nsp; k++) {
//             ctml::addFloat(srf, componentName(k), s[k], "", "coverage",
//                 0.0, 1.0);
//         }
//     }





