/**
 * @file Solid1D.cpp
 */

/*
 * $Author: dggoodwin $
 * $Revision: 1.4 $
 * $Date: 2004/08/28 16:12:41 $
 */

// Copyright 2003  California Institute of Technology


// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include <stdlib.h>
#include <time.h>

#include "Solid1D.h"
#include "../ArrayViewer.h"
#include "../ctml.h"
#include "MultiJac.h"

using namespace ctml;

namespace Cantera {


    int Solid1D::c_T_loc = 0;
    int Solid1D::c_C_loc = 1;

    Solid1D::Solid1D(ThermoPhase* ph, int points) : 
        Domain1D(1, points),
        m_kin(0),
	m_trans(0),
        m_jac(0),
	m_ok(false)
    {
        m_type = cSolidType;
        m_points = points;
        m_thermo = ph;

        if (ph == 0) { m_nsp = 1; return; }// used to create a dummy object

        m_nsp = m_thermo->nSpecies();
        Domain1D::resize(m_nsp+1, points);

        // make a local copy of the species molecular weight vector
        m_wt = m_thermo->molecularWeights();

        m_nv = m_nsp + 1; 

        // turn off the energy equation at all points
        m_do_energy.resize(m_points,false);
        m_do_species.resize(m_nsp,false);

        m_diff.resize(m_nsp*m_points);
        m_flux.resize(m_nsp,m_points);
        m_wdot.resize(m_nsp,m_points, 0.0);
        m_cbar.resize(m_nsp);


        //-------------- default solution bounds --------------------

        vector_fp vmin(m_nv), vmax(m_nv);
        
        // temperature bounds
        vmin[c_T_loc] = 200.0;
        vmax[c_T_loc]= 1.e9;

        // concentration bounds
        int k;
        for (k = 0; k < m_nsp; k++) {
            vmin[c_C_loc + k] = -1.0e-5;
            vmax[c_C_loc + k] = 1.0e5;
        }
        setBounds(vmin.size(), vmin.begin(), vmax.size(), vmax.begin());


        //-------------------- default error tolerances ----------------
        vector_fp rtol(m_nv, 1.0e-8);
        vector_fp atol(m_nv, 1.0e-15);
        setTolerances(rtol.size(), rtol.begin(), atol.size(), atol.begin(),false);
        setTolerances(rtol.size(), rtol.begin(), atol.size(), atol.begin(),true);

        //-------------------- grid refinement -------------------------
        m_refiner->setActive(c_T_loc, false);

        vector_fp gr;
        for (int ng = 0; ng < m_points; ng++) gr.push_back(1.0*ng/m_points);
        setupGrid(m_points, gr.begin());
        setID("solid");
    }


    /**
     * Change the grid size. Called after grid refinement.
     */
    void Solid1D::resize(int points) {
        Domain1D::resize(m_nv, points);

        m_rho.resize(m_points, 0.0);
        m_wtm.resize(m_points, 0.0);
        m_cp.resize(m_points, 0.0);
        m_tcon.resize(m_points, 0.0);
        m_diff.resize(m_nsp*m_points);

        m_flux.resize(m_nsp,m_points);
        m_wdot.resize(m_nsp,m_points, 0.0);

        m_do_energy.resize(m_points,false);
        m_fixedtemp.resize(m_points);        

        m_dz.resize(m_points-1);
        m_z.resize(m_points);
    }        
        


    void Solid1D::setupGrid(int n, const doublereal* z) {        
        resize(n);
        int j;
        m_z[0] = z[0];
        for (j = 1; j < m_points; j++) {
            m_z[j] = z[j];
            m_dz[j-1] = m_z[j] - m_z[j-1];
        }
    }


    /**
     * Install a transport manager.
     */
    void Solid1D::setTransport(Transport& trans) {
        m_trans = &trans;

        if (m_trans->model() != cSolidTransport) {
            throw CanteraError("setTransport","unknown transport model.");
    }


    /**
     * Set the solid object state to be consistent with the solution at
     * point j.
     */
    void Solid1D::setThermoState(const doublereal* x,int j) {
        m_thermo->setTemperature(T(x,j));
        const doublereal* yy = x + m_nv*j + 1;
        m_thermo->setConcentrations(yy);
    }


    /**
     * Set the state to be consistent with the solution at the
     * midpoint between j and j + 1.
     */
    void Solid1D::setStateAtMidpoint(const doublereal* x,int j) {
        m_thermo->setTemperature(0.5*(T(x,j)+T(x,j+1)));
        const doublereal* ccj = x + m_nv*j + 1;
        const doublereal* ccjp = x + m_nv*(j+1) + 1;
        for (int k = 0; k < m_nsp; k++)
            m_ybar[k] = 0.5*(ccj[k] + ccjp[k]);
        m_thermo->setConcentrations(m_cbar.begin());
    }



    void Solid1D::eval(int jg, doublereal* xg, 
        doublereal* rg, integer* diagg, doublereal rdt) {

        // if evaluating a Jacobian, and the global point is outside
        // the domain of influence for this domain, then skip
        // evaluating the residual
        if (jg >=0 && (jg < firstPoint() - 1 || jg > lastPoint() + 1)) return;

        // if evaluating a Jacobian, compute the steady-state residual
        if (jg >= 0) rdt = 0.0;

        // start of local part of global arrays
        doublereal* x = xg + loc();
        doublereal* rsd = rg + loc();
        integer* diag = diagg + loc();
        
        int jmin, jmax, jpt;
        jpt = jg - firstPoint();

        if (jg < 0) {      // evaluate all points
            jmin = 0;
            jmax = m_points - 1;
        }
        else {            // evaluate points for Jacobian
            jmin = max(jpt-1, 0);
            jmax = min(jpt+1,m_points-1);
        }

        // properties are computed for grid points from j0 to j1
        int j0 = max(jmin-1,0);
        int j1 = min(jmax+1,m_points-1);


        int j, k;

        //----------------------------------------------------- 
        //              update properties
        //-----------------------------------------------------

        // thermodynamic properties only if a Jacobian is
        // not being evaluated
        if (jpt < 0) 
            updateThermo(x, j0, j1);

        // update transport properties only if a Jacobian is
        // not being evaluated
        if (jpt < 0) 
            updateTransport(x, j0, j1);

        // update the species diffusive mass fluxes whether or not a
        // Jacobian is being evaluated
        updateDiffFluxes(x, j0, j1);

        for (j = j0; j <= j1; j++) {
            setThermoState(j);
        }

        //----------------------------------------------------
        // evaluate the residual equations at all required
        // grid points
        //----------------------------------------------------

        for (j = jmin; j <= jmax; j++) {


            //----------------------------------------------
            //         left boundary
            //----------------------------------------------

            if (j == 0) {
                rsd[index(c_T_loc,0)] = T(x,0);

                 // The default boundary condition for species is zero
                 // flux. However, the boundary object may modify
                 // this.
                 for (k = 0; k < m_nsp; k++) {
                     rsd[index(c_C_loc + k, 0)] = - m_flux(k,0);
                 }
            }


            //----------------------------------------------
            //
            //         right boundary
            //
            //----------------------------------------------

            else if (j == m_points - 1) {
                rsd[index(c_T_loc,j)] = T(x,j);
                for (k = 0; k < m_nsp; k++) {
                    rsd[index(k+c_C_loc,j)] = m_flux(k,j-1);
                }
            }


            //------------------------------------------
            //     interior points                 
            //------------------------------------------
            
            else {

                //-------------------------------------------------
                //    Species equations
                //                
                //   \rho u dY_k/dz + dJ_k/dz + M_k\omega_k
                //
                //-------------------------------------------------
                getWdot(x,j);

                doublereal diffus;
                for (k = 0; k < m_nsp; k++) {
                    diffus = 2.0*(m_flux(k,j) - m_flux(k,j-1))
                             /(z(j+1) - z(j-1));
                    rsd[index(c_C_loc + k, j)]   
                        = wdot(k,j) - diffus
                        - rdt*(C(x,k,j) - C_prev(k,j));
                    diag[index(c_C_loc + k, j)] = 1;
                }

                //-----------------------------------------------
                //    energy equation
                //-----------------------------------------------

                if (m_do_energy[j]) {

                    rsd[index(c_T_loc, j)]   = - divHeatFlux(x,j);
                    rsd[index(c_T_loc, j)] /= (m_rho[j]*m_cp[j]);

                    rsd[index(c_T_loc, j)] -= rdt*(T(x,j) - T_prev(j));
                    diag[index(c_T_loc, j)] = 1;
                }
            }

            // residual equations if the energy or species equations
            // are disabled

            if (!m_do_energy[j]) {
                rsd[index(c_T_loc, j)] = T(x,j) - T_fixed(j);
                diag[index(c_T_loc, j)] = 0;
            }
        }
    }



    /**
     * Update the transport properties at grid points in the range
     * from j0 to j1, based on solution x.
     */
    void Surf1D::updateTransport(doublereal* x,int j0, int j1) {
        int j;
        for (j = j0; j < j1; j++) {
            setStateAtMidpoint(x,j);
            m_trans->getMixDiffCoeffs(m_diff.begin() + j*m_nsp);
            m_tcon[j] = m_trans->thermalConductivity();
        }
    }


    /**
     * Print the solution.
     */    
    void Solid1D::showSolution(const doublereal* x) {
        int nn = m_nv/5;
        int i, j, n;
        char* buf = new char[100];

        // The mean molecular weight is needed to convert
        updateThermo(x, 0, m_points-1);

        for (i = 0; i < nn; i++) {
            drawline();
            sprintf(buf, "\n        z   ");
            writelog(buf);
            for (n = 0; n < 5; n++) { 
                sprintf(buf, " %10s ",componentName(i*5 + n).c_str());
                writelog(buf);
            }
            drawline();
            for (j = 0; j < m_points; j++) {
                sprintf(buf, "\n %10.4g ",m_z[j]);
                writelog(buf);
                for (n = 0; n < 5; n++) { 
                    sprintf(buf, " %10.4g ",component(x, i*5+n,j));
                    writelog(buf);
                }
            }
            writelog("\n");
        }
        int nrem = m_nv - 5*nn;
        drawline();
        sprintf(buf, "\n        z   ");
        writelog(buf);
        for (n = 0; n < nrem; n++) {
            sprintf(buf, " %10s ", componentName(nn*5 + n).c_str());
            writelog(buf);
        }
        drawline();
        for (j = 0; j < m_points; j++) {
            sprintf(buf, "\n %10.4g ",m_z[j]);
            writelog(buf);
            for (n = 0; n < nrem; n++) { 
                sprintf(buf, " %10.4g ",component(x, nn*5+n,j));
                writelog(buf);
            }
        }
        writelog("\n");
    }


    /**
     * Update the diffusive mass fluxes.
     */
    void Solid1D::updateDiffFluxes(const doublereal* x, int j0, int j1) {
        int j, k, m;
        doublereal sum, wtm, rho, dz, gradlogT, s;
        doublereal dphidz, a1;
        for (j = j0; j < j1; j++) {
            sum = 0.0;
            rho = density(j);
            dz = z(j+1) - z(j);
            for (k = 0; k < m_nsp; k++) {
                m_flux(k,j) = m_diff[k+m_nsp*j] * 
                              (C(x,k,j) - C(x,k,j+1))/dz;
                sum -= m_flux(k,j);
            }
            for (k = 0; k < m_nsp; k++) m_flux(k,j) += C(x,k,j)*sum;
        } 
        break;
    }


    void Solid1D::outputTEC(ostream &s, const doublereal* x, 
        string title, int zone) {
        int j,k;
        s << "TITLE     = \"" + title + "\"" << endl;
        s << "VARIABLES = \"Z (m)\"" << endl;
        s << "\"T (K)\"" << endl;

        for (k = 0; k < m_nsp; k++) {
            s << "\"" << m_thermo->speciesName(k) << "\"" << endl;
        }
        s << "ZONE T=\"c" << zone << "\"" << endl;
        s << " I=" << m_points << ",J=1,K=1,F=POINT" << endl;
        s << "DT=(SINGLE";
        for (k = 0; k < m_nsp; k++) s << " SINGLE";
        s << " )" << endl;
        for (j = 0; j < m_points; j++) {
            s << z(j) << " ";
            for (k = 0; k < m_nv; k++) {
                s << component(x, k, j) << " ";
            }
            s << endl;
        }
    }


    string Solid1D::componentName(int n) const {
        switch(n) {
        case c_T_loc: return "T";
        default:
            if (n >= (int) 1 && n < (int) (c_C_loc + m_nsp)) {
                    return m_thermo->speciesName(n - 1);
            }
            else 
                return "<unknown>";
        }
    }


    void Solid1D::restore(XML_Node& dom, doublereal* soln) {

        vector<string> ignored;
        int nsp = m_thermo->nSpecies();
        vector_int did_species(nsp, 0);

        vector<XML_Node*> str;
        dom.getChildren("string",str);
        int nstr = str.size();
        for (int istr = 0; istr < nstr; istr++) {
            XML_Node& nd = *str[istr];
            writelog(nd["title"]+": "+nd.value()+"\n");
        }

        map<string, double> params;
        getFloats(dom, params);

        vector<XML_Node*> d;
        dom.child("grid_data").getChildren("floatArray",d);
        int nd = d.size();

        vector_fp x;
        int n, np, j, ks, k;
        string nm;
        bool readgrid = false, wrote_header = false;
        for (n = 0; n < nd; n++) {
            XML_Node& fa = *d[n];
            nm = fa["title"];
            if (nm == "z") {
                getFloatArray(fa,x,false);
                np = x.size();
                writelog("Grid contains "+int2str(np)+
                    " points.\n");
                readgrid = true;

                // note that setupGrid also resizes the domain.
                setupGrid(np, x.begin());
            }
        }
        if (!readgrid) {
            throw CanteraError("Solid1D::restore",
                "domain contains no grid points.");
        }

        writelog("Importing datasets:\n");
        for (n = 0; n < nd; n++) {
            XML_Node& fa = *d[n];
            nm = fa["title"];
            getFloatArray(fa,x,false);
            if (nm == "z") {
                ;   // already read grid
            }
            else if (nm == "T") {
                writelog("temperature   ");
                if ((int) x.size() == np) {
                    for (j = 0; j < np; j++)
                        soln[index(c_T_loc,j)] = x[j];

                    // For fixed-temperature simulations, use the imported temperature profile by default. 
                    // If this is not desired, call setFixedTempProfile *after* restoring the solution.
                    vector_fp zz(np);
                    for (int jj = 0; jj < np; jj++) zz[jj] = (grid(jj) - zmin())/(zmax() - zmin());
                    setFixedTempProfile(zz, x);
                }
                else goto error;
            }
            else if (m_thermo->speciesIndex(nm) >= 0) {
                writelog(nm+"   ");
                if ((int) x.size() == np) {
                    k = m_thermo->speciesIndex(nm);
                    did_species[k] = 1;
                    for (j = 0; j < np; j++) 
                        soln[index(k+c_C_loc,j)] = x[j];
                }
            }
            else
                ignored.push_back(nm);
        }

        if (ignored.size() != 0) {
            writelog("\n\n");
            writelog("Ignoring datasets:\n");
            int nn = ignored.size();
            for (int n = 0; n < nn; n++) {
                writelog(ignored[n]+"   ");
            }
        }

        for (ks = 0; ks < nsp; ks++) {
            if (did_species[ks] == 0) {
                if (!wrote_header) {
                    writelog("Missing data for species:\n");
                    wrote_header = true;
                }
                writelog(m_thermo->speciesName(ks)+" ");
            }
        }

        return;
 error:
        throw CanteraError("Solid1D::restore","Data size error");
    }



    void Solid1D::save(XML_Node& o, doublereal* sol) {
        int k;

        ArrayViewer soln(m_nv, m_points, sol + loc());

        XML_Node& flow = (XML_Node&)o.addChild("domain");
        flow.addAttribute("type",flowType());
        flow.addAttribute("id",m_id);
        flow.addAttribute("points",m_points);
        flow.addAttribute("components",m_nv);

        if (m_desc != "") addString(flow,"description",m_desc);
        XML_Node& gv = flow.addChild("grid_data");
        addFloatArray(gv,"z",m_z.size(),m_z.begin(),
            "m","length");
        vector_fp x(soln.nColumns());

        soln.getRow(c_T_loc,x.begin());
        addFloatArray(gv,"T",x.size(),x.begin(),"K","temperature",0.0);

        for (k = 0; k < m_nsp; k++) {
            soln.getRow(c_C_loc+k,x.begin());
            addFloatArray(gv,m_thermo->speciesName(k),
                x.size(),x.begin(),"","concentration",0.0,1.0);
        }
    }


    void Solid1D::setJac(MultiJac* jac) {
        m_jac = jac;
    }


}
