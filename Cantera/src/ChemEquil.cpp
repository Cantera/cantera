/**
 *
 *  @file ChemEquil.cpp
 *
 *  Chemical equilibrium. 
 *  Implementation file for class ChemEquil
 *
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2001 California Institute of Technology
 *
 */

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include <vector>
using namespace std;

#include "ChemEquil.h"
#include "DenseMatrix.h"
#include "recipes.h"
#include "sort.h"
#include "PropertyCalculator.h"
#include "ctexceptions.h"
#include "vec_functions.h"
#include "stringUtils.h"


namespace Cantera {

    int _equilflag(const char* xy) {
        string flag = string(xy);
        if (flag == "TP") return TP;
        else if (flag == "TV") return TV;
        else if (flag == "HP") return HP;
        else if (flag == "UV") return UV;
        else if (flag == "SP") return SP;
        else if (flag == "SV") return SV;
        else if (flag == "UP") return UP;
        else throw CanteraError("_equilflag","unknown property pair "+flag);
        }


    //-----------------------------------------------------------
    //                 construction / destruction 
    //-----------------------------------------------------------


    ///  Default Constructor.
    ChemEquil::ChemEquil() : m_skip(-1), m_p1(0), m_p2(0), m_p0(OneAtm), m_eloc(-1),
                             m_abscharge(Tiny)
    {}


    /// Destructor
    ChemEquil::~ChemEquil(){
        delete m_p1;
        delete m_p2;
    }


    /**
     *  Prepare for equilibrium calculations with a specified
     *  mixture.
     *  @param s  mixture
     */
    void ChemEquil::initialize(thermo_t& s) 
    {
        // store a pointer to s and some of its properties locally
        m_thermo = &s;
        m_phase = &s;

        m_p0 = s.refPressure(); 
        m_kk = m_phase->nSpecies();
        m_mm = m_phase->nElements();
        if (m_kk < m_mm) {
            throw CanteraError("ChemEquil::initialize",
                "number of species cannot be less than the number of elements.");
        }
        
        // allocate space in internal work arrays
        m_molefractions.resize(m_kk);
        m_lambda.resize(m_mm, -100.0);
        m_elementmolefracs.resize(m_mm);
        m_comp.resize(m_mm * m_kk);
        m_jwork1.resize(m_mm+2);
        m_jwork2.resize(m_mm+2);
        m_startSoln.resize(m_mm+1);
        m_grt.resize(m_kk);
        m_mu_RT.resize(m_kk);

        // set up elemental composition matrix
        int m, k, mneg = -1;
        doublereal na, ewt;
        for (m = 0; m < m_mm; m++) {
            for (k = 0; k < m_kk; k++) {
                na = m_phase->nAtoms(k,m);

                // handle the case of negative atom numbers (used to 
                // represent positive ions)
                if (na < 0.0) {

                    // if negative atom numbers have already been specified
                    // for some element other than this one, throw
                    // an exception
                    if (mneg >= 0 && mneg != m) 
                        throw CanteraError("ChemEquil::initialize",
                            "negative atom numbers allowed for only one element"); 
                    mneg = m;
                    ewt = m_phase->atomicWeight(m);

                    // the element should be an electron... if it isn't 
                    // print a warning.
                    if (ewt > 1.0e-3) 
                        writelog(string("WARNING: species "
                                +m_phase->speciesName(k)
                                +" has "+fp2str(m_phase->nAtoms(k,m))
                                +" atoms of element "
                                +m_phase->elementName(m)+
                                ", but this element is not an electron.\n"));
                }
            }
        }
        m_eloc = mneg;

        // set up the elemental composition matrix
        for (k = 0; k < m_kk; k++) {
            for (m = 0; m < m_mm; m++) {
                m_comp[k*m_mm + m] = m_phase->nAtoms(k,m);
            }
        }
    }


    /** 
     * Set mixture to an equilibrium state consistent with specified 
     * element potentials and temperature.
     *
     * @param lambda_RT vector of non-dimensional element potentials
     * \f[ \lambda_m/RT \f].
     * @param t temperature in K.
     *
     */
    void ChemEquil::setToEquilState(thermo_t& s, 
        const vector_fp& lambda_RT, doublereal t) 
    {
        // compute the chemical potentials by summing element potentials
        fill(m_mu_RT.begin(), m_mu_RT.end(), 0.0);
        for (int k = 0; k < m_kk; k++)
            for (int m = 0; m < m_mm; m++) 
                m_mu_RT[k] += lambda_RT[m]*nAtoms(k,m);

        // set the temperature
        s.setTemperature(t);

        // call the phase-specific method to set the phase to the
        // equilibrium state with the specified species chemical
        // potentials.
        s.setToEquilState(m_mu_RT.begin());
        update(s);
    }


    /** 
     *  update internally stored state information.
     * @todo argument not used.
     */
    void ChemEquil::update(const thermo_t& s) {

        // get the mole fractions, temperature, and density
        m_phase->getMoleFractions(m_molefractions.begin());
        m_temp = m_phase->temperature();
        m_dens = m_phase->density();

        // compute the elemental mole fractions
        doublereal sum = 0.0;
        int m, k;
        for (m = 0; m < m_mm; m++) {
            m_elementmolefracs[m] = 0.0;
            for (k = 0; k < m_kk; k++) {
                m_elementmolefracs[m] += nAtoms(k,m) * m_molefractions[k];
                if (m_molefractions[k] < 0.0) {
                    throw CanteraError("update",
                        "negative mole fraction for "+m_phase->speciesName(k)+
                        ": "+fp2str(m_molefractions[k]));
                }
            }
            sum += m_elementmolefracs[m];
        }

        // normalize the element mole fractions
        for (m = 0; m < m_mm; m++) m_elementmolefracs[m] /= sum;
    }

        
    /**
     *  Estimate the initial mole fractions.  Uses the Simplex method
     *  to estimate the initial number of moles of each species.  The
     *  linear Gibbs minimization problem is solved, neglecting the
     *  free energy of mixing terms. This procedure produces a good
     *  estimate of the low-temperature equilibrium composition.
     *
     *  @param s             phase object
     *  @param elementMoles  vector of elemental moles
     */
    int ChemEquil::setInitialMoles(thermo_t& s, 
        vector_fp& elementMoles) 
    {
        int m, n;
        double pres = s.pressure();
        double lp = log(pres/m_p0);
        integer mm = m_phase->nElements();
        integer kksp = m_phase->nSpecies();
        
        DenseMatrix aa(mm+2, kksp+1, 0.0);
        
        // first column contains fixed element moles
        for (m = 0; m < mm; m++) {
            aa(m+1,0) = elementMoles[m] + 0.01;
        }

        // get the array of non-dimensional Gibbs functions for the pure 
        // species
        s.getGibbs_RT(m_grt.begin());
        
        int kpp = 0;
        for (int k = 0; k < kksp; k++) {
            kpp++;
            aa(0, kpp) =  -m_grt[k];
            aa(0, kpp) -= lp;               // ideal gas
            for (int q = 0; q < mm; q++)                           
                aa(q+1, kpp) = -nAtoms(k, q);
        }

        integer mp = mm+2;               // parameters for SIMPLX
        integer np = kksp+1;
        integer m1 = 0;
        integer m2 = 0;
        integer m3 = mm;
        integer icase=0;

        vector_int iposv(mm);
        vector_int izrov(kksp);
    
        //  solve the linear programming problem

        simplx_(&aa(0,0), &mm, &kksp, &mp, &np, &m1, &m2, &m3, 
            &icase, izrov.begin(), iposv.begin());
        
        fill(m_molefractions.begin(), m_molefractions.end(), 0.0);
        for (n = 0; n < mm; n++) {
            int ksp = 0;
            int ip = iposv[n] - 1;
            for (int k = 0; k < kksp; k++) { 
                if (ip == ksp) {
                    m_molefractions[k] = aa(n+1, 0);
                }
                ksp++;
            }
        }
        s.setState_PX(pres, m_molefractions.begin());
        update(s);
        return icase;
    }

 
    /**
     * Generate a starting estimate for the element potentials.
     */
    int ChemEquil::estimateElementPotentials(thermo_t& s, vector_fp& lambda) 
    {
        int k, ksp, m, n;
        for (k = 0; k < m_kk; k++) {
            if (m_molefractions[k] > 0.0) {
                m_molefractions[k] = fmaxx(m_molefractions[k], 0.05);
            }
            //else
            //     m_molefractions[k] = 0.001;
        }
        s.setState_PX(s.pressure(), m_molefractions.begin());

        // sort mole fractions
        vector_fp mol(m_kk, 0.0);
        vector_int index(m_kk, 0);
        for (k = 0; k < m_kk; k++) {
            mol[k] = m_molefractions[k];
            index[k] = k;
        }
        heapsort(mol, index);

        DenseMatrix aa(m_mm, m_mm, 0.0);
        vector_fp b(m_mm, -999.0);
        vector_fp ipvt(m_mm, 0);
        
        // find a set of constituents
        vector_int kc(m_mm, 0);
        vector_fp tmp(m_mm, 0.0);
        vector_fp mu_RT(m_kk, 0.0);

        s.getChemPotentials(mu_RT.begin());
        doublereal rrt = 1.0/(GasConstant*m_phase->temperature());
        scale(mu_RT.begin(), mu_RT.end(), mu_RT.begin(), rrt);
        int j = 0;
        for (k = m_kk - 1; k >= 0; k--) {
            ksp = index[k];
            if ( mol[k] > 0.0 ) {
                kc[j] = ksp;
                j++;
                if (j == m_mm) break;
            }
        }
        //if (j < m_mm) 
        //    return -1;
        //throw CanteraError("estimateElementPotentials",
        //      "too few species (" + int2str(j) + ").");

        for (m = 0; m < j; m++) {

            for (n = 0; n < m_mm; n++) {
                aa(m,n) = nAtoms(kc[m], n);
            }
            b[m] = mu_RT[kc[m]];
        }
        for (m = j+1; m < m_mm; m++) {
            aa(m,m) = 1.0;
            b[m] = lambda[m];
        }

        int info;
        try {
            info = solve(aa, b.begin());
        }
        catch (CanteraError) {
            return -2; 
        }

        if (info == 0) {
            for (m = 0; m < m_mm; m++) {
                lambda[m] = b[m];
            }
        }
        return info;
    }


    /**
     * Equilibrate a phase, holding the elemental composition fixed
     * at the initial vaollue.
     */
    int ChemEquil::equilibrate(thermo_t& s, int XY) {
        vector_fp emol(s.nElements());
        initialize(s);
        update(s);
        copy(m_elementmolefracs.begin(), m_elementmolefracs.end(),
            emol.begin());
        return equilibrate(s, XY, emol);
    }
    

    /**
     *   compute the equilibrium composition for 2 specified
     *   properties and specified element moles.
     */
    int ChemEquil::equilibrate(thermo_t& s, int XY, vector_fp& elMoles)
    {
        doublereal xval, yval;
        int fail = 0;

        delete m_p1;
        delete m_p2;
        bool tempFixed = true;
        initialize(s);
        update(s);
        switch (XY) {
        case TP: case PT:
            m_p1 = new TemperatureCalculator<thermo_t>;
            m_p2 = new PressureCalculator<thermo_t>; break;
        case HP: case PH:
            tempFixed = false;
            m_p1 = new EnthalpyCalculator<thermo_t>;
            m_p2 = new PressureCalculator<thermo_t>; break;
        case SP: case PS:
            tempFixed = false;
            m_p1 = new EntropyCalculator<thermo_t>;
            m_p2 = new PressureCalculator<thermo_t>; break;
        case SV: case VS:
            tempFixed = false;
            m_p1 = new EntropyCalculator<thermo_t>;
            m_p2 = new DensityCalculator<thermo_t>; break;
        case TV: case VT:
            m_p1 = new TemperatureCalculator<thermo_t>;
            m_p2 = new DensityCalculator<thermo_t>; break;
        case UV: case VU:
            tempFixed = false;
            m_p1 = new IntEnergyCalculator<thermo_t>;
            m_p2 = new DensityCalculator<thermo_t>; break;
        default:
            throw CanteraError("equilibrate","illegal property pair.");
        }

        // If the temperature is one of the specified variables, and
        // it is outside the valid range, throw an exception.
        if (tempFixed) {
            double tfixed = s.temperature();
            if (tfixed > s.maxTemp() + 1.0 || tfixed < s.minTemp() - 1.0) {
                throw CanteraError("ChemEquil","Specified temperature ("
                    +fp2str(m_thermo->temperature())+" K) outside "
                    "valid range of "+fp2str(m_thermo->minTemp())+" K to "
                    +fp2str(m_thermo->maxTemp())+" K\n");
            }                
        }
        xval = m_p1->value(s);
        yval = m_p2->value(s);
        int mm = m_mm;
        
        int m;
        int nvar = mm + 1;
        
        DenseMatrix jac(nvar, nvar);       // jacobian
        vector_fp x(nvar, -102.0);          // solution vector
    
        vector_fp res_trial(nvar);
        vector_fp elementMol(mm, 0.0);
        double perturb;
        for (m = 0; m < mm; m++) {
            if (m_skip < 0 && elMoles[m] > 0.0 ) m_skip = m;
#define PERTURB_ELEMENT_MOLES
#ifdef PERTURB_ELEMENT_MOLES
            perturb = Cutoff*(1.0 + rand());
#else
            perturb = 0.0;
#endif
            elementMol[m] = elMoles[m] + perturb;
        }

        update(s);


        // loop to estimate T
        if (!tempFixed) {
            
            doublereal tmax = m_thermo->maxTemp();
            doublereal tmin = m_thermo->minTemp();
            doublereal slope, phigh, plow, pval, dt;


            // first get the property values at the upper and lower
            // temperature limits. Since p1 (h, s, or u) is monotonic
            // in T, these values determine the upper and lower
            // bounnds (phigh, plow) for p1.

            m_phase->setTemperature(tmax);
            setInitialMoles(s, elementMol);
            phigh = m_p1->value(s);

            m_phase->setTemperature(tmin);
            setInitialMoles(s, elementMol);
            plow = m_p1->value(s);

            // start with T at the midpoint of the range
            doublereal t0 = 0.5*(tmin + tmax);
            m_phase->setTemperature(t0);

            // loop up to 5 times
            for (int it = 0; it < 5; it++) {

                // set the composition and get p1
                setInitialMoles(s, elementMol);
                pval = m_p1->value(s);


                // If this value of p1 is greater than the specified
                // property value, then the current temperature is too
                // high. Use it as the new upper bound. Otherwise, it
                // is too low, so use it as the new lower bound.
                if (pval > xval) { 
                    tmax = t0; 
                    phigh = pval; 
                }
                else { 
                    tmin = t0; 
                    plow = pval; 
                }

                // Determine the new T estimate by linearly intepolation
                // between the upper and lower bounds
                slope = (phigh - plow)/(tmax - tmin);
                dt = (xval - plow)/slope;

                // If within 100 K, terminate the search
                if (fabs(dt) < 100.0) break;

                // update the T estimate
                t0 = tmin + dt;
                m_phase->setTemperature(t0);
            }
        }

        if (m_lambda[0] == -100.0) {
            setInitialMoles(s, elementMol);
            for (int ii = 0; ii < m_mm; ii++) x[ii] = -101.0;
            estimateElementPotentials(s, x);
        }
        else {
            doublereal rt = GasConstant * m_phase->temperature();
            for (int ii = 0; ii < m_mm; ii++) x[ii] = m_lambda[ii]/rt;
        }

        x[m_mm] = log(m_phase->temperature());

        vector_fp above(nvar);
        vector_fp below(nvar);

        for (m = 0; m < mm; m++) {
            above[m] = 30.0;
            below[m] = -2000.0;
            if (elMoles[m] < Cutoff && m != m_eloc) x[m] = -1000.0;
            //if (m == m_eloc) x[m] = -10.0;
        }
        above[mm] = log(m_thermo->maxTemp() + 1.0);  //log(1.e4);
        below[mm] = log(m_thermo->minTemp() - 1.0); //log(10.0);

        vector_fp grad(nvar, 0.0);        // gradient of f = F*F/2
        vector_fp oldx(nvar, 0.0);        // old solution
        vector_fp prevx(nvar, 0.0);       // old solution
        vector_fp oldresid(nvar, 0.0);
        doublereal f, oldf;

        int iter = 0;
        int info=0;
        doublereal fctr = 1.0, newval;

        goto converge;

 next:


        // if the problem involves charged species, then the
        // "electron" element equation is a charge balance. Compute
        // the sum of the absolute values of the charge to use as the
        // normalizing factor.
        if (m_eloc >= 0) {
            m_abscharge = 0.0;
            int k;
            for (k = 0; k < m_kk; k++) 
                m_abscharge += fabs(m_phase->charge(k)*m_molefractions[k]);
        }


        iter++;

        // compute the residual and the jacobian using the current
        // solution vector
        equilResidual(s, x, elMoles, res_trial, XY, xval, yval);
        f = 0.5*dot(res_trial.begin(), res_trial.end(), res_trial.begin());
        equilJacobian(s, x, elMoles, jac, XY, xval, yval);
    
        // compute grad f = F*J
        jac.leftMult(res_trial.begin(), grad.begin());
        copy(x.begin(), x.end(), oldx.begin());
        copy(oldx.begin(), oldx.end(), prevx.begin());
        oldf = f;
        scale(res_trial.begin(), res_trial.end(), res_trial.begin(), -1.0);
        try {
            info = solve(jac, res_trial.begin());
        }
        catch (CanteraError) {
            cout << x << endl;
            //cout << res_trial << endl;
            //cout << grad << endl;
            cout << elMoles << endl;
            //cout << jac << endl;

            //cout << "m_skip = " << m_skip << endl;
            throw CanteraError("equilibrate",
                "Jacobian is singular. \nTry adding more species, "
                "changing the elemental composition slightly, \nor removing "
                "unused elements.");
            return -3;
        }

        // find the factor by which the Newton step can be multiplied
        // to keep the solution within bounds.
        fctr = 1.0;
        for (m = 0; m < nvar; m++) {
            newval = x[m] + res_trial[m];
            if (newval > above[m]) {
                fctr = fmaxx( 0.0, fminn( fctr, 
                    0.4*(above[m] - x[m])/(newval - x[m])));
            }
            else if (newval < below[m]) {
                fctr = fminn(fctr, 0.4*(x[m] - below[m])/(x[m] - newval));
            }
        }

        // multiply the step by the scaing factor
        scale(res_trial.begin(), res_trial.end(), res_trial.begin(), fctr);

        if (!dampStep(s, oldx, oldf, grad, res_trial, 
            x, f, elMoles , XY, xval, yval))
        {
            fail++;
            if (fail > 3) {
                throw CanteraError("equilibrate",
                    "Cannot find an acceptable Newton damping coefficient.");
                return -4;
            }
        }
        else fail = 0;

converge:    

        //  check for convergence.
    
        equilResidual(s, x, elMoles, res_trial, XY, xval, yval);
        f = 0.5*dot(res_trial.begin(), res_trial.end(), res_trial.begin());
        doublereal xx, yy, deltax, deltay;
        xx = m_p1->value(s);
        yy = m_p2->value(s);
        deltax = (xx - xval)/xval;
        deltay = (yy - yval)/yval;

        if (absmax(res_trial.begin(), res_trial.end()) < options.relTolerance 
            && fabs(deltax) < options.relTolerance 
            && fabs(deltay) < options.relTolerance) {
            options.iterations = iter;

            m_lambda.resize(m_mm);
            doublereal rt = GasConstant*m_thermo->temperature();
            for (m = 0; m < m_mm; m++) {
                m_lambda[m] = x[m]*rt;
            }

            if (m_thermo->temperature() > m_thermo->maxTemp() + 1.0 ||
                m_thermo->temperature() < m_thermo->minTemp() - 1.0 ) {
                writelog("Warning: Temperature ("
                    +fp2str(m_thermo->temperature())+" K) outside "
                    "valid range of "+fp2str(m_thermo->minTemp())+" K to "
                    +fp2str(m_thermo->maxTemp())+" K\n");
            }
            
            return 0;
        }
    
        // no convergence
    
        if (iter > options.maxIterations) {
            throw CanteraError("equilibrate",
                "no convergence in "+int2str(options.maxIterations)
                +" iterations.");
            return -1;
        }
        goto next;
    }   


    int ChemEquil::dampStep(thermo_t& mix, vector_fp& oldx, 
        double oldf, vector_fp& grad, vector_fp& step, vector_fp& x, 
        double& f, vector_fp& elmols, int XY, double xval, double yval )
    {
        int nvar = x.size();
    
        double slope;
        double f2 = 0.0;
        double oldf2 = 0.0;
        double alpha = 1.e-4;
        double tmpdamp = 0.0; 
        double rhs1;
        double rhs2;
        double damp = 1.0;
        double damp2=0.0;
        double a;
        double bb;
        double disc;
        double minDamp = 0.0;
        double xTol = 1.e-7;

        static vector_fp res_new(nvar); // fix

        //slope = grad * step;
        slope = dot(grad.begin(), grad.end(), step.begin());
        double temp, test = 0.0;
  
        for (int i=0; i<nvar; i++)
        {
            temp = fabs(step[i]/fmaxx(fabs(oldx[i]), 1.0));
            if (temp > test) test = temp;
        }
        minDamp = xTol/test;
    
 retry:

        x = step;
        scale(x, damp);
        add_each(x, oldx);
    

        equilResidual(mix, x, elmols, res_new, XY, xval, yval);
        //f = 0.5*(res_new*res_new);
        f = 0.5*dot(res_new.begin(), res_new.end(), res_new.begin());
        if (damp < minDamp && damp < 1.0)
        {
            return 0;          // check that this is not a spurious min of f
        }
        else if (f <= oldf + alpha * damp * slope)
        {
            return 1;          // good damping coefficient
        }
        else
        {
            if (damp == 1.0)   // first time
            {
                tmpdamp = -slope/(2.0*(f - oldf - slope));
            }
            else
            {
                rhs1 = f - oldf - damp*slope;
                rhs2 = f2 - oldf2 - damp2*slope;
                a = (rhs1/(damp*damp) - rhs2/(damp2*damp2))/(damp - damp2);
                bb = (-damp2*rhs1/(damp*damp) + damp*rhs2/(damp2*damp2))
                     /(damp - damp2);
	  
                if (a == 0.0)
                    tmpdamp = -slope/(2.0*bb);
                else
                {
                    disc = bb*bb - 3.0*a*slope;
                    if (disc < 0.0)
                        tmpdamp = -slope/(2.0*bb);
                    else
                        tmpdamp = (-bb +sqrt(disc))/(3.0*a);
                }
                if (tmpdamp > 0.5*damp) tmpdamp = 0.5*damp;
            }
      
            damp2 = damp;
            f2 = f;
            oldf2 = oldf;
            damp = fmaxx(tmpdamp, 0.1*damp);
            goto retry;
        }
    }


    /**
     *  evaluates the residual vector F, of length mm 
     */
    void ChemEquil::equilResidual(thermo_t& mix, const vector_fp& x, 
        const vector_fp& elmtotal, vector_fp& resid, 
        int XY, doublereal xval, doublereal yval)
    {
        int n;
        doublereal xx, yy;
        doublereal temp = exp(x[m_mm]);
        setToEquilState(mix, x, temp);

        // residuals are the total element moles
        vector_fp& elm = m_elementmolefracs; 
        for (n=0; n < m_mm; n++)
        {
            // drive element potential for absent elements to -1000
            if (elmtotal[n] < Cutoff && n != m_eloc)
                resid[n] = x[n] + 1000.0;
            else
                //                resid[n] = elmtotal[n] - elm[n]; // log( (1.0 + elmtotal[n]) / (1.0 + elm[n]) );
                resid[n] = log( (1.0 + elmtotal[n]) / (1.0 + elm[n]) );
        }
        if (m_eloc >= 0) {
            doublereal chrg, sumnet = 0.0, sumabs = 0.0;
            for (int k = 0; k < m_kk; k++) {
                chrg = m_molefractions[k]*m_phase->charge(k);
                sumnet += chrg;
                sumabs += fabs(chrg);
            }
            resid[m_eloc] = sumnet/m_abscharge; // log((1.0 + sumnet/sumabs));
        }
        xx = m_p1->value(mix);
        yy = m_p2->value(mix);
        resid[m_mm] = xx/xval - 1.0; 
        resid[m_skip] = yy/yval - 1.0; 
    }


    //-------------------- Jacobian evaluation ---------------------------

    void ChemEquil::equilJacobian(thermo_t& mix, vector_fp& x,  
        const vector_fp& elmols, DenseMatrix& jac, 
        int XY, doublereal xval, doublereal yval)
    {
        int len = x.size();
        vector_fp& r0 = m_jwork1;
        vector_fp& r1 = m_jwork2;
        r0.resize(len);
        r1.resize(len);

        int n, m;
        doublereal rdx, dx, xsave;
        doublereal atol = 1.e-10;
    
        equilResidual(mix, x, elmols, r0, XY, xval, yval);
    
        for (n = 0; n < len; n++)
        {
            // perturb x(n)
        
            xsave = x[n];
            dx = atol;
            x[n] = xsave + dx;
            dx = x[n] - xsave;
            rdx = 1.0/dx;

            // calculate perturbed residual
        
            equilResidual(mix, x, elmols, r1, XY, xval, yval);
        
            // compute nth column of Jacobian
        
            for (m = 0; m < len; m++) {
                jac(m, n) = (r1[m] - r0[m])*rdx;
            }        
            x[n] = xsave;
        }
    }

} // namespace


// $Log: ChemEquil.cpp,v
