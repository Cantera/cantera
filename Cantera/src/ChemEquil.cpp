/**
 *
 *  @file ChemEquil.cpp
 *
 *  Chemical equilibrium.  Implementation file for class
 *  ChemEquil. 
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

#include "sort.h"
#include "PropertyCalculator.h"
#include "ctexceptions.h"
#include "vec_functions.h"
#include "stringUtils.h"
#include "MultiPhase.h"

namespace Cantera {

    /// map property strings to integers
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
     *  Prepare for equilibrium calculations.
     *  @param s object representing the solution phase.
     */
    void ChemEquil::initialize(thermo_t& s) 
    {
        // store a pointer to s and some of its properties locally.
        // Note: the use of two pointers is a historical artifact.
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
        m_component.resize(m_mm,-2);

        // set up elemental composition matrix
        int m, k, mneg = -1;
        doublereal na, ewt;
        for (m = 0; m < m_mm; m++) {
            for (k = 0; k < m_kk; k++) {
                na = m_phase->nAtoms(k,m);

                // handle the case of negative atom numbers (used to
                // represent positive ions, where the 'element' is an
                // electron
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
        // construct the chemical potentials by summing element potentials
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

    /// Estimate the initial mole numbers. This version borrows from the
    /// MultiPhaseEquil solver. 
    int ChemEquil::setInitialMoles(thermo_t& s) {
        MultiPhase* mp = 0;
        MultiPhaseEquil* e = 0;
        int iok = 0;
        beginLogGroup("ChemEquil::setInitialMoles");
        try {
            mp = new MultiPhase;
            mp->addPhase(&s, 1.0);
            mp->init();
            e = new MultiPhaseEquil(mp, true);
            e->setInitialMixMoles();

            // store component indices
            for (int m = 0; m < m_mm; m++) {
                m_component[m] = e->componentIndex(m);
            }
            for (int k = 0; k < m_kk; k++) {
                if (m_phase->moleFraction(k) > 0.0) {
                    addLogEntry(m_phase->speciesName(k),
                        m_phase->moleFraction(k));
                }
            }

            update(s);
            delete e;
            delete mp;
            iok = 0;
        }
        catch (CanteraError) {
            delete e;
            delete mp;
            iok = -1;
        }
        endLogGroup();
        return iok;
    }

 
    /**
     * Generate a starting estimate for the element potentials.
     */
    int ChemEquil::estimateElementPotentials(thermo_t& s, vector_fp& lambda) 
    {
        int m, n;
        beginLogGroup("estimateElementPotentials");
        //for (k = 0; k < m_kk; k++) {
        //    if (m_molefractions[k] > 0.0) {
        //        m_molefractions[k] = fmaxx(m_molefractions[k], 0.05);
        //    }
        //}
        //s.setState_PX(s.pressure(), m_molefractions.begin());


        DenseMatrix aa(m_mm, m_mm, 0.0);
        vector_fp b(m_mm, -999.0);

        vector_fp mu_RT(m_kk, 0.0);

        s.getChemPotentials(mu_RT.begin());
        doublereal rrt = 1.0/(GasConstant*m_phase->temperature());
        scale(mu_RT.begin(), mu_RT.end(), mu_RT.begin(), rrt);

        for (m = 0; m < m_mm; m++) {
            for (n = 0; n < m_mm; n++) {
                aa(m,n) = nAtoms(m_component[m], n);
            }
            b[m] = mu_RT[m_component[m]];
        }

        int info;
        try {
            info = solve(aa, b.begin());
        }
        catch (CanteraError) {
            addLogEntry("failed to estimate initial element potentials.");
            info = -2; 
        }

        if (info == 0) {
            for (m = 0; m < m_mm; m++) {
                lambda[m] = b[m];
                addLogEntry(m_phase->elementName(m),b[m]);
            }
        }
        endLogGroup();
        return info;
    }



    /**
     * Equilibrate a phase, holding the elemental composition fixed
     * at the initial vaollue.
     */
    int ChemEquil::equilibrate(thermo_t& s, const char* XY) {
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
    int ChemEquil::equilibrate(thermo_t& s, const char* XYstr, vector_fp& elMoles)
    {
        doublereal xval, yval;
        int fail = 0;

        delete m_p1;
        delete m_p2;
        bool tempFixed = true;
        int XY = _equilflag(XYstr);

        vector_fp state;
        s.saveState(state);

        beginLogGroup("ChemEquil::equilibrate");
        initialize(s);
        update(s);
        switch (XY) {
        case TP: case PT:
            m_p1 = new TemperatureCalculator<thermo_t>;
            m_p2 = new PressureCalculator<thermo_t>;
            break;
        case HP: case PH:
            tempFixed = false;
            m_p1 = new EnthalpyCalculator<thermo_t>;
            m_p2 = new PressureCalculator<thermo_t>;
            break; 
        case SP: case PS:
            tempFixed = false;
            m_p1 = new EntropyCalculator<thermo_t>;
            m_p2 = new PressureCalculator<thermo_t>;
            break; 
        case SV: case VS:
            tempFixed = false;
            m_p1 = new EntropyCalculator<thermo_t>;
            m_p2 = new DensityCalculator<thermo_t>;
            break; 
        case TV: case VT:
            m_p1 = new TemperatureCalculator<thermo_t>;
            m_p2 = new DensityCalculator<thermo_t>;
            break; 
        case UV: case VU:
            tempFixed = false;
            m_p1 = new IntEnergyCalculator<thermo_t>;
            m_p2 = new DensityCalculator<thermo_t>;
            break; 
        default:
            endLogGroup();
            throw CanteraError("equilibrate","illegal property pair.");
        }

        addLogEntry("Problem type","fixed "+m_p1->symbol()+", "+m_p2->symbol());
        addLogEntry(m_p1->symbol(), m_p1->value(s));
        addLogEntry(m_p2->symbol(), m_p2->value(s));

        // If the temperature is one of the specified variables, and
        // it is outside the valid range, throw an exception.
        if (tempFixed) {
            double tfixed = s.temperature();
            if (tfixed > s.maxTemp() + 1.0 || tfixed < s.minTemp() - 1.0) {
                endLogGroup();
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

        for (m = 0; m < mm; m++) {
            if (m_skip < 0 && elMoles[m] > 0.0 ) m_skip = m;
        }

        // start with a composition with everything non-zero. Note
        // that since we have already save the target element moles,
        // changing the composition at this point only affects the
        // starting point, not the final solution.
        vector_fp xmm(m_kk,0.0);
        for (int k = 0; k < m_kk; k++) {
            xmm[k] = m_phase->moleFraction(k) + Cutoff;
        }
        m_phase->setMoleFractions(xmm.begin());

        update(s);


        // loop to estimate T
        if (!tempFixed) {

            beginLogGroup("Initial T Estimate");

            doublereal tmax = m_thermo->maxTemp();
            doublereal tmin = m_thermo->minTemp();
            doublereal slope, phigh, plow, pval, dt;

            // first get the property values at the upper and lower
            // temperature limits. Since p1 (h, s, or u) is monotonic
            // in T, these values determine the upper and lower
            // bounnds (phigh, plow) for p1.

            m_phase->setTemperature(tmax);
            setInitialMoles(s); 
            phigh = m_p1->value(s);

            m_phase->setTemperature(tmin);
            setInitialMoles(s); 
            plow = m_p1->value(s);

            // start with T at the midpoint of the range
            doublereal t0 = 0.5*(tmin + tmax);
            m_phase->setTemperature(t0);

            // loop up to 5 times
            for (int it = 0; it < 5; it++) {

                // set the composition and get p1
                setInitialMoles(s);
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
                addLogEntry("new T estimate", t0);

                m_phase->setTemperature(t0);
            }
            endLogGroup(); // initial T estimate
        }

        //if (m_lambda[0] == -100.0) {
        setInitialMoles(s);
        for (int ii = 0; ii < m_mm; ii++) x[ii] = -101.0;
        estimateElementPotentials(s, x);
        //}
        //else {
        //    doublereal rt = GasConstant * m_phase->temperature();
        //    for (int ii = 0; ii < m_mm; ii++) x[ii] = m_lambda[ii]/rt;
        //}

        x[m_mm] = log(m_phase->temperature());

        vector_fp above(nvar);
        vector_fp below(nvar);

        for (m = 0; m < mm; m++) {
            above[m] = 200.0; 
            below[m] = -2000.0;
            if (elMoles[m] < Cutoff && m != m_eloc) x[m] = -1000.0;
        }
        above[mm] = log(m_thermo->maxTemp() + 1.0);
        below[mm] = log(m_thermo->minTemp() - 1.0);

        vector_fp grad(nvar, 0.0);        // gradient of f = F*F/2
        vector_fp oldx(nvar, 0.0);        // old solution
        //vector_fp prevx(nvar, 0.0);       // old solution
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
        if (iter > 1) endLogGroup(); // iteration
        beginLogGroup("Iteration "+int2str(iter));

        // compute the residual and the jacobian using the current
        // solution vector
        equilResidual(s, x, elMoles, res_trial, xval, yval);

        f = 0.5*dot(res_trial.begin(), res_trial.end(), res_trial.begin());
        addLogEntry("Residual norm", f);

        equilJacobian(s, x, elMoles, jac, xval, yval);
    
        // compute grad f = F*J
        jac.leftMult(res_trial.begin(), grad.begin());
        copy(x.begin(), x.end(), oldx.begin());

        oldf = f;
        scale(res_trial.begin(), res_trial.end(), res_trial.begin(), -1.0);
        try {
            info = solve(jac, res_trial.begin());
        }
        catch (CanteraError) {
            addLogEntry("Jacobian is singular.");
            endLogGroup(); // iteration
            endLogGroup(); // equilibrate
            s.restoreState(state);

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
                    0.8*(above[m] - x[m])/(newval - x[m])));
            }
            else if (newval < below[m]) {
                fctr = fminn(fctr, 0.8*(x[m] - below[m])/(x[m] - newval));
            }
        }
        if (fctr != 1.0) addLogEntry("factor to keep solution in bounds",
            fctr);

        // multiply the step by the scaling factor
        scale(res_trial.begin(), res_trial.end(), res_trial.begin(), fctr);

        if (!dampStep(s, oldx, oldf, grad, res_trial, 
            x, f, elMoles , xval, yval))
        {
            fail++;
            if (fail > 3) {
                addLogEntry("dampStep","Failed 3 times. Giving up.");
                endLogGroup(); // iteration
                endLogGroup(); // equilibrate
                s.restoreState(state);
                throw CanteraError("equilibrate",
                    "Cannot find an acceptable Newton damping coefficient.");
                return -4;
            }
        }
        else fail = 0;

converge:    

        //  check for convergence.
    
        equilResidual(s, x, elMoles, res_trial, xval, yval);
        f = 0.5*dot(res_trial.begin(), res_trial.end(), res_trial.begin());
        doublereal xx, yy, deltax, deltay;
        xx = m_p1->value(s);
        yy = m_p2->value(s);
        deltax = (xx - xval)/xval;
        deltay = (yy - yval)/yval;
        doublereal rmax = absmax(res_trial.begin(), res_trial.end());
        if (iter > 0 && rmax < options.relTolerance 
            && fabs(deltax) < options.relTolerance 
            && fabs(deltay) < options.relTolerance) {
            options.iterations = iter;
            endLogGroup(); // iteration
            m_lambda.resize(m_mm);
            beginLogGroup("Converged solution");
            addLogEntry("Iterations",iter);
            addLogEntry("Relative error in "+m_p1->symbol(),deltax);
            addLogEntry("Relative error in "+m_p2->symbol(),deltay);
            addLogEntry("Max residual",rmax);
            beginLogGroup("Element potentials");
            doublereal rt = GasConstant*m_thermo->temperature();
            for (m = 0; m < m_mm; m++) {
                m_lambda[m] = x[m]*rt;
                addLogEntry("element "+m_phase->elementName(m), fp2str(x[m]));
            }
            endLogGroup(); // element potentials

            if (m_thermo->temperature() > m_thermo->maxTemp() + 1.0 ||
                m_thermo->temperature() < m_thermo->minTemp() - 1.0 ) {
                writelog("Warning: Temperature ("
                    +fp2str(m_thermo->temperature())+" K) outside "
                    "valid range of "+fp2str(m_thermo->minTemp())+" K to "
                    +fp2str(m_thermo->maxTemp())+" K\n");
            }
            endLogGroup();  // converged solution

            endLogGroup();  // equilibrate
            return 0;
        }
    
        // no convergence
    
        if (iter > options.maxIterations) {
            addLogEntry("equilibrate","no convergence");
            endLogGroup(); // iteration
            endLogGroup(); // equilibrate
            s.restoreState(state);
            throw CanteraError("equilibrate",
                "no convergence in "+int2str(options.maxIterations)
                +" iterations.");
            return -1;
        }
        goto next;
    }   


    int ChemEquil::dampStep(thermo_t& mix, vector_fp& oldx, 
        double oldf, vector_fp& grad, vector_fp& step, vector_fp& x, 
        double& f, vector_fp& elmols, double xval, double yval )
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

        vector_fp res_new(nvar); // fix

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
    

        equilResidual(mix, x, elmols, res_new, xval, yval);
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
        doublereal xval, doublereal yval)
    {
        beginLogGroup("ChemEquil::equilResidual");
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
                resid[n] = log( (1.0 + elmtotal[n]) / (1.0 + elm[n]) );
            addLogEntry(m_phase->elementName(n),fp2str(elm[n])+"  ("
                +fp2str(elmtotal[n])+")");
        }
        if (m_eloc >= 0) {
            doublereal chrg, sumnet = 0.0, sumabs = 0.0;
            for (int k = 0; k < m_kk; k++) {
                chrg = m_molefractions[k]*m_phase->charge(k);
                sumnet += chrg;
                sumabs += fabs(chrg);
            }
            addLogEntry("net charge",sumnet);
            resid[m_eloc] = sumnet/m_abscharge; // log((1.0 + sumnet/sumabs));
        }
        xx = m_p1->value(mix);
        yy = m_p2->value(mix);
        resid[m_mm] = xx/xval - 1.0; 
        resid[m_skip] = yy/yval - 1.0;
	string xstr = fp2str(xx)+"  ("+fp2str(xval)+")";
        addLogEntry(m_p1->symbol(), xstr);
	string ystr = fp2str(yy)+"  ("+fp2str(yval)+")";
        addLogEntry(m_p2->symbol(), ystr);
        endLogGroup();
    }


    //-------------------- Jacobian evaluation ---------------------------

    void ChemEquil::equilJacobian(thermo_t& mix, vector_fp& x,  
        const vector_fp& elmols, DenseMatrix& jac, 
        doublereal xval, doublereal yval)
    {
        beginLogGroup("equilJacobian",0);

        int len = x.size();
        vector_fp& r0 = m_jwork1;
        vector_fp& r1 = m_jwork2;
        r0.resize(len);
        r1.resize(len);

        int n, m;
        doublereal rdx, dx, xsave;
        doublereal atol = 1.e-10;
    
        equilResidual(mix, x, elmols, r0, xval, yval);
    
        for (n = 0; n < len; n++)
        {
            // perturb x(n)
        
            xsave = x[n];
            dx = atol;
            x[n] = xsave + dx;
            dx = x[n] - xsave;
            rdx = 1.0/dx;

            // calculate perturbed residual
        
            equilResidual(mix, x, elmols, r1, xval, yval);
        
            // compute nth column of Jacobian
        
            for (m = 0; m < len; m++) {
                jac(m, n) = (r1[m] - r0[m])*rdx;
            }        
            x[n] = xsave;
        }
        endLogGroup();
    }




} // namespace

 
// $Log: ChemEquil.cpp,v
