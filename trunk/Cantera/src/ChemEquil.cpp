/**
 *
 *  @file ChemEquil.cpp
 *
 *  Chemical equilibrium.  Implementation file for class
 *  ChemEquil. 
 *
 *
 *  $Id: ChemEquil.cpp,v 1.26 2006/10/23 00:56:31 dggoodwin Exp $
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

#ifdef DEBUG_HKM_EPEQUIL
#include "stdio.h"
int debug_prnt_lvl = 0;
#endif

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
  ChemEquil::ChemEquil() : m_skip(-1), m_p1(0), m_p2(0), m_elementTotalSum(1.0), 
			   m_p0(OneAtm), m_eloc(-1),
			   m_abscharge(Tiny),
			   m_elemFracCutoff(1.0E-100),
			   m_doResPerturb(false)
  {}

  /// Destructor
  ChemEquil::~ChemEquil(){
    if (m_p1) {
      delete m_p1;
    }
    if (m_p2) {
      delete m_p2;
    }
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
        
    // allocate space in internal work arrays within the ChemEquil object
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
    // Construct the chemical potentials by summing element potentials
    fill(m_mu_RT.begin(), m_mu_RT.end(), 0.0);
    for (int k = 0; k < m_kk; k++)
      for (int m = 0; m < m_mm; m++) 
	m_mu_RT[k] += lambda_RT[m]*nAtoms(k,m);

    // Set the temperature
    s.setTemperature(t);

    // Call the phase-specific method to set the phase to the
    // equilibrium state with the specified species chemical
    // potentials.
    s.setToEquilState(DATA_PTR(m_mu_RT));
    update(s);
  }


  /** 
   *  update internally stored state information.
   * @todo argument not used.
   */
  void ChemEquil::update(const thermo_t& s) {

    // get the mole fractions, temperature, and density
    m_phase->getMoleFractions(DATA_PTR(m_molefractions));
    m_temp = m_phase->temperature();
    m_dens = m_phase->density();

    // compute the elemental mole fractions
    double sum = 0.0;
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
    // Store the sum for later use
    m_elementTotalSum = sum;
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
      /*
       * Update the current values of the temp, density, and
       * mole fraction, and element abundance vectors kept
       * within tine ChemEquil object.
       */
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


    vector_fp xMF_est(m_kk, 0.0);

    s.getMoleFractions(DATA_PTR(xMF_est));
    for (n = 0; n < s.nSpecies(); n++) {
      if (xMF_est[n] < 1.0E-20) {
	xMF_est[n] = 1.0E-20;
      }
    }
    s.setMoleFractions(DATA_PTR(xMF_est));
	


    s.getChemPotentials(DATA_PTR(mu_RT));
    doublereal rrt = 1.0/(GasConstant*m_phase->temperature());
    scale(mu_RT.begin(), mu_RT.end(), mu_RT.begin(), rrt);

#ifdef DEBUG_HKM_EPEQUIL
    if (debug_prnt_lvl > 0) {
      for (m = 0; m < m_mm; m++) {
	int isp = m_component[m];
	string nnn = s.speciesName(isp);
	printf("isp = %d, %s\n", isp, nnn.c_str());
      }
      double pres = s.pressure();
      double temp = s.temperature();
      printf("Pressure = %g\n", pres);
      printf("Temperature = %g\n", temp);
      printf("  id       Name     MF     mu/RT \n");

      s.getMoleFractions(DATA_PTR(xMF_est));

      for (n = 0; n < s.nSpecies(); n++) {
	string nnn = s.speciesName(n);
	printf("%10d %15s %10.5g %10.5g\n",
	       n, nnn.c_str(), xMF_est[n], mu_RT[n]);
      }
    }
#endif

    for (m = 0; m < m_mm; m++) {
      for (n = 0; n < m_mm; n++) {
	aa(m,n) = nAtoms(m_component[m], n);
      }
      b[m] = mu_RT[m_component[m]];
    }

    int info;
    try {
      info = solve(aa, DATA_PTR(b));
    }
    catch (CanteraError) {
      addLogEntry("failed to estimate initial element potentials.");
      info = -2; 
    }
    for (m = 0; m < m_mm; m++) {
      lambda[m] = b[m];
    }
    if (info == 0) {
      for (m = 0; m < m_mm; m++) {
	addLogEntry(m_phase->elementName(m),lambda[m]);
      }
    }

#ifdef DEBUG_HKM_EPEQUIL
    if (debug_prnt_lvl > 0) {
      for (m = 0; m < m_mm; m++) {
	int isp = m_component[m];
	double tmp = 0.0;
	for (n = 0; n < m_mm; n++) {
	  tmp += nAtoms(isp, n) * lambda[n];
	}
	printf("%3d   %10.5g   %10.5g   %10.5g\n",
	       m, mu_RT[isp], tmp, tmp - mu_RT[isp]);

      }
    }
#endif
    endLogGroup();
    return info;
  }


  /**
   * Equilibrate a phase, holding the elemental composition fixed
   * at the initial value found within the ThermoPhase object.
   *
   *   The value of 2 specified properties are obtained by querying the
   *   ThermoPhase object. The properties must be already contained
   *   within the current thermodynamic state of the system.
   */
  int ChemEquil::equilibrate(thermo_t& s, const char* XY, 
			     bool useThermoPhaseElementPotentials) {
    vector_fp emol(s.nElements());
    initialize(s);
    update(s);
    copy(m_elementmolefracs.begin(), m_elementmolefracs.end(),
	 emol.begin());
    return equilibrate(s, XY, emol, useThermoPhaseElementPotentials);
  }


  /**
   *   Compute the equilibrium composition for 2 specified
   *   properties and the specified element moles.
   *
   *   elMoles = specified vector of element abundances.
   *
   *   The 2 specified properties are obtained by querying the
   *   ThermoPhase object. The properties must be already contained
   *   within the current thermodynamic state of the system.
   *
   *   Return variable:
   *     Successful returns are indicated by a return value of 0.
   *     Unsuccessful returns are indicated by a return value of -1 for
   *     lack of convergence or -3 for a singular jacobian.
   */
  int ChemEquil::equilibrate(thermo_t& s, const char* XYstr, vector_fp& elMoles, 
			     bool useThermoPhaseElementPotentials)
  {
    doublereal xval, yval, tmp;
    int fail = 0;
    int m;

    if (m_p1) delete m_p1;
    if (m_p2) delete m_p2;
    bool tempFixed = true;
    int XY = _equilflag(XYstr);

    vector_fp state;
    s.saveState(state);

#ifdef DEBUG_HKM_EPEQUIL
    int n;
    const vector<string>& eNames = s.elementNames();
#endif
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
      endLogGroup("ChemEquil::equilibrate");
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
	endLogGroup("ChemEquil::equilibrate");
	throw CanteraError("ChemEquil","Specified temperature ("
			   +fp2str(m_thermo->temperature())+" K) outside "
			   "valid range of "+fp2str(m_thermo->minTemp())+" K to "
			   +fp2str(m_thermo->maxTemp())+" K\n");
      }                
    }

    /*
     * Before we do anything to change the ThermoPhase object,
     * we calculate and store the two specified thermodynamic
     * properties that we are after.
     */
    xval = m_p1->value(s);
    yval = m_p2->value(s);

    int mm = m_mm;
    int nvar = mm + 1;        
    DenseMatrix jac(nvar, nvar);       // jacobian
    vector_fp x(nvar, -102.0);         // solution vector
    vector_fp res_trial(nvar, 0.0);    // residual

    /*
     * Replace one of the element abundance fraction equations
     * with the specified property calculation.
     *
     * We choose the equation of the element with the highest element 
     * abundance.
     */
    tmp = -1.0;
    for (m = 0; m < mm; m++) {
      if (elMoles[m] > tmp ) {
	m_skip = m;
	tmp = elMoles[m];
      }
    }
    if (tmp <= 0.0) {
      throw CanteraError("ChemEquil",
			 "Element Abundance Vector is zeroed");
    }

    // start with a composition with everything non-zero. Note
    // that since we have already save the target element moles,
    // changing the composition at this point only affects the
    // starting point, not the final solution.
    vector_fp xmm(m_kk,0.0);
    for (int k = 0; k < m_kk; k++) {
      xmm[k] = m_phase->moleFraction(k) + 1.0E-32;
    }
    m_phase->setMoleFractions(DATA_PTR(xmm));

    /*
     * Update the internally storred values of m_temp,
     * m_dens, and the element mole fractions.
     */
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
      endLogGroup("Initial T Estimate"); // initial T estimate
    }

	
    setInitialMoles(s);
  
    /*
     * If requested, get the initial estimate for the
     * chemical potentials from the ThermoPhase object
     * itself. Or else, create our own estimate.
     */
    if (useThermoPhaseElementPotentials) {
      bool haveEm = s.getElementPotentials(DATA_PTR(x));
      if (haveEm) {
	doublereal rt = GasConstant * m_thermo->temperature();
	for (m = 0; m < m_mm; m++) {
	  x[m] /= rt;
	}
      } else {
	estimateElementPotentials(s, x);
      }
    } else {
      /*
       * Calculate initial estimates of the element potentials.
       * This algorithm uese the MultiPhaseEquil object's
       * initialization capabilities to calculate an initial
       * estimate of the mole fractions for a set of linearly
       * independent component species. Then, the element
       * potentials are solved for based on the chemical
       * potentials of the component species.
       */
      estimateElementPotentials(s, x);
    }

    /*
     * Do a better estimate of the element potentials.
     * We have found that the current estimate may not be good
     * enough to avoid drastic numerical issues associated with
     * the use of a numerically generated jacobian.
     *
     * The Brinkley algorithm assumes a constant T, P system
     * and uses a linearized analytical Jacobian that turns out
     * to be very stable.
     */
    int info = estimateEP_Brinkley(s, x, elMoles);
    if (info != 0) {
      if (info == 1) {
	addLogEntry("estimateEP_Brinkley didn't converge in given max interations");
      } else if (info == -3) {
	addLogEntry("estimateEP_Brinkley had a singular Jacobian. Continuing anyway");
      }
    } else {
      setToEquilState(s, x, m_phase->temperature());
      // Tempting -> However, nonideal is a problem. Turn on if not worried
      //             about nonideality and you are having problems with the main
      //             algorithm.
      //if (XY == TP) {
      // endLogGroup("ChemEquil::equilibrate"); 
      // return 0;
      //}
    }
 
    /*
     * Install the log(temp) into the last solution unknown
     * slot.
     */
    x[m_mm] = log(m_phase->temperature());

    /*
     * Setting the max and min values for x[]. Also, if element
     * abundance vector is zero, setting x[] to -1000. This 
     * effectively zeroes out all species containing that element.
     */
    vector_fp above(nvar);
    vector_fp below(nvar);
    for (m = 0; m < mm; m++) {
      above[m] = 200.0; 
      below[m] = -2000.0;
      if (elMoles[m] < m_elemFracCutoff && m != m_eloc) x[m] = -1000.0;
    }
    above[mm] = log(m_thermo->maxTemp() + 1.0);
    below[mm] = log(m_thermo->minTemp() - 1.0);

    vector_fp grad(nvar, 0.0);        // gradient of f = F*F/2
    vector_fp oldx(nvar, 0.0);        // old solution
    vector_fp oldresid(nvar, 0.0);
    doublereal f, oldf;

    int iter = 0;
    doublereal fctr = 1.0, newval;

    goto converge;
  next:
    // If the problem involves charged species, then the
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
    if (iter > 1) endLogGroup("Iteration "+int2str(iter-1)); // iteration
    beginLogGroup("Iteration "+int2str(iter));

    // compute the residual and the jacobian using the current
    // solution vector
    equilResidual(s, x, elMoles, res_trial, xval, yval);
    f = 0.5*dot(res_trial.begin(), res_trial.end(), res_trial.begin());
    addLogEntry("Residual norm", f);

    // Compute the Jacobian matrix
    equilJacobian(s, x, elMoles, jac, xval, yval);
    
#ifdef DEBUG_HKM_EPEQUIL
    if (debug_prnt_lvl > 0) {
      printf("Jacobian matrix %d:\n", iter);
      for (m = 0; m <= m_mm; m++) {
	printf("      [ ");
	for (n = 0; n <= m_mm; n++) {
	  printf("%10.5g ", jac(m,n));
	}
	printf(" ]");
	char xName[32];
	if (m < m_mm) {
	  string nnn = eNames[m];
	  sprintf(xName, "x_%-10s", nnn.c_str());
	} else {
	  sprintf(xName, "x_XX");
	}
	if (m_eloc == m) {
	  sprintf(xName, "x_ELOC");
	}
	if (m == m_skip) {
	  sprintf(xName, "x_YY");
	}
	printf("%-12s", xName);
	printf(" =  - (%10.5g)\n", res_trial[m]);
      }
    }
#endif


    // compute grad f = F*J
    jac.leftMult(DATA_PTR(res_trial), DATA_PTR(grad));
    copy(x.begin(), x.end(), oldx.begin());
    oldf = f;
    scale(res_trial.begin(), res_trial.end(), res_trial.begin(), -1.0);

    /*
     * Solve the system
     */
    try {
      info = solve(jac, DATA_PTR(res_trial));
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
    doublereal rmax = 0.0;
    bool passThis = true;
    for (m = 0; m < nvar; m++) {
      double tval =  options.relTolerance;
      if (m < mm) {
	tval = elMoles[m] * options.relTolerance + options.absElemTol;
      }
      if (fabs(res_trial[m]) > tval) {
	passThis = false;
      }
    }
    if (iter > 0 && passThis
	&& fabs(deltax) < options.relTolerance 
	&& fabs(deltay) < options.relTolerance) {
      options.iterations = iter;
      endLogGroup("Iteration "+int2str(iter)); // iteration
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
      /*
       * Save the calculated and converged element potentials
       * to the original ThermoPhase object.
       */
      s.setElementPotentials(m_lambda);
      addLogEntry("Saving Element Potentials to ThermoPhase Object");
      endLogGroup("Element potentials");

      if (m_thermo->temperature() > m_thermo->maxTemp() + 1.0 ||
	  m_thermo->temperature() < m_thermo->minTemp() - 1.0 ) {
	writelog("Warning: Temperature ("
		 +fp2str(m_thermo->temperature())+" K) outside "
		 "valid range of "+fp2str(m_thermo->minTemp())+" K to "
		 +fp2str(m_thermo->maxTemp())+" K\n");
      }
      endLogGroup("Converged solution");  
      endLogGroup("ChemEquil::equilibrate"); 
      return 0;
    }
    
    // no convergence
    
    if (iter > options.maxIterations) {
      addLogEntry("equilibrate","no convergence");
      endLogGroup("Iteration "+int2str(iter)); 
      endLogGroup("ChemEquil::equilibrate");
      s.restoreState(state);
      throw CanteraError("equilibrate",
			 "no convergence in "+int2str(options.maxIterations)
			 +" iterations.");
      return -1;
    }
    goto next;
  }   


  /*
   * dampStep: Come up with an acceptable step size. The original implementation
   *           employed a line search technique that enforced a reduction in the
   *           norm of the residual at every successful step. Unfortunately,
   *           this method created false convergence errors near the end of 
   *           a significant number of steps, usually special conditions where
   *           there were stoichiometric constraints. 
   *   
   *  This new method just does a delta damping approach, based on limiting
   *  the jump in the dimensionless element potentials. Mole fractions are
   *  limited to a factor of 2 jump in the values from this method.
   *  Near convergence, the delta damping gets out of the way.
   */
  int ChemEquil::dampStep(thermo_t& mix, vector_fp& oldx, 
			  double oldf, vector_fp& grad, vector_fp& step, vector_fp& x,
			  double& f, vector_fp& elmols, double xval, double yval )
  {
    int nvar = x.size();
    int m;
    double damp;
    
    /*
     * Carry out a delta damping approach on the dimensionless element potentials.
     */
    damp = 1.0;
    for (m = 0; m < m_mm; m++) {
      if (step[m] > 0.75) {
	damp = 0.75 /step[m];
      }
      if (step[m] < -0.75) {
	damp = -0.75 / step[m];
      }
    }

    /*
     * Update the solution unknown
     */
    for (m = 0; m < nvar; m++) {
      x[m] = oldx[m] + damp * step[m];
    }
#ifdef DEBUG_HKM_EPEQUIL
    if (debug_prnt_lvl > 0) {
      printf("Solution Unknowns: damp = %g\n", damp);
      printf("            X_new      X_old       Step\n");
      for (m = 0; m < nvar; m++) {
	printf("     %10.5g   %10.5g    %10.5g\n", x[m], oldx[m], step[m]);
      }
    }
#endif
    return 1;
  }


  /**
   *  evaluates the residual vector F, of length mm 
   */
  void ChemEquil::equilResidual(thermo_t& mix, const vector_fp& x, 
				const vector_fp& elmFracGoal, vector_fp& resid, 
				doublereal xval, doublereal yval)
  {
    beginLogGroup("ChemEquil::equilResidual");
    int n;
    doublereal xx, yy;
    doublereal temp = exp(x[m_mm]);
    setToEquilState(mix, x, temp);

    // residuals are the total element moles
    vector_fp& elmFrac = m_elementmolefracs; 
    for (n = 0; n < m_mm; n++)
      {
	// drive element potential for absent elements to -1000
	if (elmFracGoal[n] < m_elemFracCutoff && n != m_eloc)
	  resid[n] = x[n] + 1000.0;
	else {
	  /*
	   * Change the calculation for small element number, using
	   * L'Hopital's rule.
	   * The log formulation is unstable.
	   */
	  if (elmFracGoal[n] < 1.0E-10 || elmFrac[n] < 1.0E-10) {
	    resid[n] = elmFracGoal[n] - elmFrac[n];
	  } else {
	    resid[n] = log( (1.0 + elmFracGoal[n]) / (1.0 + elmFrac[n]) );
	  }
	}
	addLogEntry(m_phase->elementName(n),fp2str(elmFrac[n])+"  ("
		    +fp2str(elmFracGoal[n])+")");
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
    endLogGroup("ChemEquil::equilResidual");

#ifdef DEBUG_HKM_EPEQUIL
    if (debug_prnt_lvl > 0 && !m_doResPerturb) {
      printf("Residual: ElFracGoal   ElFracCurrent   Resid\n");
      for (n = 0; n < m_mm; n++) {
	printf("    %14.9g %14.9g %10.5g\n", elmFracGoal[n], elmFrac[n], resid[n]);
      }
    }
#endif
  }


  //-------------------- Jacobian evaluation ---------------------------

  void ChemEquil::equilJacobian(thermo_t& s, vector_fp& x,  
				const vector_fp& elmols, DenseMatrix& jac, 
				doublereal xval, doublereal yval)
  {
    beginLogGroup("equilJacobian");
    int len = x.size();
    vector_fp& r0 = m_jwork1;
    vector_fp& r1 = m_jwork2;
    r0.resize(len);
    r1.resize(len);
    int n, m;
    doublereal rdx, dx, xsave, dx2;
    doublereal atol = 1.e-10;
    
    equilResidual(s, x, elmols, r0, xval, yval);
    
    m_doResPerturb = true;
    for (n = 0; n < len; n++) {
      xsave = x[n];
      dx = atol;
      dx2 = fabs(xsave) * 1.0E-7;
      if (dx2 > dx) dx = dx2;
      x[n] = xsave + dx;
      dx = x[n] - xsave;
      rdx = 1.0/dx;

      // calculate perturbed residual
        
      equilResidual(s, x, elmols, r1, xval, yval);
        
      // compute nth column of Jacobian
        
      for (m = 0; m < len; m++) {
	jac(m, n) = (r1[m] - r0[m])*rdx;
      }        
      x[n] = xsave;
    }
    m_doResPerturb = false;
    endLogGroup("equilJacobian");
  }

  /**
   *  Do a calculation of the element potentials using
   *  the Brinkley method, p. 129 Smith and Missen.
   *
   * We have found that the previous estimate may not be good
   * enough to avoid drastic numerical issues associated with
   * the use of a numerically generated jacobian used in the
   * main algorithm.
   *
   * The Brinkley algorithm, here, assumes a constant T, P system
   * and uses a linearized analytical Jacobian that turns out
   * to be very stable even given bad initial guesses.
   *
   * The pressure and temperature to be used are in the
   * ThermoPhase object input into the routine.
   *
   * The initial guess for the element potentials
   * used by this routine is taken from the
   * input vector, x.
   *
   * elMoles is the input element abundance vector to be matched.
   *
   * Nonideal phases are handled in principle. This is done by 
   * calculating the activity coefficients and adding them
   * into the formula in the correct position. However, 
   * these are treated as a rhs contribution only. Therefore,
   * convergence might be a problem. This has not been tested.
   * Also molality based unit systems aren't handled.
   *
   * On return, int return value contains the success code:
   *    0 - successful
   *    1 - unsuccessful, max num iterations exceeded
   *   -3 - unsuccessful, singular jacobian
   *
   * NOTE: update for activity coefficients.
   */
  int ChemEquil::estimateEP_Brinkley(thermo_t& s, vector_fp& x, 
				     vector_fp& elMoles) {
    /*
     * Before we do anything, we will save the state of the solution.
     * Then, if things go drastically wrong, we will restore the 
     * saved state.
     */
    vector_fp state;
    s.saveState(state);
    double tmp, sum;
    bool modifiedMatrix = false;
    int neq = m_mm+1;
    int retn = 1;
    int m, n, k, info;
    DenseMatrix a1(neq, neq, 0.0);
    vector_fp b(neq, 0.0);
    vector_fp muSS_RT(m_kk, 0.0);
    vector_fp n_i(m_kk,0.0);
    vector_fp n_i_calc(m_kk,0.0);
    vector_fp actCoeff(m_kk, 1.0);
    vector_fp muSS_RT_mod(m_kk, 0.0); 
    double beta = 1.0;

    s.getMoleFractions(DATA_PTR(n_i));

    vector_fp x_old(m_mm+1, 0.0);
    vector_fp resid(m_mm+1, 0.0);
    vector_int lumpSum(m_mm+1, 0);

    /*
     * Get the nondimensional Gibbs functions for the species
     * at their standard states of solution at the current T and P
     * of the solution.
     */
    s.getGibbs_RT(DATA_PTR(muSS_RT));
    copy(muSS_RT.begin(), muSS_RT.end(), muSS_RT_mod.begin());


    vector_fp eMolesCalc(m_mm, 0.0);
    vector_fp eMolesFix(m_mm, 0.0);
    double elMolesTotal = 0.0;
    for (m = 0; m < m_mm; m++) {
      elMolesTotal += elMoles[m];
      for (k = 0; k < m_kk; k++) {
	eMolesFix[m] +=  nAtoms(k,m) * n_i[k];
      }
    }

    for (m = 0; m < m_mm; m++) {
      if (x[m] > 50.0) {
	x[m] = 50.;
      }
      if (elMoles[m] > 1.0E-70) {
	if (x[m] < -100) {
	  x[m] = -100.;
	}
      } else {
	if (x[m] < -1000.) {
	  x[m] = -1000.;
	}
      }
    }


    double n_t = 0.0;
    double sum2 = 0.0;
    double nAtomsMax = 1.0;
    s.setMoleFractions(DATA_PTR(x));
    s.getActivityCoefficients(DATA_PTR(actCoeff));
    for (k = 0; k < m_kk; k++) {
      tmp = - (muSS_RT[k] + log(actCoeff[k]));
      sum2 = 0.0;
      for (m = 0; m < m_mm; m++) {
	sum = nAtoms(k,m);
	tmp += sum * x[m];
	sum2 += sum;
	if (sum2 > nAtomsMax) {
	  nAtomsMax = sum2;
	}
      }
      if (tmp > 100.) {
	n_t += 2.8E43;
      } else {
	n_t += exp(tmp);
      }
    }
   

#ifdef DEBUG_HKM_EPEQUIL
    const vector<string>& eNames = s.elementNames();
    if (debug_prnt_lvl > 0) {
      printf("estimateEP_Brinkley::\n\n");
      double temp = s.temperature();
      double pres = s.pressure();
      printf("temp = %g\n", temp);
      printf("pres = %g\n", pres);
      printf("Initial mole numbers and mu_SS:\n"); 
      printf("         Name           MoleNum        mu_SS   actCoeff\n");
      for (k = 0; k < m_kk; k++) {
	string nnn = s.speciesName(k);
	printf("%15s  %13.5g  %13.5g %13.5g\n",
	       nnn.c_str(), n_i[k], muSS_RT[k], actCoeff[k]);
      }
      printf("Initial n_t = %10.5g\n", n_t);
      printf("Comparison of Goal Element Abundance with Initial Guess:\n");
      printf("  eName       eCurrent       eGoal\n");
      for (m = 0; m < m_mm; m++) {
	string nnn = s.elementName(m);
	printf("%5s   %13.5g  %13.5g\n",nnn.c_str(), eMolesFix[m], elMoles[m]); 
      }
    }
#endif
    for (m = 0; m < m_mm; m++) {
      if (m != m_eloc) {
	if (elMoles[m] <= options.absElemTol) {
	  x[m] = -200.;
	}
      }
    }
    /*
     * -------------------------------------------------------------------
     * Main Loop.
     */
    for (int iter = 0; iter < 2* options.maxIterations; iter++) {
      /*
       * Save the old solution
       */
      for (m = 0; m < m_mm; m++) {
	x_old[m] = x[m];
      }
      x_old[m_mm] = n_t;
      /*
       * Calculate the mole numbers of species
       */
#ifdef DEBUG_HKM_EPEQUIL
      if (debug_prnt_lvl > 0) {
	printf("START ITERATION %d:\n", iter);
      }
#endif
    
      double n_t_calc = 0.0;
      s.setMoleFractions(DATA_PTR(x));
      s.getActivityCoefficients(DATA_PTR(actCoeff));
      sum2 = 0.0;
      for (k = 0; k < m_kk; k++) {
	tmp = - (muSS_RT[k] + log(actCoeff[k]));
	for (m = 0; m < m_mm; m++) {
	  tmp += nAtoms(k,m) * x[m];
	}
	if (tmp > 100.) tmp = 100.;
	if (tmp < -300.) {
	  n_i_calc[k] = 0.0;
	} else {
	  n_i_calc[k] = n_t * exp(tmp);
	}
	n_t_calc +=  n_i_calc[k];
#ifdef DEBUG_HKM_EPEQUIL
	if (debug_prnt_lvl > 0) {
	  string nnn = s.speciesName(k);
	  printf("%15s: %10.5g (%10.5g)\n", nnn.c_str(), 
		 n_i_calc[k], tmp);
	}
#endif
      }
#ifdef DEBUG_HKM_EPEQUIL
      if (debug_prnt_lvl > 0) {
	printf("%15s: %10.5g\n", "Total Molar Sum", n_t_calc);
      }
#endif
      for (m = 0; m < m_mm; m++) {
	eMolesCalc[m] = 0.0;
	for (k = 0; k < m_kk; k++) {
	  eMolesCalc[m] += nAtoms(k,m) * n_i_calc[k];
	}
      }
#ifdef DEBUG_HKM_EPEQUIL
      if (debug_prnt_lvl > 0) {
	printf("(iter %d) element moles bal:   Goal  Calculated\n", iter);
	for (m = 0; m < m_mm; m++) {
	  string nnn = eNames[m];
	  printf("              %8s: %10.5g %10.5g \n", nnn.c_str(), elMoles[m], eMolesCalc[m]);
	}
      }
#endif

      double nCutoff;

      bool normalStep = true;
      /*
       * Decide if we are to do a normal step or a modified step
       */
      int iM = -1;
      for (m = 0; m < m_mm; m++) {
	if (elMoles[m] > 0.001 * elMolesTotal) {
	  if (eMolesCalc[m] > 1000. * elMoles[m]) {
	    normalStep = false;
	    iM = m;
	  }
	  if (1000 * eMolesCalc[m] < elMoles[m]) {
	    normalStep = false;
	    iM = m;
	  }
	}
      }
#ifdef DEBUG_HKM_EPEQUIL
      if (debug_prnt_lvl > 0) {
	if (!normalStep) {
	  printf(" NOTE: iter(%d) Doing an abnormal step due to row %d\n", iter, iM); 
	}
      }
#endif 
      if (!normalStep) {
	beta = 1.0;
	resid[m_mm] = 0.0;
	for (m = 0; m < m_mm; m++) {
	  resid[m] = 0.0;
	  if (elMoles[m] > 0.001 * elMolesTotal) {
	    if (eMolesCalc[m] > 1000. * elMoles[m]) {
	      resid[m] = -0.5;
	      resid[m_mm] -= 0.5;
	    }
	    if (1000 * eMolesCalc[m] < elMoles[m]) {
	      resid[m] = 0.5;
	      resid[m_mm] += 0.5;
	    }
	  }
	}
	if (n_t < (elMolesTotal / nAtomsMax)) {
	  if (resid[m_mm] < 0.0) {
	    resid[m_mm] = 0.1;
	  }
	} else if (n_t > elMolesTotal) {
	  if (resid[m_mm] > 0.0) {
	    resid[m_mm] = 0.0;
	  }
	}
	goto updateSolnVector;
      }
   

      /*
       * Determine whether the matrix should be dumbed down because 
       * the coefficient matrix of species (with significant concentrations)
       * is rank deficient.
       *
       * The basic idea is that at any time during the calculation only a
       * small subset of species with sufficient concentration matters.
       * If the rank of the element coefficient matrix for that subset of species 
       * is less than the number of elements, then the matrix created by
       * the Brinkley method below may become singular.
       *
       * The logic below looks for obvious cases where the current element
       * coefficient matrix is rank deficient.
       * 
       * The way around rank-deficiency is to lump-sum the corresponding row 
       * of the matrix. Note, lump-summing seems to work very well in terms of 
       * its stability properties, i.e., it heads in the right direction,
       * albeit with lousy convergence rates.
       *
       * NOTE: This probably should be extended to a full blown Gauss-Jordon
       *       factorization scheme in the future. For Example
       *       the scheme below would fail for the set: HCl  NH4Cl, NH3.
       *       Hopefully, it's caught by the equal rows logic below.
       */
      for (m = 0; m < m_mm; m++) {
	lumpSum[m] = 1;
      }

      nCutoff = 1.0E-4 * n_t_calc;
#ifdef DEBUG_HKM_EPEQUIL
      if (debug_prnt_lvl > 0) {
        printf(" Lump Sum Elements Calculation: \n");
      }
#endif
      for (m = 0; m < m_mm; m++) {
	int kMSp = -1;
	int kMSp2 = -1;
	int nSpeciesWithElem  = 0;
	for (k = 0; k < m_kk; k++) {
	  if (n_i_calc[k] > nCutoff) {
	    if (nAtoms(k,m) > 0) {
	      nSpeciesWithElem++;
	      if (kMSp != -1) {
		kMSp2 = k;
		double factor = nAtoms(kMSp,m) / nAtoms(kMSp2,m);
		for (n = 0; n < m_mm; n++) {
		  if (fabs(factor *  nAtoms(kMSp2,n) -  nAtoms(kMSp,n)) > 1.0E-8) {
		    lumpSum[m] = 0;
		    break;
		  }
		}
	      } else {
		kMSp = k;
	      }
	    }
	  }
	}
#ifdef DEBUG_HKM_EPEQUIL
	if (debug_prnt_lvl > 0) {
	  string nnn = eNames[m];
	  printf("               %5s %3d : %5d  %5d\n",nnn.c_str(), lumpSum[m], kMSp, kMSp2); 
	  
	}
#endif
      }



      for (m = 0; m < m_mm; m++) {
	for (n = 0; n < m_mm; n++) {
	  a1(m,n) = 0.0;
	  for (k = 0; k < m_kk; k++) {
	    a1(m,n) += nAtoms(k,m) * nAtoms(k,n) * n_i_calc[k];
	  }
	}
	a1(m,m_mm) = eMolesCalc[m];
	a1(m_mm, m) = eMolesCalc[m];
	a1(m_mm, m_mm) = 0.0;      
      }
 
      sum = 0.0;
      for (m = 0; m < m_mm; m++) {
	resid[m] = elMoles[m] - eMolesCalc[m];
	tmp = resid[m] / (elMoles[m] + options.absElemTol);
	sum += tmp * tmp;
      }

      for (m = 0; m < m_mm; m++) {
	if (a1(m,m) < 1.0E-50) {
#ifdef DEBUG_HKM_EPEQUIL
	  if (debug_prnt_lvl > 0) {
	    printf(" NOTE: Diagonalizing the analytical Jac row %d\n", m); 
	  }
#endif
	  for (n = 0; n < m_mm; n++) {
	    a1(m,n) = 0.0;
	  }
	  a1(m,m) = 1.0;
	  if (resid[m] > 0.0) {
	    resid[m] = 1.0;
	  } else if (resid[m] < 0.0) {
	    resid[m] = -1.0;
	  } else {
	    resid[m] = 0.0;
	  }
	}
      }


      resid[m_mm] = n_t - n_t_calc;

#ifdef DEBUG_HKM_EPEQUIL
      if (debug_prnt_lvl > 0) {
	printf("Matrix:\n");
	for (m = 0; m <= m_mm; m++) {
	  printf("       [");
	  for (n = 0; n <= m_mm; n++) {
	    printf(" %10.5g", a1(m,n));
	  }
	  printf("]  =   %10.5g\n", resid[m]);
	}
      }
#endif
      
      tmp = resid[m_mm] /(n_t + 1.0E-15);
      sum += tmp * tmp;
#ifdef DEBUG_HKM_EPEQUIL
      if (debug_prnt_lvl > 0) {
	printf("(it %d) Convergence = %g\n", iter, sum);
      }
#endif
      if (sum < 100. * options.relTolerance) {
	retn = 0;
	goto exit;
      }


      /*
       * Row Sum scaling
       */
      for (m = 0; m <= m_mm; m++) {
	tmp = 0.0;
	for (n = 0; n <= m_mm; n++) {
	  tmp += fabs(a1(m,n));
	}
	if (m < m_mm && tmp < 1.0E-30) {
#ifdef DEBUG_HKM_EPEQUIL
	  if (debug_prnt_lvl > 0) {
	    printf(" NOTE: Diagonalizing row %d\n", m); 
	  }
#endif
	  for (n = 0; n <= m_mm; n++) {
	    if (n != m) {
	      a1(m,n) = 0.0;
	      a1(n,m) = 0.0;
	    }
	  }
	}
	tmp = 1.0/tmp;
	for (n = 0; n <= m_mm; n++) {
	  a1(m,n) *= tmp;
	}
	resid[m] *= tmp;
      }

#ifdef DEBUG_HKM_EPEQUIL
      if (debug_prnt_lvl > 0) {
	printf("Row Summed Matrix:\n");
	for (m = 0; m <= m_mm; m++) {
	  printf("       [");
	  for (n = 0; n <= m_mm; n++) {
	    printf(" %10.5g", a1(m,n));
	  }
	  printf("]  =   %10.5g\n", resid[m]);
	}
      }
#endif

      /*
       * Next Step: We have row-summed the equations.
       * However, there are some degenerate cases where two
       * rows will be multiplies of each other in terms of 
       * 0 < m, 0 < m part of the matrix. This occurs on a case
       * by case basis, and depends upon the current state of the
       * element potential values, which affect the concentrations
       * of species. 
       * So, the way we have found to eliminate this problem is to
       * lump-sum one of the rows of the matrix, except for the
       * last column, and stick it all on the diagonal. 
       * Then, we at least have a non-singular matrix, and the
       * modified equation moves the corresponding unknown in the
       * correct direction.
       * The previous row-sum operation has made the identification
       * of identical rows much simpler.
       *
       * Note at least 6E-4 is necessary for the comparison.
       * I'm guessing 1.0E-3. If two rows are anywhere close to being
       * equivalent, the algorithm can get stuck in an oscillatory mode.
       */
      modifiedMatrix = false;
      for (m = 0; m < m_mm; m++) {
	int sameAsRow = -1;
	for (int im = 0; im < m; im++) {
	  bool theSame = true;
	  for (n = 0; n < m_mm; n++) {
	    if (fabs(a1(m,n) - a1(im,n)) > 1.0E-3) {
	      theSame = false;
	      break;
	    }
	  }
	  if (theSame) {
	    sameAsRow = im;
	  }
	}
	if (sameAsRow >= 0 || lumpSum[m]) {
#ifdef DEBUG_HKM_EPEQUIL
	  if (debug_prnt_lvl > 0) {
	    if (lumpSum[m]) {
	      printf("Lump summing row %d, due to rank deficiency analysis\n", m);
	    } else if (sameAsRow >= 0) {
	      printf("Identified that rows %d and %d are the same\n", m, sameAsRow);
	    }
	  }
#endif
	  modifiedMatrix = true;
	  for (n = 0; n < m_mm; n++) {
	    if (n != m) {
	      a1(m,m) += fabs(a1(m,n));
	      a1(m,n) = 0.0;
	    }
	  }
	}
      }

#ifdef DEBUG_HKM_EPEQUIL
      if (debug_prnt_lvl > 0 && modifiedMatrix) {
	printf("Row Summed, MODIFIED Matrix:\n");
	for (m = 0; m <= m_mm; m++) {
	  printf("       [");
	  for (n = 0; n <= m_mm; n++) {
	    printf(" %10.5g", a1(m,n));
	  }
	  printf("]  =   %10.5g\n", resid[m]);
	}
      }
#endif

      try {
	info = solve(a1, DATA_PTR(resid));
      }
      catch (CanteraError) {
	addLogEntry("estimateEP_Brinkley:Jacobian is singular.");
#ifdef DEBUG_HKM_EPEQUIL
	if (debug_prnt_lvl > 0) {
	  printf("Matrix is SINGULAR.ERROR\n");
	}
#endif
	s.restoreState(state);
	throw CanteraError("equilibrate:estimateEP_Brinkley()",
			   "Jacobian is singular. \nTry adding more species, "
			   "changing the elemental composition slightly, \nor removing "
			   "unused elements.");
	return -3;
      }

      /*
       * Figure out the damping coefficient: Use a delta damping
       * coefficient formulation: magnitude of change is capped
       * to exp(1).
       */
      beta = 1.0;
      for (m = 0; m < m_mm; m++) {
	if (resid[m] > 1.0) {
	  double betat = 1.0 / resid[m];
	  if (betat < beta) {
	    beta = betat;
	  }
	}
	if (resid[m] < -1.0) {
	  double betat = -1.0 / resid[m];
	  if (betat < beta) {
	    beta = betat;
	  }
	}	
      }
#ifdef DEBUG_HKM_EPEQUIL
      if (debug_prnt_lvl > 0) {
	if (beta != 1.0) {
	  printf("(it %d) Beta = %g\n", iter, beta);
	}
      }
#endif

      /*
       * Update the solution vector
       */
    updateSolnVector:
      for (m = 0; m < m_mm; m++) {
	x[m] += beta * resid[m];
      }
      n_t *= exp(beta * resid[m_mm]);  
    

#ifdef DEBUG_HKM_EPEQUIL
      if (debug_prnt_lvl > 0) {
	printf("(it %d)    OLD_SOLUTION  NEW SOLUTION    (undamped updated)\n", iter);
	for (m = 0; m < m_mm; m++) {
	  string eee = eNames[m];
	  printf("     %5s   %10.5g   %10.5g   %10.5g\n", eee.c_str(), x_old[m], x[m], resid[m]);
	}
        printf("       n_t    %10.5g   %10.5g  %10.5g \n",  x_old[m_mm], n_t, exp(resid[m_mm]));
      }
#endif
    }
  exit:
    return retn;
  }

} // namespace
