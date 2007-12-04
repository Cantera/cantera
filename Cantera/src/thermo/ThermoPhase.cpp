/**
 *  @file ThermoPhase.cpp
 * Definition file for class ThermoPhase, the base class for phases with
 * thermodynamic properties
 * (see class \link Cantera::ThermoPhase ThermoPhase\endlink).
 */

/*
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright 2002 California Institute of Technology
 *
 */

// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ThermoPhase.h"

#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif

using namespace std;

namespace Cantera {

    //! Constructor. Note that ThermoPhase is meant to be used as
    //! a base class, so this constructor should not be called
    //! explicitly.
    ThermoPhase::ThermoPhase() :
         Phase(), 
         m_spthermo(0), m_speciesData(0),
         m_index(-1), 
         m_phi(0.0), 
         m_hasElementPotentials(false),
         m_chargeNeutralityNecessary(false)
    {
    }

   ThermoPhase::~ThermoPhase() 
   {
     delete m_spthermo;
   }
    /**
     * Copy Constructor for the ThermoPhase object. 
     *
     * Currently, this is implemented, but not tested. If called it will
     * throw an exception until fully tested.
     */
    ThermoPhase::ThermoPhase(const ThermoPhase &right)  :
	Phase(),
	m_spthermo(0), 
	m_speciesData(0),
	m_index(-1),
	m_phi(0.0), 
	m_hasElementPotentials(false),
        m_chargeNeutralityNecessary(false)
    {
      /*
       * Call the assignment operator
       */
      *this = operator=(right);
    }
  
    /*
     * operator=()
     *
     *  Note this stuff will not work until the underlying phase
     *  has a working assignment operator
     */
    ThermoPhase& ThermoPhase::
    operator=(const ThermoPhase &right) {
	/*
         * Check for self assignment.
         */
        if (this == &right) return *this;
        /*
         * Call the base class assignment operator
         */
	(void)Phase::operator=(right);

	/*
	 * Pointer to the species thermodynamic property manager
	 * We own this, so we need to do a deep copy
	 */
	if (m_spthermo) {
	  delete m_spthermo;
	}
        m_spthermo = (right.m_spthermo)->duplMyselfAsSpeciesThermo();

	// We don't do a deep copy here, because we don't own this
        m_speciesData = right.m_speciesData;
      
	m_index = right.m_index;
        m_phi = right.m_phi;
	m_lambdaRRT = right.m_lambdaRRT;
	m_hasElementPotentials = right.m_hasElementPotentials;
	m_chargeNeutralityNecessary = right.m_chargeNeutralityNecessary;

	return *this;
    }

    /*
     * Duplication routine for objects which inherit from 
     * ThermoPhase.
     *
     *  This virtual routine can be used to duplicate thermophase objects
     *  inherited from ThermoPhase even if the application only has
     *  a pointer to ThermoPhase to work with.
     * 
     *  Currently, this is not fully implemented. If called, an
     *  exception will be called by the ThermoPhase copy constructor.
     */
    ThermoPhase *ThermoPhase::duplMyselfAsThermoPhase() const {
	ThermoPhase* tp = new ThermoPhase(*this);
	return tp;
    }

    int ThermoPhase::activityConvention() const {
	return cAC_CONVENTION_MOLAR;
    }

    void ThermoPhase::getActivities(doublereal* a) const {
        getActivityConcentrations(a);
        int nsp = nSpecies();
        int k;
        for (k = 0; k < nsp; k++) a[k] /= standardConcentration(k);
    }

    void ThermoPhase::setState_TPX(doublereal t, doublereal p, 
        const doublereal* x) {
        setMoleFractions(x); setTemperature(t); setPressure(p);
    }

    void ThermoPhase::setState_TPX(doublereal t, doublereal p, 
        compositionMap& x) {
        setMoleFractionsByName(x); setTemperature(t); setPressure(p);
    }

    void ThermoPhase::setState_TPX(doublereal t, doublereal p, 
        const std::string& x) {
        compositionMap xx;
        int kk = nSpecies();
        for (int k = 0; k < kk; k++) xx[speciesName(k)] = -1.0;
        try {
            parseCompString(x, xx);
        }
        catch (CanteraError) {
            throw CanteraError("setState_TPX",
                "Unknown species in composition map: "+ x);
        }
        setMoleFractionsByName(xx); setTemperature(t); setPressure(p);
    }        

    void ThermoPhase::setState_TPY(doublereal t, doublereal p, 
        const doublereal* y) {
        setMassFractions(y); setTemperature(t); setPressure(p);
    }

    void ThermoPhase::setState_TPY(doublereal t, doublereal p, 
        compositionMap& y) {
        setMassFractionsByName(y); setTemperature(t); setPressure(p);
    }
        
    void ThermoPhase::setState_TPY(doublereal t, doublereal p, 
        const std::string& y) {
        compositionMap yy;
        int kk = nSpecies();
        for (int k = 0; k < kk; k++) yy[speciesName(k)] = -1.0;
        try {
            parseCompString(y, yy);
        }
        catch (CanteraError) {
            throw CanteraError("setState_TPY",
                "Unknown species in composition map: "+ y);
        }
        setMassFractionsByName(yy); setTemperature(t); setPressure(p);
    }

    void ThermoPhase::setState_TP(doublereal t, doublereal p) {
        setTemperature(t); setPressure(p);
    }

    void ThermoPhase::setState_PX(doublereal p, doublereal* x) {
        setMoleFractions(x); setPressure(p);
    }

    void ThermoPhase::setState_PY(doublereal p, doublereal* y) {
        setMassFractions(y); setPressure(p);
    }

  void ThermoPhase::setState_HP(doublereal Htarget, doublereal p, 
				doublereal dTtol) {
    setState_HPorUV(Htarget, p, dTtol, false);
  }

  void ThermoPhase::setState_UV(doublereal u, doublereal v, 
				doublereal dTtol) {
    setState_HPorUV(u, v, dTtol, true);
  }

  void ThermoPhase::setState_HPorUV(doublereal Htarget, doublereal p, 
				    doublereal dTtol, bool doUV) {
    doublereal dt;
    doublereal Hmax = 0.0, Hmin = 0.0;;
    doublereal v = 0.0;
    if (doUV) {
      v = p;
      setDensity(1.0/v);
    } else {
      setPressure(p);
    }
    double Tmax = maxTemp() + 0.1;
    double Tmin = minTemp() - 0.1;

    // Make sure we are within the temperature bounds at the start
    // of the iteration
    double Tnew = temperature();
    if (Tnew > Tmax) {
      Tnew = Tmax - 1.0;
      if (doUV) {
	setTemperature(Tnew);
      } else {
	setState_TP(Tnew, p);
      }
    }
    if (Tnew < Tmin) {
      Tnew = Tmin + 1.0;
      if (doUV) {
	setTemperature(Tnew);
      } else {
	setState_TP(Tnew, p);
      }
    }

    double Hnew = 0.0;
    double Cpnew = 0.0;
    if (doUV) {
      Hnew = intEnergy_mass();
      Cpnew = cv_mass();
    } else {
      Hnew = enthalpy_mass();
      Cpnew = cp_mass();
    }
    double Htop = Hnew;
    double Ttop = Tnew;
    double Hbot = Hnew;
    double Tbot = Tnew;
    double Told = Tnew;
    double Hold = Hnew;

    bool ignoreBounds = false;
    // Unstable phases are those for which
    // cp < 0.0. These are possible for cases where
    // we have passed the spinodal curve.
    bool unstablePhase = false;
    bool unstablePhaseNew = false;
   

    // Newton iteration
    for (int n = 0; n < 500; n++) {
      Told = Tnew;
      Hold = Hnew;
      double cpd = Cpnew;
      if (cpd < 0.0) {
	unstablePhase = true;
      }
      dt = (Htarget - Hold)/cpd;

      // limit step size to 210 K
      if (dt > 100.0)       dt =  100.0;
      else if (dt < -100.0) dt = -100.0; 

      // Calculate the new T
      Tnew = Told + dt;

      // Limit the step size so that we are convergent
      // This is the step that makes it different from a 
      // Newton's algorithm
      if (dt > 0.0) {
	if (!unstablePhase) {
	  if (Htop > Htarget) {
	    if (Tnew > (0.75 * Ttop + 0.25 * Told)) {
	      dt = 0.75 * (Ttop - Told);
	      Tnew = Told + dt;
	    }
	  }
	} else {
	  if (Hbot < Htarget) {
	    if (Tnew < (0.75 * Tbot + 0.25 * Told)) {
	      dt = 0.75 * (Tbot - Told);
	      Tnew = Told + dt;
	    }
	  }
	}
      } else {
	if (!unstablePhase) {
	  if (Hbot < Htarget) {
	    if (Tnew < (0.75 * Tbot + 0.25 * Told)) {
	      dt = 0.75 * (Tbot - Told);
	      Tnew = Told + dt;
	    }
	  }
	} else {
	  if (Htop > Htarget) {
	    if (Tnew > (0.75 * Ttop + 0.25 * Told)) {
	      dt = 0.75 * (Ttop - Told);
	      Tnew = Told + dt;
	    }
	  }
	}
      }
      // Check Max and Min values
      if (Tnew > Tmax) {
	if (!ignoreBounds) {
	  if (doUV) {
	    setTemperature(Tmax);
	    Hmax = intEnergy_mass();
	  } else {
	    setState_TP(Tmax, p);
	    Hmax = enthalpy_mass();
	  }
	  if (Hmax >= Htarget) {
	    if (Htop < Htarget) {
	      Ttop = Tmax;
	      Htop = Hmax;
	    }
	  } else {
	    Tnew = Tmax + 1.0;
	    ignoreBounds = true;
	  }
	}
      }
      if (Tnew < Tmin) {
	if (!ignoreBounds) {
	  if (doUV) {
	    setTemperature(Tmin);
	    Hmin = intEnergy_mass();
	  } else {
	    setState_TP(Tmin, p);
	    Hmin = enthalpy_mass();
	  }
	  if (Hmin <= Htarget) {
	    if (Hbot > Htarget) {
	      Tbot = Tmin;
	      Hbot = Hmin;
	    }
	  } else {
	    Tnew = Tmin - 1.0;
	    ignoreBounds = true;
	  }
	}
      }
 
      // Try to keep phase within its region of stability
      // -> Could do a lot better if I calculate the
      //    spinodal value of H.
      for (int its = 0; its < 10; its++) {
	Tnew = Told + dt;
	if (doUV) {
	  setTemperature(Tnew);
	  Hnew = intEnergy_mass();
	  Cpnew = cv_mass();
	} else {
	  setState_TP(Tnew, p);
	  Hnew = enthalpy_mass();
	  Cpnew = cp_mass();
	}
	if (Cpnew < 0.0) {
	  unstablePhaseNew = true;
	} else {
	  break;
	  unstablePhaseNew = false;
	}
	if (unstablePhase == false) {
	  if (unstablePhaseNew == true) {
	    dt *= 0.25;
	  }
	}
      }

      if (Hnew == Htarget) {
	return;
      } else if (Hnew > Htarget) {
	if ((Htop < Htarget) || (Hnew < Htop)) {
	  Htop = Hnew;
	  Ttop = Tnew;
	} 
      } else if (Hnew < Htarget) {
	if ((Hbot > Htarget) || (Hnew > Hbot)) {
	  Hbot = Hnew;
	  Tbot = Tnew;
	}
      }
      // Convergence in H
      double Herr = Htarget - Hnew;
      double acpd = MAX(fabs(cpd), 1.0E-5);
      double denom = MAX(fabs(Htarget), acpd * dTtol);
      double HConvErr = fabs((Herr)/denom);
      if (HConvErr < 0.00001 *dTtol) {
	return;
      }
      if (fabs(dt) < dTtol) {
	return;
      }
    }
    throw CanteraError("setState_HPorUV","No convergence. dt = " + fp2str(dt));
  }

  void ThermoPhase::setState_SP(doublereal Starget, doublereal p, 
				doublereal dTtol) {
    setState_SPorSV(Starget, p, dTtol, false);
  }

  void ThermoPhase::setState_SV(doublereal Starget, doublereal v, 
				doublereal dTtol) {
    setState_SPorSV(Starget, v, dTtol, true);
  }

  void ThermoPhase::setState_SPorSV(doublereal Starget, doublereal p, 
				    doublereal dTtol, bool doSV) {
    doublereal v = 0.0;
    doublereal dt;
    if (doSV) {
      v = p;
      setDensity(1.0/v); 
    } else {
      setPressure(p);
    }
    double Tmax = maxTemp() + 0.1;
    double Tmin = minTemp() - 0.1;

    // Make sure we are within the temperature bounds at the start
    // of the iteration
    double Tnew = temperature();
    if (Tnew > Tmax) {
      Tnew = Tmax - 1.0;
      if (doSV) {
	setTemperature(Tnew);
      } else {
	setState_TP(Tnew, p);
      }
    }
    if (Tnew < Tmin) {
      Tnew = Tmin + 1.0;
      if (doSV) {
	setTemperature(Tnew);
      } else {
	setState_TP(Tnew, p);
      }
    }

    double Snew = entropy_mass();
    double Cpnew = 0.0;
    if (doSV) {
      Cpnew = cv_mass();
    } else {
      Cpnew = cp_mass();
    }

    double Stop = Snew;
    double Ttop = Tnew;
    double Sbot = Snew;
    double Tbot = Tnew;
    double Told = Tnew;
    double Sold = Snew;

    bool ignoreBounds = false;
    // Unstable phases are those for which
    // cp < 0.0. These are possible for cases where
    // we have passed the spinodal curve.
    bool unstablePhase = false;
    bool unstablePhaseNew = false;
   

    // Newton iteration
    for (int n = 0; n < 500; n++) {
      Told = Tnew;
      Sold = Snew;
      double cpd = Cpnew;
      if (cpd < 0.0) {
	unstablePhase = true;
      }
      dt = (Starget - Sold)*Told/cpd;

      // limit step size to 200 K
      if (dt > 100.0)       dt =  100.0;
      else if (dt < -100.0) dt = -100.0; 
      Tnew = Told + dt;
      // Limit the step size so that we are convergent
      if (dt > 0.0) {
	if (!unstablePhase) {
	  if (Stop > Starget) {
	    if (Tnew > Ttop) {
	      dt = 0.75 * (Ttop - Told);
	    Tnew = Told + dt;
	    }
	  }
	} else {
	  if (Sbot < Starget) {
	    if (Tnew < Tbot) {
	      dt = 0.75 * (Tbot - Told);
	      Tnew = Told + dt;
	    }
	  }
	}
      } else {
	if (!unstablePhase) {
	  if (Sbot < Starget) {
	    if (Tnew < Tbot) {
	      dt = 0.75 * (Tbot - Told);
	      Tnew = Told + dt;
	    }
	  }
	} else {
	  if (Stop > Starget) {
	    if (Tnew > Ttop) {
	      dt = 0.75 * (Ttop - Told);
	      Tnew = Told + dt;
	    }
	  }
	}
      }
      // Check Max and Min values
      if (Tnew > Tmax) {
	if (!ignoreBounds) {
	  if (doSV) {
	    setTemperature(Tmax);
	  } else {
	    setState_TP(Tmax, p);
	  }
	  double Smax = entropy_mass();
	  if (Smax >= Starget) {
	    if (Stop < Starget) {
	      Ttop = Tmax;
	      Stop = Smax;
	    }
	  } else {
	    Tnew = Tmax + 1.0;
	    ignoreBounds = true;
	  }
	}
      }
      if (Tnew < Tmin) {
	if (!ignoreBounds) {
	  if (doSV) {
	    setTemperature(Tmin);
	  } else {
	    setState_TP(Tmin, p);
	  }
	  double Smin = enthalpy_mass();
	  if (Smin <= Starget) {
	    if (Sbot > Starget) {
	      Sbot = Tmin;
	      Sbot = Smin;
	    }
	  } else {
	    Tnew = Tmin - 1.0;
	    ignoreBounds = true;
	  }
	}
      }
 
      // Try to keep phase within its region of stability
      // -> Could do a lot better if I calculate the
      //    spinodal value of H.
      for (int its = 0; its < 10; its++) {
	Tnew = Told + dt;
	if (doSV) {
	  setTemperature(Tnew);
	  Cpnew = cv_mass();
	} else {
	  setState_TP(Tnew, p);
	  Cpnew = cp_mass();
	}
	Snew = entropy_mass();
	if (Cpnew < 0.0) {
	  unstablePhaseNew = true;
	} else {
	  break;
	  unstablePhaseNew = false;
	}
	if (unstablePhase == false) {
	  if (unstablePhaseNew == true) {
	    dt *= 0.25;
	  }
	}
      }

      if (Snew == Starget) {
	return;
      } else if (Snew > Starget) {
	if ((Stop < Starget) || (Snew < Stop)) {
	  Stop = Snew;
	  Ttop = Tnew;
	} 
      } else if (Snew < Starget) {
	if ((Sbot > Starget) || (Snew > Sbot)) {
	  Sbot = Snew;
	  Tbot = Tnew;
	}
      }
      // Convergence in S
      double Serr = Starget - Snew;
      double acpd = MAX(fabs(cpd), 1.0E-5);
      double denom = MAX(fabs(Starget), acpd * dTtol);
      double SConvErr = fabs((Serr * Tnew)/denom);
      if (SConvErr < 0.00001 *dTtol) {
	return;
      }
      if (fabs(dt) < dTtol) {
	return;
      }
    }
    throw CanteraError("setState_SPorSV","No convergence. dt = " + fp2str(dt));
  }

  doublereal ThermoPhase::err(std::string msg) const {
    throw CanteraError("ThermoPhase","Base class method "
		       +msg+" called. Equation of state type: "+int2str(eosType()));
    return 0;
  }

  /*
   * Returns the units of the standard and general concentrations
   * Note they have the same units, as their divisor is 
   * defined to be equal to the activity of the kth species
   * in the solution, which is unitless.
   *
   * This routine is used in print out applications where the
   * units are needed. Usually, MKS units are assumed throughout
   * the program and in the XML input files. 
   *
   * On return uA contains the powers of the units (MKS assumed)
   * of the standard concentrations and generalized concentrations
   * for the kth species.
   *
   * The base %ThermoPhase class assigns thedefault quantities
   * of (kmol/m3).
   * Inherited classes are responsible for overriding the default 
   * values if necessary.
   *
   *  uA[0] = kmol units - default  = 1
   *  uA[1] = m    units - default  = -nDim(), the number of spatial
   *                                dimensions in the Phase class.
   *  uA[2] = kg   units - default  = 0;
   *  uA[3] = Pa(pressure) units - default = 0;
   *  uA[4] = Temperature units - default = 0;
   *  uA[5] = time units - default = 0
   */
  void ThermoPhase::getUnitsStandardConc(double *uA, int k, int sizeUA) const {
    for (int i = 0; i < sizeUA; i++) {
      if (i == 0) uA[0] = 1.0;
      if (i == 1) uA[1] = -nDim();
      if (i == 2) uA[2] = 0.0;
	  if (i == 3) uA[3] = 0.0;
	  if (i == 4) uA[4] = 0.0;
	  if (i == 5) uA[5] = 0.0;
    }
  }

    /*
     * initThermoFile():
     *
     * Initialization of a phase using an xml file.
     *
     * This routine is a precursor to initThermoXML(XML_Node*)
     * routine, which does most of the work. 
     *
     * @param infile XML file containing the description of the
     *        phase
     *
     * @param id  Optional parameter identifying the name of the
     *            phase. If none is given, the first XML
     *            phase element will be used.
     */
    void ThermoPhase::initThermoFile(std::string inputFile, std::string id) {

	if (inputFile.size() == 0) {
	  throw CanteraError("ThermoPhase::initThermoFile",
			     "input file is null");
	}
	string path = findInputFile(inputFile);
	ifstream fin(path.c_str());
	if (!fin) {
	  throw CanteraError("initThermoFile","could not open "
			     +path+" for reading.");
	}
	/*
	 * The phase object automatically constructs an XML object.
	 * Use this object to store information.
	 */
	XML_Node &phaseNode_XML = xml();
        XML_Node *fxml = new XML_Node();
	fxml->build(fin);
	XML_Node *fxml_phase = findXMLPhase(fxml, id);
	if (!fxml_phase) {
	  throw CanteraError("ThermoPhase::initThermo",
			     "ERROR: Can not find phase named " +
			     id + " in file named " + inputFile);
	}
	fxml_phase->copy(&phaseNode_XML);	
	initThermoXML(*fxml_phase, id);
	delete fxml;
    }

    /*
     *   Import and initialize a ThermoPhase object
     *
     *   This function is called from importPhase() 
     *   after the elements and the
     *   species are initialized with default ideal solution
     *   level data.
     *
     * @param phaseNode This object must be the phase node of a
     *             complete XML tree
     *             description of the phase, including all of the
     *             species data. In other words while "phase" must
     *             point to an XML phase object, it must have
     *             sibling nodes "speciesData" that describe
     *             the species in the phase.
     * @param id   ID of the phase. If nonnull, a check is done
     *             to see if phaseNode is pointing to the phase
     *             with the correct id. 
     */
    void ThermoPhase::initThermoXML(XML_Node& phaseNode, std::string id) {
      /*
       * The default implementation just calls initThermo(), which
       * inheriting classes may override.
       */
      initThermo();
      /*
       * and sets the state
       */
      if (phaseNode.hasChild("state")) {
	XML_Node& stateNode = phaseNode.child("state");
	setStateFromXML(stateNode);
      }
    }
  
    /*
     * Initialize. 
     *
     * This method is provided to allow
     * subclasses to perform any initialization required after all
     * species have been added. For example, it might be used to
     * resize internal work arrays that must have an entry for
     * each species.  The base class implementation does nothing,
     * and subclasses that do not require initialization do not
     * need to overload this method.  When importing a CTML phase
     * description, this method is called just prior to returning
     * from function importPhase.
     *
     * @see importCTML.cpp
     */
    void ThermoPhase::initThermo() {
      // Check to see that there is at least one species defined in the phase
      if (m_kk <= 0) {
        throw CanteraError("ThermoPhase::initThermo()",
                           "Number of species is less than or equal to zero");
      }
    }
    
  /**
   * Set the thermodynamic state.
   */
    void ThermoPhase::setStateFromXML(const XML_Node& state) {

      string comp = getString(state,"moleFractions");
      if (comp != "") 
	setMoleFractionsByName(comp);
      else {
	comp = getString(state,"massFractions");
	if (comp != "") 
	  setMassFractionsByName(comp);
      }
      if (state.hasChild("temperature")) {
	double t = getFloat(state, "temperature", "temperature");
	setTemperature(t);
      }
      if (state.hasChild("pressure")) {
	double p = getFloat(state, "pressure", "pressure");
	setPressure(p);
      }
      if (state.hasChild("density")) {
	double rho = getFloat(state, "density", "density");
	setDensity(rho);
      }
    }



  /*
   * Called by function 'equilibrate' in ChemEquil.h to transfer
   * the element potentials to this object after every successful
   *  equilibration routine.
   * The element potentials are storred in their dimensionless
   * forms, calculated by dividing by RT.
   *    @param lambda vector containing the element potentials.
   *           Length = nElements. Units are Joules/kmol.
   */
  void ThermoPhase::setElementPotentials(const vector_fp& lambda) {
    doublereal rrt = 1.0/(GasConstant* temperature());
    int mm = nElements();
    if (lambda.size() < (size_t) mm) {
      throw CanteraError("setElementPotentials", "lambda too small");
    }
    if (!m_hasElementPotentials) {
      m_lambdaRRT.resize(mm);
    }
    for (int m = 0; m < mm; m++) {
      m_lambdaRRT[m] = lambda[m] * rrt;
    }
    m_hasElementPotentials = true;
  }

  /*
   * Returns the storred element potentials.
   * The element potentials are retrieved from their storred
   * dimensionless forms by multiplying by RT.
   * @param lambda Vector containing the element potentials.
   *        Length = nElements. Units are Joules/kmol.
   */
   bool ThermoPhase::getElementPotentials(doublereal* lambda) const {
    doublereal rt = GasConstant* temperature();
    int mm = nElements();
    if (m_hasElementPotentials) {
      for (int m = 0; m < mm; m++) {
	lambda[m] =  m_lambdaRRT[m] * rt;
      }
    }
    return (m_hasElementPotentials);
  }

}
