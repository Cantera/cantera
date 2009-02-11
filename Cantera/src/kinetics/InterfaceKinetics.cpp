/**
 *  @file InterfaceKinetics.cpp 
 *
 */

// Copyright 2002  California Institute of Technology


// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "InterfaceKinetics.h"
#include "EdgeKinetics.h"
#include "SurfPhase.h"

#include "ReactionData.h"
#include "RateCoeffMgr.h"

#include "ImplicitSurfChem.h"

using namespace std;

namespace Cantera {

    //////////////////////////////////////////////////////////////////

    /**
     * Construct an empty InterfaceKinetics reaction mechanism.
     *  @param thermo This is an optional parameter that may be
     *         used to initialize the inherited Kinetics class with
     *         one ThermoPhase class object -> in other words it's
     *         useful for initialization of homogeneous kinetics 
     *         mechanisms.
     */    
    InterfaceKinetics::
    InterfaceKinetics(thermo_t* thermo) :
        Kinetics(),
        m_kk(0), 
        m_redo_rates(false),
        m_nirrev(0), 
        m_nrev(0),
        m_surf(0),
        m_integrator(0),
        m_finalized(false),
        m_has_coverage_dependence(false),
        m_has_electrochem_rxns(false),
	m_ioFlag(0)
    {
        if (thermo != 0) addPhase(*thermo);
        m_kdata = new InterfaceKineticsData;
        m_kdata->m_temp = 0.0;
    }

    /**
     * Destructor
     */
    InterfaceKinetics::
    ~InterfaceKinetics(){
        delete m_kdata;
	if (m_integrator) {
	  delete m_integrator; 
	}
    }

  //  Copy Constructor for the %InterfaceKinetics object.
  /*
   * Currently, this is not fully implemented. If called it will
   * throw an exception.
   */
  InterfaceKinetics::InterfaceKinetics(const InterfaceKinetics &right) :
    Kinetics(),
    m_kk(0), 
    m_redo_rates(false),
    m_nirrev(0), 
    m_nrev(0),
    m_surf(0),
    m_integrator(0),
    m_finalized(false),
    m_has_coverage_dependence(false),
    m_has_electrochem_rxns(false),
    m_ioFlag(0)
  {  
    m_kdata = new InterfaceKineticsData;
    m_kdata->m_temp = 0.0;
    /*
     * Call the assignment operator
     */
    *this = operator=(right);
  }
  
  // Assignment operator
  /*
   *  This is NOT a virtual function.
   *
   * @param right    Reference to %Kinetics object to be copied into the
   *                 current one.
   */
  InterfaceKinetics& InterfaceKinetics::
  operator=(const InterfaceKinetics &right) {
    /*
     * Check for self assignment.
     */
    if (this == &right) return *this;

    Kinetics::operator=(right);

    m_kk                   = right.m_kk;
    m_revindex             = right.m_revindex;
    m_rates                = right.m_rates;
    m_redo_rates           = right.m_redo_rates;
    m_index                = right.m_index;
    m_irrev                = right.m_irrev; 
    m_rxnstoich            = right.m_rxnstoich;
    m_nirrev               = right.m_nirrev;
    m_nrev                 = right.m_nrev;
    m_rrxn                 = right.m_rrxn;
    m_prxn                 = right.m_prxn;
    m_rxneqn               = right.m_rxneqn;
    *m_kdata               = *right.m_kdata;  // needs to be developed
    m_mu0                  = right.m_mu0;
    m_phi                  = right.m_phi;
    m_pot                  = right.m_pot;
    m_rwork                = right.m_rwork;
    m_E                    = right.m_E;
    m_surf                 = right.m_surf;  //DANGER - shallow copy
    m_integrator           = right.m_integrator;  //DANGER - shallow copy
    m_beta                 = right.m_beta;
    m_ctrxn                = right.m_ctrxn;
    m_finalized            = right.m_finalized;
    m_has_coverage_dependence = right.m_has_coverage_dependence;
    m_has_electrochem_rxns = right.m_has_electrochem_rxns;
    m_ioFlag               = right.m_ioFlag;

    return *this;
  }


  // Duplication routine for objects which inherit from
  // Kinetics
  /*
   *  This virtual routine can be used to duplicate %Kinetics objects
   *  inherited from %Kinetics even if the application only has
   *  a pointer to %Kinetics to work with.
   *
   *  These routines are basically wrappers around the derived copy
   *  constructor.
   */
  Kinetics *InterfaceKinetics::duplMyselfAsKinetics() const {
    InterfaceKinetics* tp = new InterfaceKinetics(*this);
    return dynamic_cast<Kinetics *>(tp);
  }


    /**
     * Update properties that depend on temperature
     *
     */ 
    void InterfaceKinetics::
    _update_rates_T() {
        _update_rates_phi();
        if (m_has_coverage_dependence) {
            m_surf->getCoverages(DATA_PTR(m_conc));
            m_rates.update_C(DATA_PTR(m_conc));
            m_redo_rates = true;
        }
        doublereal T = thermo(surfacePhaseIndex()).temperature();
        if (T != m_kdata->m_temp || m_redo_rates) {
            m_kdata->m_logtemp = log(T);
            m_rates.update(T, m_kdata->m_logtemp, DATA_PTR(m_kdata->m_rfn));
            if (m_has_electrochem_rxns)
                applyButlerVolmerCorrection(DATA_PTR(m_kdata->m_rfn));
            m_kdata->m_temp = T;
            updateKc();
            m_kdata->m_ROP_ok = false;
            m_redo_rates = false;
        }
    }

    void InterfaceKinetics::
    _update_rates_phi() {
        int np = nPhases();
        for (int n = 0; n < np; n++) {
            if (thermo(n).electricPotential() != m_phi[n]) {
                m_phi[n] = thermo(n).electricPotential();
                m_redo_rates = true;
            }
        }
    }


  /**
   * Update properties that depend on concentrations. This method
   * fills out the array of generalized concentrations by calling
   * method getActivityConcentrations for each phase, which classes
   * representing phases should overload to return the appropriate
   * quantities.
   */ 
  void InterfaceKinetics::
  _update_rates_C() {
    int n;

    int np = nPhases();
    for (n = 0; n < np; n++) {
      /*
       * We call the getActivityConcentrations function of each
       * ThermoPhase class that makes up this kinetics object to 
       * obtain the generalized concentrations for species within that 
       * class. This is collected in the vector m_conc. m_start[]
       * are integer indecises for that vector denoting the start of the
       * species for each phase.
       */
      thermo(n).getActivityConcentrations(DATA_PTR(m_conc) + m_start[n]);
    }
    m_kdata->m_ROP_ok = false;
  }


  /**
   * Update the equilibrium constants in molar units for all
   * reversible reactions. Irreversible reactions have their 
   * equilibrium constant set to zero.
   */
  void InterfaceKinetics::updateKc() {
    int i, irxn;

    vector_fp& m_rkc = m_kdata->m_rkcn;
    fill(m_rkc.begin(), m_rkc.end(), 0.0);

    //static vector_fp mu(nTotalSpecies());
    if (m_nrev > 0) {

      int n, nsp, k, ik = 0;
      doublereal rt = GasConstant*thermo(0).temperature();
      doublereal rrt = 1.0 / rt;
      int np = nPhases();
      for (n = 0; n < np; n++) {
	thermo(n).getStandardChemPotentials(DATA_PTR(m_mu0) + m_start[n]);
	nsp = thermo(n).nSpecies();
	for (k = 0; k < nsp; k++) {
	  m_mu0[ik] -= rt * thermo(n).logStandardConc(k);
	  m_mu0[ik] += Faraday * m_phi[n] * thermo(n).charge(k);
	  ik++;
	}
      }

      // compute Delta mu^0 for all reversible reactions
      m_rxnstoich.getRevReactionDelta(m_ii, DATA_PTR(m_mu0), 
				      DATA_PTR(m_rkc));

      for (i = 0; i < m_nrev; i++) {
	irxn = m_revindex[i];
	if (irxn < 0 || irxn >= nReactions()) {
	  throw CanteraError("InterfaceKinetics",
			     "illegal value: irxn = "+int2str(irxn));
	}
	m_rkc[irxn] = exp(m_rkc[irxn]*rrt);
      }
      for (i = 0; i != m_nirrev; ++i) {
	m_rkc[ m_irrev[i] ] = 0.0;
      }
    }
  }


    void InterfaceKinetics::checkPartialEquil() {
        int i, irxn;
        vector_fp dmu(nTotalSpecies(), 0.0);
        vector_fp rmu(nReactions(), 0.0);
        vector_fp frop(nReactions(), 0.0);
        vector_fp rrop(nReactions(), 0.0);
        vector_fp netrop(nReactions(), 0.0);
        if (m_nrev > 0) {
            doublereal rt = GasConstant*thermo(0).temperature();
            cout << "T = " << thermo(0).temperature() << " " << rt << endl;
            int n, nsp, k, ik=0;
            //doublereal rt = GasConstant*thermo(0).temperature();
            //            doublereal rrt = 1.0/rt;
            int np = nPhases();
            doublereal delta;
            for (n = 0; n < np; n++) {
                thermo(n).getChemPotentials(DATA_PTR(dmu) + m_start[n]);
                nsp = thermo(n).nSpecies();
                for (k = 0; k < nsp; k++) {
                    delta = Faraday * m_phi[n] * thermo(n).charge(k);
                    //cout << thermo(n).speciesName(k) << "   " << (delta+dmu[ik])/rt << " " << dmu[ik]/rt << endl;
                    dmu[ik] += delta;
                    ik++;
                }
            }

            // compute Delta mu^ for all reversible reactions
            m_rxnstoich.getRevReactionDelta(m_ii, DATA_PTR(dmu), DATA_PTR(rmu));
            getFwdRatesOfProgress(DATA_PTR(frop));
            getRevRatesOfProgress(DATA_PTR(rrop));
            getNetRatesOfProgress(DATA_PTR(netrop));
            for (i = 0; i < m_nrev; i++) {
                irxn = m_revindex[i];
                cout << "Reaction " << reactionString(irxn) 
                     << "  " << rmu[irxn]/rt << endl;
                printf("%12.6e  %12.6e  %12.6e  %12.6e \n", 
                    frop[irxn], rrop[irxn], netrop[irxn], 
                    netrop[irxn]/(frop[irxn] + rrop[irxn]));
            }
        }
    }


    /**
     * Get the equilibrium constants of all reactions, whether
     * reversible or not.
     */
    void InterfaceKinetics::getEquilibriumConstants(doublereal* kc) {
        int i;

        int n, nsp, k, ik=0;
        doublereal rt = GasConstant*thermo(0).temperature();
        doublereal rrt = 1.0/rt;
        int np = nPhases();
        for (n = 0; n < np; n++) {
            thermo(n).getStandardChemPotentials(DATA_PTR(m_mu0) + m_start[n]);
            nsp = thermo(n).nSpecies();
            for (k = 0; k < nsp; k++) {
                m_mu0[ik] -= rt*thermo(n).logStandardConc(k);
                m_mu0[ik] += Faraday * m_phi[n] * thermo(n).charge(k);
                ik++;
            }
        }

        fill(kc, kc + m_ii, 0.0);

        m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(m_mu0), kc);

        for (i = 0; i < m_ii; i++) {
            kc[i] = exp(-kc[i]*rrt);
        }
    }

  // Returns the Species creation rates [kmol/m^2/s]. 
  /*
   *   Return the species
   * creation rates in array cdot, which must be
   * dimensioned at least as large as the total number of
   * species in all phases of the kinetics
   * model
   *  
   *  @param cdot Vector containing the creation rates.
   *              length = m_kk. units = kmol/m^2/s
   */
  void InterfaceKinetics::getCreationRates(doublereal* cdot) {
    updateROP();
    m_rxnstoich.getCreationRates(m_kk, &m_kdata->m_ropf[0],
				 &m_kdata->m_ropr[0], cdot);
  }
   
  //  Return the Species destruction rates [kmol/m^2/s].
  /*
   *  Return the species destruction rates in array ddot, which must be
   *  dimensioned at least as large as the total number of
   *  species in all phases of the kinetics model
   */
  void InterfaceKinetics::getDestructionRates(doublereal* ddot) {
    updateROP();
    m_rxnstoich.getDestructionRates(m_kk, &m_kdata->m_ropf[0],
				    &m_kdata->m_ropr[0], ddot);
  }

  // Return the species net production rates
  /*
   * Species net production rates [kmol/m^2/s]. Return the species
   * net production rates (creation - destruction) in array
   * wdot, which must be dimensioned at least as large as the
   * total number of species in all phases of the kinetics
   * model
   *
   * @param net  Vector of species production rates.
   *             units kmol m-d s-1, where d is dimension.
   */
  void InterfaceKinetics::getNetProductionRates(doublereal* net) {
    updateROP();
    m_rxnstoich.getNetProductionRates(m_kk,
				      &m_kdata->m_ropnet[0],
				      net);
  }

  /**
   * For reactions that transfer charge across a potential difference,
   * the activation energies are modified by the potential difference.
   * (see, for example, ...). This method applies this correction.
   */
  void InterfaceKinetics::applyButlerVolmerCorrection(doublereal* kf) {
    int i;

    int n, nsp, k, ik=0;
    doublereal rt = GasConstant*thermo(0).temperature();
    doublereal rrt = 1.0/rt;
    int np = nPhases();

    // compute the electrical potential energy of each species
    for (n = 0; n < np; n++) {
      nsp = thermo(n).nSpecies();
      for (k = 0; k < nsp; k++) {
	m_pot[ik] = Faraday*thermo(n).charge(k)*m_phi[n];
	ik++;
      }
    }

    // Compute the change in electrical potential energy for each
    // reaction. This will only be non-zero if a potential
    // difference is present.
    m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(m_pot), 
				 DATA_PTR(m_rwork));

    // Modify the reaction rates. Only modify those with a
    // non-zero activation energy, and do not decrease the
    // activation energy below zero.
    doublereal eamod;
#ifdef DEBUG_MODE
    double ea;
#endif
    int nct = m_beta.size();
    int irxn;
    for (i = 0; i < nct; i++) {
      irxn = m_ctrxn[i];
      eamod = m_beta[i]*m_rwork[irxn];
      if (eamod != 0.0 && m_E[irxn] != 0.0) {
#ifdef DEBUG_MODE
	ea = GasConstant * m_E[irxn];
	if (eamod + ea < 0.0) {
	  writelog("Warning: act energy mod too large!\n");
	  writelog("  Delta phi = "+fp2str(m_rwork[irxn]/Faraday)+"\n");
	  writelog("  Delta Ea = "+fp2str(eamod)+"\n");
	  writelog("  Ea = "+fp2str(ea)+"\n");
	  for (n = 0; n < np; n++) {
	    writelog("Phase "+int2str(n)+": phi = "
		     +fp2str(m_phi[n])+"\n");
	  }
	}
#endif
	kf[irxn] *= exp(-eamod*rrt);
      }
    }
  }




    /**
     * Update the rates of progress of the reactions in the reaciton
     * mechanism. This routine operates on internal data.
     */
    void InterfaceKinetics::getFwdRateConstants(doublereal* kfwd) {

        _update_rates_T();
        _update_rates_C();

        const vector_fp& rf = m_kdata->m_rfn;

        // copy rate coefficients into kfwd
        copy(rf.begin(), rf.end(), kfwd);

        // multiply by perturbation factor
        multiply_each(kfwd, kfwd + nReactions(), m_perturb.begin());
           
    }


    /**
     * Update the rates of progress of the reactions in the reaciton
     * mechanism. This routine operates on internal data.
     */
    void InterfaceKinetics::getRevRateConstants(doublereal* krev, bool doIrreversible) {
        getFwdRateConstants(krev);
        if (doIrreversible) {
            doublereal *tmpKc = DATA_PTR(m_kdata->m_ropnet);
            getEquilibriumConstants(tmpKc);
            for (int i = 0; i < m_ii; i++) {
                krev[i] /=  tmpKc[i];
            }
        }
        else {
            const vector_fp& rkc = m_kdata->m_rkcn;
            multiply_each(krev, krev + nReactions(), rkc.begin());
        }
    }


    void InterfaceKinetics::getActivationEnergies(doublereal *E) {
        copy(m_E.begin(), m_E.end(), E);
    }

    /**
     * Update the rates of progress of the reactions in the reaciton
     * mechanism. This routine operates on internal data.
     */
    void InterfaceKinetics::updateROP() {

        _update_rates_T();
        _update_rates_C();

        if (m_kdata->m_ROP_ok) return;

        const vector_fp& rf = m_kdata->m_rfn;
        const vector_fp& m_rkc = m_kdata->m_rkcn;
        array_fp& ropf = m_kdata->m_ropf;
        array_fp& ropr = m_kdata->m_ropr;
        array_fp& ropnet = m_kdata->m_ropnet;

        // copy rate coefficients into ropf
        copy(rf.begin(), rf.end(), ropf.begin());

        // multiply by perturbation factor
        multiply_each(ropf.begin(), ropf.end(), m_perturb.begin());
           
        // copy the forward rates to the reverse rates                
        copy(ropf.begin(), ropf.end(), ropr.begin());
        
        // for reverse rates computed from thermochemistry, multiply
        // the forward rates copied into m_ropr by the reciprocals of
        // the equilibrium constants
        multiply_each(ropr.begin(), ropr.end(), m_rkc.begin());

        // multiply ropf by concentration products
        m_rxnstoich.multiplyReactants(DATA_PTR(m_conc), DATA_PTR(ropf)); 
        //m_reactantStoich.multiply(m_conc.begin(), ropf.begin()); 

        // for reversible reactions, multiply ropr by concentration
        // products
        m_rxnstoich.multiplyRevProducts(DATA_PTR(m_conc), 
            DATA_PTR(ropr)); 
        //m_revProductStoich.multiply(m_conc.begin(), ropr.begin());

        // do global reactions
        //m_globalReactantStoich.power(m_conc.begin(), ropf.begin());

        for (int j = 0; j != m_ii; ++j) {
            ropnet[j] = ropf[j] - ropr[j];
        }        

        m_kdata->m_ROP_ok = true;
    }


    /**
     *
     * getDeltaGibbs():
     *
     * Return the vector of values for the reaction gibbs free energy
     * change
     * These values depend upon the concentration
     * of the ideal gas.
     *
     *  units = J kmol-1
     */
    void InterfaceKinetics::getDeltaGibbs(doublereal* deltaG) {
	/*
	 * Get the chemical potentials of the species in the 
	 * ideal gas solution.
	 */
        int np = nPhases();
        int n;
        for (n = 0; n < np; n++) {
            thermo(n).getChemPotentials(DATA_PTR(m_grt) + m_start[n]);
        }
        //for (n = 0; n < m_grt.size(); n++) {
        //    cout << n << "G_RT = " << m_grt[n] << endl;
        //}

	/*
	 * Use the stoichiometric manager to find deltaG for each
	 * reaction.
	 */
	m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(m_grt), deltaG);
    }
    
    /**
     *
     * getDeltaEnthalpy():
     * 
     * Return the vector of values for the reactions change in
     * enthalpy.
     * These values depend upon the concentration
     * of the solution.
     *
     *  units = J kmol-1
     */
    void InterfaceKinetics::getDeltaEnthalpy(doublereal* deltaH) {
	/*
	 * Get the partial molar enthalpy of all species in the 
	 * ideal gas.
	 */
        int np = nPhases();
        int n;
        for (n = 0; n < np; n++) {
            thermo(n).getPartialMolarEnthalpies(DATA_PTR(m_grt) + m_start[n]);
        }
	/*
	 * Use the stoichiometric manager to find deltaG for each
	 * reaction.
	 */
	m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(m_grt), deltaH);
    }

   
  // Return the vector of values for the change in
  // entropy due to each reaction
  /*
   * These values depend upon the concentration
   * of the solution.
   *
   *  units = J kmol-1 Kelvin-1
   *
   * @param deltaS vector of Enthalpy changes 
   *        Length = m_ii, number of reactions
   *         
   */
  void InterfaceKinetics::getDeltaEntropy(doublereal* deltaS) {
    /*
     * Get the partial molar entropy of all species in all of
     * the phases
     */
    int np = nPhases();
    int n;
    for (n = 0; n < np; n++) {
      thermo(n).getPartialMolarEntropies(DATA_PTR(m_grt) + m_start[n]);
    }
    /*
     * Use the stoichiometric manager to find deltaS for each
     * reaction.
     */
    m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(m_grt), deltaS);
  }

    /**
     *
     * getDeltaSSGibbs():
     *
     * Return the vector of values for the reaction 
     * standard state gibbs free energy change.
     * These values don't depend upon the concentration
     * of the solution.
     *
     *  units = J kmol-1
     */
    void InterfaceKinetics::getDeltaSSGibbs(doublereal* deltaG) {
	/*
	 *  Get the standard state chemical potentials of the species.
         *  This is the array of chemical potentials at unit activity 
	 *  We define these here as the chemical potentials of the pure
	 *  species at the temperature and pressure of the solution.
	 */
        int np = nPhases();
        int n;
        for (n = 0; n < np; n++) {
            thermo(n).getStandardChemPotentials(DATA_PTR(m_grt) + m_start[n]);
        }
	/*
	 * Use the stoichiometric manager to find deltaG for each
	 * reaction.
	 */
	m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(m_grt), deltaG);
    }

    /**
     *
     * getDeltaSSEnthalpy():
     *
     * Return the vector of values for the change in the
     * standard state enthalpies of reaction.
     * These values don't depend upon the concentration
     * of the solution.
     *
     *  units = J kmol-1
     */
    void InterfaceKinetics::getDeltaSSEnthalpy(doublereal* deltaH) {
	/*
	 *  Get the standard state enthalpies of the species.
         *  This is the array of chemical potentials at unit activity 
	 *  We define these here as the enthalpies of the pure
	 *  species at the temperature and pressure of the solution.
	 */
        int np = nPhases();
        int n;
        for (n = 0; n < np; n++) {
            thermo(n).getEnthalpy_RT(DATA_PTR(m_grt) + m_start[n]);
        }
	doublereal RT = thermo().temperature() * GasConstant;
	for (int k = 0; k < m_kk; k++) {
	  m_grt[k] *= RT;
	}
	/*
	 * Use the stoichiometric manager to find deltaG for each
	 * reaction.
	 */
	m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(m_grt), deltaH);
    }

    /*********************************************************************
     *
     * getDeltaSSEntropy():
     *
     * Return the vector of values for the change in the
     * standard state entropies for each reaction.
     * These values don't depend upon the concentration
     * of the solution.
     *
     *  units = J kmol-1 Kelvin-1
     */
    void InterfaceKinetics::getDeltaSSEntropy(doublereal* deltaS) {
	/*
	 *  Get the standard state entropy of the species.
	 *  We define these here as the entropies of the pure
	 *  species at the temperature and pressure of the solution.
	 */
        int np = nPhases();
        int n;
        for (n = 0; n < np; n++) {
            thermo(n).getEntropy_R(DATA_PTR(m_grt) + m_start[n]);
        }
	doublereal R = GasConstant;
	for (int k = 0; k < m_kk; k++) {
	  m_grt[k] *= R;
	}
	/*
	 * Use the stoichiometric manager to find deltaS for each
	 * reaction.
	 */
	m_rxnstoich.getReactionDelta(m_ii, DATA_PTR(m_grt), deltaS);
    }


    /**
     * Add a single reaction to the mechanism. This routine
     * must be called after init() and before finalize().
     * This function branches on the types of reactions allowed
     * by the interfaceKinetics manager in order to install
     * the reaction correctly in the manager.
     * The manager allows the following reaction types
     *  Elementary
     *  Surface
     *  Global  
     * There is no difference between elementary and surface 
     * reactions.
     */
    void InterfaceKinetics::
    addReaction(const ReactionData& r) {

        addElementaryReaction(r);

        // operations common to all reaction types
        installReagents( r );
        //installGroups(reactionNumber(), r.rgroups, r.pgroups);
        incrementRxnCount();
        m_rxneqn.push_back(r.equation);
    }


    void InterfaceKinetics::
    addElementaryReaction(const ReactionData& r) {
        int iloc;

        // install rate coeff calculator

        vector_fp rp = r.rateCoeffParameters;
        int ncov = r.cov.size();
        if (ncov > 3) {
            m_has_coverage_dependence = true;
        }
        for (int m = 0; m < ncov; m++) rp.push_back(r.cov[m]);

        iloc = m_rates.install( reactionNumber(),
				r.rateCoeffType, rp.size(), 
				DATA_PTR(rp) );

        // store activation energy
        m_E.push_back(r.rateCoeffParameters[2]);

        if (r.beta > 0.0) {
            m_has_electrochem_rxns = true;
            m_beta.push_back(r.beta);
            m_ctrxn.push_back(reactionNumber());
        }

        // add constant term to rate coeff value vector
        m_kdata->m_rfn.push_back(r.rateCoeffParameters[0]);                
        registerReaction( reactionNumber(), ELEMENTARY_RXN, iloc);
    }


  void InterfaceKinetics::setIOFlag(int ioFlag) {
    m_ioFlag = ioFlag;
    if (m_integrator) {
      m_integrator->setIOFlag(ioFlag);
    }
  }

//     void InterfaceKinetics::
//     addGlobalReaction(const ReactionData& r) {
            
//         int iloc;
//         // install rate coeff calculator
//         vector_fp rp = r.rateCoeffParameters;
//         int ncov = r.cov.size();
//         for (int m = 0; m < ncov; m++) rp.push_back(r.cov[m]);
//         iloc = m_rates.install( reactionNumber(),
//             r.rateCoeffType, rp.size(),
//             rp.begin() );
//         // store activation energy
//         m_E.push_back(r.rateCoeffParameters[2]);
//         // add constant term to rate coeff value vector
//         m_kdata->m_rfn.push_back(r.rateCoeffParameters[0]);

//         int nr = r.order.size();
//         vector_fp ordr(nr);
//         for (int n = 0; n < nr; n++) {
//             ordr[n] = r.order[n] - r.rstoich[n];
//         }
//         m_globalReactantStoich.add( reactionNumber(),
//             r.reactants, ordr);

//         registerReaction( reactionNumber(), GLOBAL_RXN, iloc);
//     }

        
  void InterfaceKinetics::installReagents(const ReactionData& r) {

    int n, ns, m; 
    doublereal nsFlt;
    /*
     * extend temporary storage by one for this rxn.
     */
    m_kdata->m_ropf.push_back(0.0);
    m_kdata->m_ropr.push_back(0.0);
    m_kdata->m_ropnet.push_back(0.0);
    m_kdata->m_rkcn.push_back(0.0);

    /*
     * Obtain the current reaction index for the reaction that we
     * are adding. The first reaction is labeled 0.
     */
    int rnum = reactionNumber();

    // vectors rk and pk are lists of species numbers, with
    // repeated entries for species with stoichiometric
    // coefficients > 1. This allows the reaction to be defined
    // with unity reaction order for each reactant, and so the
    // faster method 'multiply' can be used to compute the rate of
    // progress instead of 'power'.

    vector_int rk;
    int nr = r.reactants.size();
    for (n = 0; n < nr; n++) {
      nsFlt = r.rstoich[n];
      ns = (int) nsFlt;
      if ((doublereal) ns != nsFlt) {
	if (ns < 1) ns = 1;
      }
      /*
       * Add to m_rrxn. m_rrxn is a vector of maps. m_rrxn has a length
       * equal to the total number of species for each species, there
       * exists a map, with the reaction number being the key, and the
       * reactant stoichiometric coefficient being the value.
       */
      m_rrxn[r.reactants[n]][rnum] = nsFlt;
      for (m = 0; m < ns; m++) {
	rk.push_back(r.reactants[n]);
      }
    }
    /*
     * Now that we have rk[], we add it into the vector<vector_int> m_reactants
     * in the rnum index spot. Thus m_reactants[rnum] yields a vector
     * of reactants for the rnum'th reaction
     */
    m_reactants.push_back(rk);
    vector_int pk;
    int np = r.products.size();
    for (n = 0; n < np; n++) {
      nsFlt = r.pstoich[n];
      ns = (int) nsFlt;
      if ((doublereal) ns != nsFlt) {
	if (ns < 1) ns = 1;
      }
      /*
       * Add to m_prxn. m_prxn is a vector of maps. m_prxn has a length
       * equal to the total number of species for each species, there
       * exists a map, with the reaction number being the key, and the
       * product stoichiometric coefficient being the value.
       */
      m_prxn[r.products[n]][rnum] = nsFlt;
      for (m = 0; m < ns; m++) {
	pk.push_back(r.products[n]);
      }
    }
    /*
     * Now that we have pk[], we add it into the vector<vector_int> m_products
     * in the rnum index spot. Thus m_products[rnum] yields a vector
     * of products for the rnum'th reaction
     */
    m_products.push_back(pk);
    /*
     * Add this reaction to the stoichiometric coefficient manager. This
     * calculates rates of species production from reaction rates of 
     * progress.
     */
    m_rxnstoich.add( reactionNumber(), r);
    /*
     * register reaction in lists of reversible and irreversible rxns.
     */
    if (r.reversible) {
      m_revindex.push_back(reactionNumber());
      m_nrev++;
    } else {
      m_irrev.push_back( reactionNumber() );
      m_nirrev++;
    }        
  }

    /**
     * Prepare the class for the addition of reactions. This function
     * must be called after instantiation of the class, but before
     * any reactions are actually added to the mechanism.
     * This function calculates m_kk the number of species in all
     * phases participating in the reaction mechanism. We don't know
     * m_kk previously, before all phases have been added. 
     */
    void InterfaceKinetics::init() {
        int n;
        m_kk = 0;
        int np = nPhases();
        for (n = 0; n < np; n++) {
            m_kk += thermo(n).nSpecies();
        }
        m_rrxn.resize(m_kk);
        m_prxn.resize(m_kk);
        m_conc.resize(m_kk);
        m_mu0.resize(m_kk);
        m_grt.resize(m_kk);
        m_pot.resize(m_kk, 0.0);
        m_phi.resize(np, 0.0);
    }

    /**
     * Finish adding reactions and prepare for use. This function
     * must be called after all reactions are entered into the mechanism
     * and before the mechanism is used to calculate reaction rates.
     *
     * Here, we resize work arrays based on the number of reactions,
     * since we don't know this number up to now.
     */
    void InterfaceKinetics::finalize() {
        m_rwork.resize(nReactions());
        int ks = reactionPhaseIndex();
        if (ks < 0) throw CanteraError("InterfaceKinetics::finalize",
             "no surface phase is present.");
        m_surf = (SurfPhase*)&thermo(ks);
        if (m_surf->nDim() != 2) 
            throw CanteraError("InterfaceKinetics::finalize",
                "expected interface dimension = 2, but got dimension = "
                +int2str(m_surf->nDim()));
        m_finalized = true;
    }


  doublereal InterfaceKinetics::electrochem_beta(int irxn) const{
    int n = m_ctrxn.size();
    for (int i = 0; i < n; i++) {
      if (m_ctrxn[i] == irxn) {
	return m_beta[i];
      }
    }
    return 0.0;
  }


    bool InterfaceKinetics::ready() const {
        return (m_finalized);
    }

  // Advance the surface coverages in time
  /*
   * @param tstep  Time value to advance the surface coverages
   */
  void InterfaceKinetics::
  advanceCoverages(doublereal tstep) {
    if (m_integrator == 0) {
      vector<InterfaceKinetics*> k;
      k.push_back(this);
      m_integrator = new ImplicitSurfChem(k);
      m_integrator->initialize();
    }
    m_integrator->integrate(0.0, tstep);
    delete m_integrator;
    m_integrator = 0;
  }

  // Solve for the pseudo steady-state of the surface problem
  /*
   * Solve for the steady state of the surface problem.
   * This is the same thing as the advanceCoverages() function,
   * but at infinite times.
   *
   * Note, a direct solve is carried out under the hood here,
   * to reduce the computational time.
   *
   * the integrator object is saved inbetween calls to 
   * reduce the computational cost of repeated calls.
   */
  void InterfaceKinetics::
  solvePseudoSteadyStateProblem(int ifuncOverride,
				doublereal timeScaleOverride) {
    // create our own solver object
    if (m_integrator == 0) {
      vector<InterfaceKinetics*> k;
      k.push_back(this);
      m_integrator = new ImplicitSurfChem(k);
      m_integrator->initialize();
    }
    m_integrator->setIOFlag(m_ioFlag);
    /*
     * New direct method to go here
     */
    m_integrator->solvePseudoSteadyStateProblem(ifuncOverride, timeScaleOverride);
  }

  void EdgeKinetics::finalize() {
    m_rwork.resize(nReactions());
    int ks = reactionPhaseIndex();
    if (ks < 0) throw CanteraError("EdgeKinetics::finalize",
				   "no edge phase is present.");
    m_surf = (SurfPhase*)&thermo(ks);
    if (m_surf->nDim() != 1) 
      throw CanteraError("EdgeKinetics::finalize",
			 "expected interface dimension = 1, but got dimension = "
			 +int2str(m_surf->nDim()));
    m_finalized = true;
  }
}


