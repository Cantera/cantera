/**
 *   @file PureFluidPhase.cpp
 *   Definitions for a ThermoPhase object for a pure fluid phase consisting 
 *   of gas, liquid, mixed-gas-liquid
 *   and supercritical fluid (see \ref thermoprops 
 *   and class \link Cantera::PureFluidPhase PureFluidPhase\endlink).
 */
/*
 *  $Id: PureFluidPhase.cpp,v 1.5 2009/03/13 03:21:34 hkmoffa Exp $
 */
#include "xml.h"
#include "PureFluidPhase.h"

#ifdef WITH_PURE_FLUIDS

#include "../../../ext/tpx/Sub.h"
#include "../../../ext/tpx/utils.h"

#include <cstdlib>

namespace Cantera {

  // Base Constructor
  PureFluidPhase::PureFluidPhase() :
    ThermoPhase(), 
    m_sub(0),
    m_subflag(0), 
    m_mw(-1.0), 
    m_verbose(false)
  {
  }

  // CopyConstructor
  PureFluidPhase::PureFluidPhase(const PureFluidPhase& right) :
    ThermoPhase(), 
    m_sub(0),
    m_subflag(0), 
    m_mw(-1.0), 
    m_verbose(false)
  {
    *this = right;
  }

  //! Assignment operator
  /*!
   * @param right Object to be copied
   */
  PureFluidPhase& PureFluidPhase::operator=(const PureFluidPhase& right) {
    if (&right != this) {
      ThermoPhase::operator=(right);
      if (m_sub) {
	delete m_sub;
      }
      m_subflag    = right.m_subflag;
      m_sub        = tpx::GetSub(m_subflag);
      m_mw         = right.m_mw;
      m_verbose    = right.m_verbose;
    }
    return *this;
  }
  
  // Duplicator from the %ThermoPhase parent class
  /*
   * Given a pointer to a %ThermoPhase object, this function will
   * duplicate the %ThermoPhase object and all underlying structures.
   * This is basically a wrapper around the copy constructor.
   *
   * @return returns a pointer to a %ThermoPhase
   */
  ThermoPhase *PureFluidPhase::duplMyselfAsThermoPhase() const {
    PureFluidPhase *igp = new PureFluidPhase(*this);
    return (ThermoPhase *) igp;
  }




    PureFluidPhase::~PureFluidPhase() {
      delete m_sub; 
    }

    void PureFluidPhase::
    initThermo() {
        if (m_sub) delete m_sub;
        m_sub = tpx::GetSub(m_subflag);
        if (m_sub == 0) {
            throw CanteraError("PureFluidPhase::initThermo",
                "could not create new substance object.");
        }
        m_mw = m_sub->MolWt();
        m_weight[0] = m_mw;
        setMolecularWeight(0,m_mw);
        double one = 1.0;
        setMoleFractions(&one);
        double cp0_R, h0_RT, s0_R, T0, p;
        T0 = 298.15;
        if (T0 < m_sub->Tcrit()) {
            m_sub->Set(tpx::TX, T0, 1.0);
            p = 0.01*m_sub->P();
        }
        else {
            p = 0.001*m_sub->Pcrit();
        }
        m_sub->Set(tpx::TP, T0, p);
        
        m_spthermo->update_one(0, T0, &cp0_R, &h0_RT, &s0_R);
        double s_R = s0_R - log(p/refPressure());
        m_sub->setStdState(h0_RT*GasConstant*298.15/m_mw,
            s_R*GasConstant/m_mw, T0, p);
        if (m_verbose) {
            writelog("PureFluidPhase::initThermo: initialized phase "
                +id()+"\n");
        }
    }

    void PureFluidPhase::
    setParametersFromXML(const XML_Node& eosdata) {
        eosdata._require("model","PureFluid");
        m_subflag = atoi(eosdata["fluid_type"].c_str());
        if (m_subflag < 0) 
            throw CanteraError("PureFluidPhase::setParametersFromXML",
                "missing or negative substance flag");
    }

    doublereal PureFluidPhase::
    enthalpy_mole() const {
        setTPXState();
        doublereal h = m_sub->h() * m_mw;
        check(h);
        return h;            
    }
        
    doublereal PureFluidPhase::
    intEnergy_mole() const {
        setTPXState();
        doublereal u = m_sub->u() * m_mw;
        check(u);
        return u;            
    }

    doublereal PureFluidPhase::
    entropy_mole() const {
        setTPXState();
        doublereal s = m_sub->s() * m_mw;
        check(s);
        return s;            
    }

    doublereal PureFluidPhase::
    gibbs_mole() const {
        setTPXState();
        doublereal g = m_sub->g() * m_mw;
        check(g);
        return g;            
    }

    doublereal PureFluidPhase::
    cp_mole() const {
        setTPXState();
        doublereal cp = m_sub->cp() * m_mw;
        check(cp);
        return cp;            
    }

    doublereal PureFluidPhase::
    cv_mole() const {
        setTPXState();
        doublereal cv = m_sub->cv() * m_mw;
        check(cv);
        return cv;
    }
        
    doublereal PureFluidPhase::
    pressure() const {
        setTPXState();
        doublereal p = m_sub->P();
        check(p);
        return p;
    }
        
    void PureFluidPhase::
    setPressure(doublereal p) {
        Set(tpx::TP, temperature(), p);
        setDensity(1.0/m_sub->v());
        check();
    }

    void PureFluidPhase::Set(int n, double x, double y) const {
        try { 
            m_sub->Set(n, x, y); 
        }
        catch(tpx::TPX_Error) {
            reportTPXError();
        }
    }
            
    void PureFluidPhase::setTPXState() const {
        Set(tpx::TV, temperature(), 1.0/density());
    }
        
    void PureFluidPhase::check(doublereal v) const {
        if (m_sub->Error() || v == tpx::Undef) {
            throw CanteraError("PureFluidPhase",string(tpx::errorMsg(
                                                           m_sub->Error())));
        }
    }

    void PureFluidPhase::reportTPXError() const {
        string msg = tpx::TPX_Error::ErrorMessage;
        string proc = "tpx::"+tpx::TPX_Error::ErrorProcedure;
        throw CanteraError(proc,msg);
    }


    doublereal PureFluidPhase::isothermalCompressibility() const {
        return m_sub->isothermalCompressibility();
    }

    doublereal PureFluidPhase::thermalExpansionCoeff() const {
        return m_sub->thermalExpansionCoeff();
    }

    tpx::Substance& PureFluidPhase::TPX_Substance() { return *m_sub; }

    /// critical temperature 
    doublereal PureFluidPhase::critTemperature() const { return m_sub->Tcrit(); }
        
    /// critical pressure
    doublereal PureFluidPhase::critPressure() const { return m_sub->Pcrit(); }
        
    /// critical density
    doublereal PureFluidPhase::critDensity() const { return 1.0/m_sub->Vcrit(); }
        
        
    /// saturation temperature
    doublereal PureFluidPhase::satTemperature(doublereal p) const { 
        try {
            doublereal ts = m_sub->Tsat(p);
            return ts;
        }
        catch(tpx::TPX_Error) {
            reportTPXError();
            return -1.0;
        }
    }
        
    void PureFluidPhase::setState_HP(doublereal h, doublereal p, 
        doublereal tol) {
        Set(tpx::HP, h, p);
        setState_TR(m_sub->Temp(), 1.0/m_sub->v());
        check();
    }

    void PureFluidPhase::setState_UV(doublereal u, doublereal v, 
        doublereal tol) {
        Set(tpx::UV, u, v);
        setState_TR(m_sub->Temp(), 1.0/m_sub->v());
        check();
    }

    void PureFluidPhase::setState_SV(doublereal s, doublereal v, 
        doublereal tol) {
        Set(tpx::SV, s, v);
        setState_TR(m_sub->Temp(), 1.0/m_sub->v());
        check();
    }

    void PureFluidPhase::setState_SP(doublereal s, doublereal p, 
        doublereal tol) {
        Set(tpx::SP, s, p);
        setState_TR(m_sub->Temp(), 1.0/m_sub->v());
        check();
    }

    /// saturation pressure
    doublereal PureFluidPhase::satPressure(doublereal t) const {
        doublereal vsv = m_sub->v();
        try {
            Set(tpx::TV,t,vsv);
            doublereal ps = m_sub->Ps();
            return ps;
        }
        catch(tpx::TPX_Error) {
            reportTPXError();
            return -1.0;
        }
    }
        
    doublereal PureFluidPhase::vaporFraction() const {
        setTPXState();
        doublereal x = m_sub->x();
        check(x);
        return x;
    }
        
    void PureFluidPhase::setState_Tsat(doublereal t, doublereal x) {
        setTemperature(t);
        setTPXState();
        Set(tpx::TX, t, x);
        setDensity(1.0/m_sub->v());
        check();
    }

    void PureFluidPhase::setState_Psat(doublereal p, doublereal x) {
        setTPXState();
        Set(tpx::PX, p, x);
        setTemperature(m_sub->Temp());
        setDensity(1.0/m_sub->v());
        check();
    }


  /**
   * Format a summary of the mixture state for output.
   */           
  std::string PureFluidPhase::report(bool show_thermo) const {


    char p[800];
    string s = "";
    try {
      if (name() != "") {
	sprintf(p, " \n  %s:\n", name().c_str());
	s += p;
      }
      sprintf(p, " \n       temperature    %12.6g  K\n", temperature());
      s += p;
      sprintf(p, "          pressure    %12.6g  Pa\n", pressure());
      s += p;
      sprintf(p, "           density    %12.6g  kg/m^3\n", density());
      s += p;
      sprintf(p, "  mean mol. weight    %12.6g  amu\n", meanMolecularWeight());
      s += p;

      if (eosType() == cPureFluid) {
      	double xx = ((PureFluidPhase *) (this))->vaporFraction();
      	sprintf(p, "    vapor fraction    %12.6g \n", 
		xx); //th.vaporFraction());
      	s += p;
      }

      doublereal phi = electricPotential();
      if (phi != 0.0) {
	sprintf(p, "         potential    %12.6g  V\n", phi);
	s += p;
      }
      if (show_thermo) {
        sprintf(p, " \n");
        s += p;
        sprintf(p, "                          1 kg            1 kmol\n");
        s += p;
        sprintf(p, "                       -----------      ------------\n");
        s += p;
        sprintf(p, "          enthalpy    %12.6g     %12.4g     J\n", 
		enthalpy_mass(), enthalpy_mole());
        s += p;
        sprintf(p, "   internal energy    %12.6g     %12.4g     J\n", 
		intEnergy_mass(), intEnergy_mole());
        s += p;
        sprintf(p, "           entropy    %12.6g     %12.4g     J/K\n", 
		entropy_mass(), entropy_mole());
        s += p;
        sprintf(p, "    Gibbs function    %12.6g     %12.4g     J\n", 
		gibbs_mass(), gibbs_mole());
        s += p;
        sprintf(p, " heat capacity c_p    %12.6g     %12.4g     J/K\n", 
		cp_mass(), cp_mole());
        s += p;
        try {
	  sprintf(p, " heat capacity c_v    %12.6g     %12.4g     J/K\n", 
		  cv_mass(), cv_mole());
	  s += p;
        }
        catch(CanteraError) {
	  sprintf(p, " heat capacity c_v    <not implemented>       \n");
	  s += p;
        }
      }

      int kk = nSpecies();
      array_fp x(kk);
      array_fp y(kk);
      array_fp mu(kk);
      getMoleFractions(&x[0]);
      getMassFractions(&y[0]);
      getChemPotentials(&mu[0]);
      doublereal rt = GasConstant * temperature(); 
      int k;
      //if (th.nSpecies() > 1) {

      if (show_thermo) {
	sprintf(p, " \n                           X     "
		"            Y          Chem. Pot. / RT    \n");
	s += p;
	sprintf(p, "                     -------------     "
		"------------     ------------\n");
	s += p;
	for (k = 0; k < kk; k++) {
	  if (x[k] > SmallNumber) {
	    sprintf(p, "%18s   %12.6g     %12.6g     %12.6g\n", 
		    speciesName(k).c_str(), x[k], y[k], mu[k]/rt);
	  }
	  else {
	    sprintf(p, "%18s   %12.6g     %12.6g     \n", 
		    speciesName(k).c_str(), x[k], y[k]);
	  }
	  s += p;
	}
      }
      else {
	sprintf(p, " \n                           X"
		"Y\n");
	s += p;
	sprintf(p, "                     -------------"
		"     ------------\n");
	s += p;
	for (k = 0; k < kk; k++) {
	  sprintf(p, "%18s   %12.6g     %12.6g\n", 
		  speciesName(k).c_str(), x[k], y[k]);
	  s += p;
	}
      }
    }
    //}
    catch (CanteraError) {
      ;
    }
    return s;
  }

 
}

#endif  // WITH_PURE_FLUIDS
