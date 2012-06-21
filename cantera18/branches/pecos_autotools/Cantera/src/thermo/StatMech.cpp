/**
 *  @file StatMech.cpp
 *  \link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType\endlink 
 */

/* $Author: hkmoffa $
 * $Revision: 279 $
 * $Date: 2009-12-05 13:08:43 -0600 (Sat, 05 Dec 2009) $
 */
// Copyright 2007  Sandia National Laboratories

#include "StatMech.h"

namespace Cantera {

 
  // Statistical mechanics 
  /*
   * @ingroup spthermo
   */
 

  //! Empty constructor
  StatMech::StatMech() 
    : m_lowT(0.0), m_highT (0.0),
      m_Pref(1.0E5), m_index (0) {}


  // constructor used in templated instantiations
  /*
   * @param n            Species index
   * @param tlow         Minimum temperature
   * @param thigh        Maximum temperature
   * @param pref         reference pressure (Pa).
   * @param coeffs       Vector of coefficients used to set the
   *                     parameters for the standard state.
   */
  StatMech::StatMech(int n, doublereal tlow, doublereal thigh, 
			 doublereal pref,
			 const doublereal* coeffs) :
    m_lowT      (tlow),
    m_highT     (thigh),
    m_Pref      (pref),
    m_index     (n)
  {
    // should error on zero -- cannot take ln(0)
    if(tlow <= 0.0){
      throw CanteraError("Error in StatMech.cpp",
			 " Cannot take 0 tmin as input. \n\n");
    }    

  }

  // copy constructor
  /*
   * @param b object to be copied
   */
  StatMech::StatMech(const StatMech& b) :
    m_lowT      (b.m_lowT),
    m_highT     (b.m_highT),
    m_Pref      (b.m_Pref),
    m_index     (b.m_index)
  {

  }

  // assignment operator
  /*
   * @param b object to be copied
   */
  StatMech& StatMech::operator=(const StatMech& b) {
    if (&b != this) {
      m_lowT   = b.m_lowT;
      m_highT  = b.m_highT;
      m_Pref   = b.m_Pref;
      m_index  = b.m_index;
    }
    return *this;
  }

  // Destructor
  StatMech::~StatMech() {
  }

  // duplicator
  SpeciesThermoInterpType *
  StatMech::duplMyselfAsSpeciesThermoInterpType() const {
    StatMech* np = new StatMech(*this);
    return (SpeciesThermoInterpType *) np;
  }

  // Returns the minimum temperature that the thermo
  // parameterization is valid
  doublereal StatMech::minTemp() const     
  { 
    return m_lowT;    
  }

  // Returns the maximum temperature that the thermo
  // parameterization is valid
  doublereal StatMech::maxTemp() const  { 
    return m_highT;
  }

  // Returns the reference pressure (Pa)
  doublereal StatMech::refPressure() const { return m_Pref; }

  // Returns an integer representing the type of parameterization
  int StatMech::reportType() const {
    return NASA9;
  }
      
  // Returns an integer representing the species index
  int StatMech::speciesIndex() const { 
    return m_index;
  }

  int StatMech::buildmap()
  {

    // build map
    name_map["Air"  ] = 1; 
    name_map["CPAir"] = 2; 
    name_map["Ar"   ] = 3;   
    name_map["Ar+"  ] = 4;  
    name_map["C"    ] = 5;    
    name_map["C+"   ] = 6;   
    name_map["C2"   ] = 7;   
    name_map["C2H"  ] = 8;
    name_map["C2H2" ] = 9;
    name_map["C3"   ] =10; 
    name_map["CF"   ] =11; 
    name_map["CF2"  ] =12; 
    name_map["CF3"  ] =13; 
    name_map["CF4"  ] =14; 
    name_map["CH"   ] =15; 
    name_map["CH2"  ] =16; 
    name_map["CH3"  ] =17; 
    name_map["CH4"  ] =18; 
    name_map["Cl"   ] =19; 
    name_map["Cl2"  ] =20; 
    name_map["CN"   ] =21; 
    name_map["CN+"  ] =22; 
    name_map["CO"   ] =23; 
    name_map["CO+"  ] =24; 
    name_map["CO2"  ] =25; 
    name_map["F"    ] =26;
    name_map["F2"   ] =27; 
    name_map["H"    ] =28; 
    name_map["H+"   ] =29; 
    name_map["H2"   ] =30; 
    name_map["H2+"  ] =31; 
    name_map["H2O"  ] =32; 
    name_map["HCl"  ] =33; 
    name_map["HCN"  ] =34; 
    name_map["He"   ] =35; 
    name_map["He+"  ] =36; 
    name_map["N"    ] =37; 
    name_map["N+"   ] =38; 
    name_map["N2"   ] =39; 
    name_map["CPN2" ] =40; 
    name_map["N2+"  ] =41; 
    name_map["Ne"   ] =42; 
    name_map["NCO"  ] =43; 
    name_map["NH"   ] =44; 
    name_map["NH+"  ] =45; 
    name_map["NH2"  ] =46; 
    name_map["NH3"  ] =47; 
    name_map["NO"   ] =48; 
    name_map["NO+"  ] =49; 
    name_map["NO2"  ] =50; 
    name_map["O"    ] =51; 
    name_map["O+"   ] =52; 
    name_map["O2"   ] =53; 
    name_map["O2+"  ] =54; 
    name_map["OH"   ] =55; 
    name_map["Si"   ] =56; 
    name_map["SiO"  ] =57; 
    name_map["e"    ] =58; 
 
    // build species struct data
    Air.cfs = 10;
    Air.mol_weight=28.96;
    Air.species_name="Air";
    Air.nvib=0;
  
    // C2
    C2.cfs=7;
    C2.mol_weight=24.022;
    C2.species_name="C2";
    C2.nvib=1;
    C2.theta=2.6687e3;

    // H2
    H2.cfs  = 20;

    // build vector
    svec[1]=Air;
    svec[2]=H2;
    svec[3]=C2;


  return 0;
  }


  // void StatMech::read_species_vibrational_table ()
  // {
  //   // ("#===========================================================================\n"
  //   //  "# LEGEND\n"
  //   //  "#===========================================================================\n"
  //   //  "#\n"
  //   //  "# Species 	-- Species name\n"
  //   //  "# Theta_v	-- Characteristic vibrational temperature (K)\n"
  //   //  "# degeneracy   -- degeneracy of the mode\n"
  //   //  "#\n"
  //   //  "# Characteristic Temperatures for Simple Harmonic Oscillator Model\n"
  //   //  "#\n"
  //   //  "# Species     Theta_v  degeneracy\n"

    

    
  //   istringstream default_species_vib_data(
  //      "C2      2.66870e+03   1\n"
  //      "C2H     5.20100e+03   1\n"
  //      "C2H     7.20000e+02   2\n"
  //      "C2H     2.66100e+03   1\n"
  //      "C2H2    4.85290e+03   1\n"
  //      "C2H2    2.84000e+03   1\n"
  //      "C2H2    4.72490e+03   1\n"
  //      "C2H2    8.81830e+02   2\n"
  //      "C2H2    1.05080e+03   2\n"
  //      "C3      1.84500e+03   1\n"
  //      "C3      7.78700e+02   2\n"
  //      "C3      3.11760e+03   1\n"
  //      "CF      1.88214e+03   1\n"
  //      "CF2     1.76120e+03   1\n"
  //      "CF2     9.56820e+02   1\n"
  //      "CF2     1.60000e+03   1\n"
  //      "CF3     1.56800e+03   1\n"
  //      "CF3     1.00900e+03   1\n"
  //      "CF3     1.81150e+03   2\n"
  //      "CF3     7.36680e+02   2\n"
  //      "CF4     1.30720e+03   1\n"
  //      "CF4     6.25892e+02   2\n"
  //      "CF4     1.84540e+03   3\n"
  //      "CF4     9.08950e+02   3\n"
  //      "CH      4.11290e+03   1\n"
  //      "CH2     4.31650e+03   1\n"
  //      "CH2     1.95972e+03   1\n"
  //      "CH2     4.60432e+03   1\n"
  //      "CH3     4.31650e+03   1\n"
  //      "CH3     8.73370e+02   1\n"
  //      "CH3     4.54960e+03   2\n"
  //      "CH3     2.01150e+03   2\n"
  //      "CH4     4.19660e+03   1\n"
  //      "CH4     2.20620e+03   2\n"
  //      "CH4     4.34450e+03   3\n"
  //      "CH4     1.88600e+03   3\n"
  //      "Cl2     8.05355e+02   1\n"
  //      "CN      2.97610e+03   1\n"
  //      "CN+     2.92520e+03   1\n"
  //      "CO      3.12200e+03   1\n"
  //      "CO+     3.18800e+03   1\n"
  //      "CO2     1.91870e+03   1\n"
  //      "CO2     9.59660e+02   2\n"
  //      "CO2     3.38210e+03   1\n"
  //      "F2      1.32020e+03   1\n"
  //      "H2      6.33140e+03   1\n"
  //      "H2+     3.34280e+03   1\n"
  //      "H2O     5.26130e+03   1\n"
  //      "H2O     2.29460e+03   1\n"
  //      "H2O     5.40395e+03   1\n"
  //      "HCl     4.30330e+03   1\n"
  //      "HCN     3.01620e+03   1\n"
  //      "HCN     1.02660e+03   2\n"
  //      "HCN     4.76450e+03   1\n"
  //      "N2      3.39500e+03   1\n"
  //      "N2+     3.17580e+03   1\n"
  //      "NCO     1.83600e+03   1\n"
  //      "NCO     7.67100e+02   2\n"
  //      "NCO     2.76800e+03   1\n"
  //      "NH      4.72240e+03   1\n"
  //      "NH3     4.78100e+03   1\n"
  //      "NH3     1.47040e+03   1\n"
  //      "NH3     4.95440e+03   2\n"
  //      "NH3     2.34070e+03   2\n"
  //      "NO      2.81700e+03   1\n"
  //      "NO+     3.42100e+03   1\n"
  //      "NO2     1.07900e+03   1\n"
  //      "NO2     1.90000e+03   1\n"
  //      "NO2     2.32700e+03   1\n"
  //      "O2      2.23900e+03   1\n"
  //      "O2+     2.74120e+03   1\n"
  //      "OH      5.37820e+03   1\n"
  //      "SiO     1.78640e+03   1\n");
    
  //   string line;
  //   string name;
  //   string ss1,ss2,ss3,ss4,sss; 
  //   int k;
  //   int i = 0;

  //   while (std::getline(default_species_vib_data, line)) 
  //     {

  // 	istringstream ss(line);      
  // 	std::getline(ss, ss1, ' ');
  // 	std::getline(ss, ss2, ' ');
  // 	std::getline(ss, ss3, ' ');       
  // 	std::getline(ss, ss4, ' ');       
  // 	name = ss1;

  // 	// now put coefficients in correct species
  // 	for (k = 0; k < m_nsp; k++) 
  // 	  {
  // 	    string sss = m_thermo->speciesName(k);

  // 	    // this is the right species index
  // 	    if(sss.compare(ss1) == 0)
  // 	      {
  // 	    	a[k] = atof(ss2.c_str());
  // 	    	b[k] = atof(ss3.c_str());
  // 	    	c[k] = atof(ss4.c_str());
    
  // 	    	// index
  // 	    	i++;
  // 	      }
  // 	    else // default to air
  // 	      {
  // 		throw CanteraError("StatMech Class",
  // 				   " Missing vibrational species data. \n\n");
  //   	      }
	    
  // 	  } // done with for loop
  //     }

  // }

  // Update the properties for this species
  /**
   *
   * \f[
   * \frac{C_p^0(T)}{R} = \frac{C_v^0(T)}{R} + 1 
   * \f]
   *
   * Where,
   * \f[
   * \frac{C_v^0(T)}{R} = \frac{C_v^{tr}(T)}{R} + \frac{C_v^{vib}(T)}{R}
   * \f]
   *
   *
   * @param tt      vector of temperature polynomials
   * @param cp_R    Vector of Dimensionless heat capacities.
   *                (length m_kk).
   * @param h_RT    Vector of Dimensionless enthalpies.
   *                (length m_kk).
   * @param s_R     Vector of Dimensionless entropies.
   *                (length m_kk).
   */
  void StatMech::updateProperties(const doublereal* tt, 
				    doublereal* cp_R, doublereal* h_RT,
				    doublereal* s_R) const {
    
    // translational + rotational specific heat
    doublereal ctr = 0.0;
    double atomicity = 0.0;

    // 5/2 * R for molecules, 3/2 * R for atoms
    if(atomicity < 1.0) // atom
      {
	ctr += 3/2 * GasConstant; 
      }
    else  // molecule
      {
	ctr += 5/2 * GasConstant; 
	doublereal theta = 1.0;
	ctr += GasConstant * theta* ( theta* exp(theta/tt[0])/(tt[0]*tt[0]))/((exp(theta/tt[0])-1) * (exp(theta/tt[0])-1));
      }   
    doublereal cpdivR = ctr/GasConstant + 1;

    // ACTUNG: fix enthalpy and entropy 
    doublereal hdivRT = 0.0;
    doublereal sdivR  = 0.0;

    // return the computed properties in the location in the output 
    // arrays for this species
    cp_R[m_index] = cpdivR;
    h_RT[m_index] = hdivRT;
    s_R[m_index]  = sdivR;
    //writelog("NASA9poly1: for species "+int2str(m_index)+", h_RT = "+
    //    fp2str(h)+"\n");
  }

 
  // Compute the reference-state property of one species
  /*
   * Given temperature T in K, this method updates the values of
   * the non-dimensional heat capacity at constant pressure,
   * enthalpy, and entropy, at the reference pressure, Pref
   * of one of the species. The species index is used
   * to reference into the cp_R, h_RT, and s_R arrays.
   *
   * Temperature Polynomial:
   *  tt[0] = t;
   *  tt[1] = t*t;
   *  tt[2] = t*t*t;
   *  tt[3] = t*t*t*t;
   *  tt[4] = 1.0/t;
   *  tt[5] = 1.0/(t*t);
   *  tt[6] = std::log(t);
   *
   * @param temp    Temperature (Kelvin)
   * @param cp_R    Vector of Dimensionless heat capacities.
   *                (length m_kk).
   * @param h_RT    Vector of Dimensionless enthalpies.
   *                (length m_kk).
   * @param s_R     Vector of Dimensionless entropies.
   *                (length m_kk).
   */
  void StatMech::updatePropertiesTemp(const doublereal temp, 
					doublereal* cp_R, doublereal* h_RT, 
					doublereal* s_R) const {
    double tPoly[7];
    tPoly[0]  = temp;
    tPoly[1]  = temp * temp;
    tPoly[2]  = tPoly[1] * temp;
    tPoly[3]  = tPoly[2] * temp;
    tPoly[4]  = 1.0 / temp;
    tPoly[5]  = tPoly[4] / temp;
    tPoly[6]  = std::log(temp);
    updateProperties(tPoly, cp_R, h_RT, s_R);
  }

  //This utility function reports back the type of 
  // parameterization and all of the parameters for the 
  // species, index.
  /*
   * All parameters are output variables
   *
   * @param n         Species index
   * @param type      Integer type of the standard type
   * @param tlow      output - Minimum temperature
   * @param thigh     output - Maximum temperature
   * @param pref      output - reference pressure (Pa).
   * @param coeffs    Vector of coefficients used to set the
   *                  parameters for the standard state.
   */
  void StatMech::reportParameters(int &n, int &type,
				    doublereal &tlow, doublereal &thigh,
				    doublereal &pref,
				    doublereal* const coeffs) const {
    n = m_index;
    type = NASA9;
    tlow = m_lowT;
    thigh = m_highT;
    pref = m_Pref;
    coeffs[0] = 1;
    coeffs[1] = m_lowT;
    coeffs[2] = m_highT;
    for (int i = 0; i < 9; i++) {
      coeffs[i+3] = 0.0;
    }

  }

  // Modify parameters for the standard state
  /*
   * @param coeffs   Vector of coefficients used to set the
   *                 parameters for the standard state.
   */
  void StatMech::modifyParameters(doublereal* coeffs) 
  {

    

  }


}

