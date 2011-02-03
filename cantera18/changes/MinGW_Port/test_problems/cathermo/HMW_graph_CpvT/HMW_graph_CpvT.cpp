/**
 *
 *  @file HMW_graph_1.cpp
 */

/*
 *  $Author: hkmoffa $
 *  $Date: 2008/12/17 17:31:12 $
 *  $Revision: 1.5 $
 */

#ifdef SRCDIRTREE
#include "ct_defs.h"
#include "logger.h"
#include "TemperatureTable.h"
#include "ThermoPhase.h"
#include "HMWSoln.h"
#include "importCTML.h"
#else
#include "Cantera.h"
#include "kernel/logger.h"
#include "thermo.h"
#include "TemperatureTable.h"
#include "HMWSoln.h"
#endif

#include <cstdio>

using namespace std;
using namespace Cantera;

class fileLog: public Logger {
public:
    fileLog(string fName) {
	m_fName = fName;
	m_fs.open(fName.c_str());
    }
    
    virtual void write(const string& msg) {
          m_fs << msg;
	  m_fs.flush();
    }

    virtual ~fileLog() {
	m_fs.close();
    }
    
    string m_fName;
    ofstream m_fs;

};

void printUsage() {
    cout << "usage: HMW_test " <<  endl;
    cout <<"                -> Everything is hardwired" << endl;
}

void pAtable(HMWSoln *HMW) {
    int nsp = HMW->nSpecies();
    double acMol[100]; 
    double mf[100];
    double activities[100];
    double moll[100];
    
    HMW->getMolalityActivityCoefficients(acMol);
    HMW->getMoleFractions(mf);
    HMW->getActivities(activities);
    HMW->getMolalities(moll);
    string sName;
    printf("            Name      Activity  ActCoeffMolal "
	   "   MoleFract      Molality\n");
    for (int k = 0; k < nsp; k++) {
      sName = HMW->speciesName(k);
      printf("%16s %13g %13g %13g %13g\n", 
	     sName.c_str(), activities[k], acMol[k], mf[k], moll[k]);
    }
  
}

int main(int argc, char **argv)
{

  int retn = 0;
  int i;
 
  try {

    char iFile[80];
    strcpy(iFile, "HMW_NaCl.xml");
    if (argc > 1) {
      strcpy(iFile, argv[1]);
    }
    double Cp0_R[20], pmCp[20];

    //fileLog *fl = new fileLog("HMW_graph_1.log");
    //setLogger(fl);

    HMWSoln *HMW = new HMWSoln(iFile, "NaCl_electrolyte");


    /*
     * Load in and initialize the 
     */
    Cantera::ThermoPhase *solid = newPhase("NaCl_Solid.xml","NaCl(S)");
 
     
    int nsp = HMW->nSpecies();
    double acMol[100]; 
    double act[100];
    double mf[100];
    double moll[100];
    for (i = 0; i < 100; i++) {
      acMol[i] = 1.0;
      act[i] = 1.0;
      mf[i] = 0.0;
      act[i] = 0.0;
    }

    HMW->getMoleFractions(mf);
    string sName;

    TemperatureTable TTable(15, false, 273.15, 25., 0, 0);


    HMW->setState_TP(298.15, 1.01325E5);
  
    int i1 = HMW->speciesIndex("Na+");
    int i2 = HMW->speciesIndex("Cl-");
    //int i3 = HMW->speciesIndex("H2O(L)");
    for (i = 0; i < nsp; i++) {
      moll[i] = 0.0;
    }
    HMW->setMolalities(moll);

    double ISQRT;
    double Is = 0.0;

    /*
     * Set the Pressure
     */
    double pres = OneAtm;

    /*
     * Fix the molality
     */
    Is = 6.146;
    ISQRT = sqrt(Is);
    moll[i1] = Is;
    moll[i2] = Is;
    HMW->setState_TPM(298.15, pres, moll);
    double Xmol[30];
    HMW->getMoleFractions(Xmol);

    /*
     * ThermoUnknowns
     */
    double T;
 
    double Cp0_NaCl = 0.0, Cp0_Naplus = 0.0, Cp0_Clminus = 0.0, Delta_Cp0s = 0.0, Cp0_H2O = 0.0;
    double Cp_NaCl = 0.0, Cp_Naplus = 0.0, Cp_Clminus = 0.0, Cp_H2O = 0.0;
    double molarCp0;
#ifdef DEBUG_HKM
    FILE *ttt = fopen("table.csv","w");
#endif
    printf("A_J/R: Comparison to Pitzer's book, p. 99, can be made.\n");
    printf("        Agreement is within 12 pc \n");
    printf("\n");

    printf("Delta_Cp0: Heat Capacity of Solution per mole of salt (standard states)\n");
    printf("           rxn for the ss heat of soln:     "
	   "NaCl(s) -> Na+(aq) + Cl-(aq)\n");

    printf("\n");
    printf("Delta_Cps: Delta heat Capacity of Solution per mole of salt\n");
    printf("          rxn for heat of soln:     "
	   " n1 H2O(l,pure) + n2 NaCl(s) -> n2 MX(aq) + n1 H2O(l) \n");
    printf("          Delta_Hs = (n1 h_H2O_bar + n2 h_MX_bar "
	   "- n1 h_H2O_0 - n2 h_MX_0)/n2\n");
    printf("\n");
    printf("phiJ:     phiJ, calculated from the program, is checked\n");
    printf("          against analytical formula in J_standalone program.\n");
    printf("          (comparison against Eq. 12, Silvester and Pitzer)\n");

    /*
     * Create a Table of NaCl Enthalpy Properties as a Function
     * of the Temperature
     */
    printf("\n\n");
    printf("            T,          Pres,          Aphi,         A_J/R,"
	   "     Delta_Cp0,"
	   "     Delta_Cps,             J,          phiJ,"
	   "     MolarCp,   MolarCp0\n");
    printf("       Kelvin,           bar, sqrt(kg/gmol), sqrt(kg/gmol),"
	   "   kJ/gmolSalt,"
	   "   kJ/gmolSalt,   kJ/gmolSoln,   kJ/gmolSalt,"
	   "       kJ/gmol,    kJ/gmol\n");
#ifdef DEBUG_HKM
    fprintf(ttt,"T, Pres, A_J/R, Delta_Cp0, Delta_Cps, J, phiJ\n");
    fprintf(ttt,"Kelvin, bar, sqrt(kg/gmol), kJ/gmolSalt, kJ/gmolSalt, kJ/gmolSoln,"
	    "kJ/gmolSalt\n");
#endif
    for (i = 0; i < TTable.NPoints + 1; i++) {
      if (i == TTable.NPoints) {
	T = 323.15;
      } else {
	T = TTable.T[i];
      }
      /*
       * RT is in units of J/kmolK
       */
      //double RT = GasConstant * T;

      /*
       * Make sure we are at the saturation pressure or above.
       */
 
      double psat = HMW->satPressure(T);
 
      pres = OneAtm;
      if (psat > pres) pres = psat;


      HMW->setState_TPM(T, pres, moll);

      solid->setState_TP(T, pres);
 
      /*
       * Get the Standard State DeltaH
       */
   
      solid->getCp_R(Cp0_R);
      Cp0_NaCl = Cp0_R[0] * GasConstant * 1.0E-6;

  
      HMW->getCp_R(Cp0_R);
      Cp0_H2O    = Cp0_R[0] * GasConstant * 1.0E-6;
      Cp0_Naplus = Cp0_R[i1] * GasConstant * 1.0E-6;
      Cp0_Clminus = Cp0_R[i2] * GasConstant * 1.0E-6;
      /*
       * Calculate the standard state heat of solution
       * for NaCl(s) -> Na+ + Cl-
       *   units: kJ/gmolSalt
       */
 
      Delta_Cp0s = Cp0_Naplus + Cp0_Clminus - Cp0_NaCl;

      pmCp[0] = solid->cp_mole();
  
      Cp_NaCl = pmCp[0] * 1.0E-6;


      HMW->getPartialMolarCp(pmCp);
      Cp_H2O     = pmCp[0]  * 1.0E-6;
      Cp_Naplus  = pmCp[i1] * 1.0E-6;
      Cp_Clminus = pmCp[i2] * 1.0E-6;

      //double Delta_Cp_Salt = Cp_NaCl - (Cp_Naplus + Cp_Clminus);

      double molarCp = HMW->cp_mole() * 1.0E-6; 

      /*
       * Calculate the heat capacity of solution for the reaction
       * NaCl(s) -> Na+ + Cl-
       */
      double Delta_Cps = (Xmol[0]  * Cp_H2O +
			 Xmol[i1] * Cp_Naplus +
			 Xmol[i2] * Cp_Clminus 
			 - Xmol[0] * Cp0_H2O
			 - Xmol[i1] * Cp_NaCl);
      Delta_Cps /= Xmol[i1];
 
   
      /*
       * Calculate the relative heat capacity, J, from the
       * partial molar quantities, units J/gmolSolutionK
       */
      double J = (Xmol[0]  * (Cp_H2O    - Cp0_H2O) +
		  Xmol[i1] * (Cp_Naplus - Cp0_Naplus) +
		  Xmol[i2] * (Cp_Clminus - Cp0_Clminus));

      /*
       * Calculate the apparent relative molal heat capacity, phiJ,
       * units of J/gmolSaltAddedK
       */
      double phiJ = J / Xmol[i1];
		   

      double Aphi = HMW->A_Debye_TP(T, pres) / 3.0;
      //double AL = HMW->ADebye_L(T,pres);
      double AJ = HMW->ADebye_J(T, pres);



      for (int k = 0; k < nsp; k++) {
	Cp0_R[k] *= GasConstant * 1.0E-6;
      }

      molarCp0 = 0.0;
      for (int k = 0; k < nsp; k++) {
	molarCp0 += Xmol[k] * Cp0_R[k];
      }

      if (i != TTable.NPoints+1) {
	printf("%13g, %13g, %13g, %13g, %13g, %13g, "
	     "%13g, %13g, %13g, %13g\n",
	     T, pres*1.0E-5,  Aphi, AJ/GasConstant, Delta_Cp0s, Delta_Cps,
	     J, phiJ, molarCp , molarCp0 );
#ifdef DEBUG_HKM
	fprintf(ttt,"%g, %g, %g, %g, %g, %g, %g\n",
		T, pres*1.0E-5, AJ/GasConstant, Delta_Cp0s, Delta_Cps, J, phiJ);
#endif
      }
 
    }

    printf("Breakdown of Heat Capacity Calculation at 323.15 K, 1atm:\n");
    
    printf(" Species     MoleFrac        Molal          Cp0      "
	   "    partCp     (partCp - Cp0)\n");
    printf("  H2O(L)");
    printf("%13g %13g %13g %13g %13g\n", Xmol[0], moll[0], Cp0_H2O , Cp_H2O,  Cp_H2O-Cp0_H2O);
    printf("  Na+   ");
    printf("%13g %13g %13g %13g %13g\n", Xmol[i1], moll[i1], 
	   Cp0_Naplus , Cp_Naplus,  Cp_Naplus -Cp0_Naplus);
    printf("  Cl-   ");
    printf("%13g %13g %13g %13g %13g\n", Xmol[i2], moll[i2], 
	   Cp0_Clminus , Cp_Clminus,  Cp_Clminus - Cp0_Clminus);

    printf(" NaCl(s)");
    printf("%13g               %13g %13g %13g\n", 1.0, 
	   Cp0_NaCl , Cp_NaCl,  Cp_NaCl - Cp0_NaCl);

    
    delete HMW;
    HMW = 0;
    delete solid;
    solid = 0;
    Cantera::appdelete();

#ifdef DEBUG_HKM
    fclose(ttt);
#endif
    return retn;

  } catch (CanteraError) {
    printf("caught error\n");
    showErrors();
    Cantera::appdelete();
    return -1;
  }
} 
