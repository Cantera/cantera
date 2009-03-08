/**
 *  @file HMWSoln_input.cpp
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 * $Id: HMWSoln_input.cpp,v 1.4 2006/07/13 20:05:11 hkmoffa Exp $
 */

#include "HMWSoln.h"
#include "importCTML.h"
#include "WaterProps.h"
#include "WaterPDSS.h"

namespace Cantera {

  /**
   * interp_est()                          (static)
   *
   * utility function to assign an integer value from a string
   * for the ElectrolyteSpeciesType field.
   */
  static int interp_est(string estString) {
    const char *cc = estString.c_str();
    if (!strcasecmp(cc, "solvent")) {
      return cEST_solvent;
    } else if (!strcasecmp(cc, "chargedspecies")) {
      return cEST_chargedSpecies;
    } else if (!strcasecmp(cc, "weakAcidAssociated")) {
      return cEST_weakAcidAssociated;
    } else if (!strcasecmp(cc, "strongAcidAssociated")) {
      return cEST_strongAcidAssociated;
    } else if (!strcasecmp(cc, "polarNeutral")) {
      return cEST_polarNeutral;
    } else if (!strcasecmp(cc, "nonpolarNeutral")) {
      return cEST_nonpolarNeutral;
    }
    int retn, rval;
    if ((retn = sscanf(cc, "%d", &rval)) != 1) {
      return -1;
    }
    return rval;
  }

  /**
   * Process an XML node called "SimpleSaltParameters. 
   * This node contains all of the parameters necessary to describe
   * the Pitzer model for that particular binary salt.
   * This function reads the XML file and writes the coefficients
   * it finds to an internal data structures.
   */
  void HMWSoln::readXMLBinarySalt(XML_Node &BinSalt) {
    double *charge = DATA_PTR(m_speciesCharge);
    string stemp;
    int nParamsFound, i;
    vector_fp vParams;
    string iName = BinSalt.attrib("cation");
    if (iName == "") {
      throw CanteraError("HMWSoln::readXMLBinarySalt", "no cation attrib");
    }
    string jName = BinSalt.attrib("anion");
    if (jName == "") {
      throw CanteraError("HMWSoln::readXMLBinarySalt", "no anion attrib");
    }
    /*
     * Find the index of the species in the current phase. It's not
     * an error to not find the species
     */
    int iSpecies = speciesIndex(iName);
    if (iSpecies < 0) {
      return;
    }
    string ispName = speciesName(iSpecies);
    if (charge[iSpecies] <= 0) {
      throw CanteraError("HMWSoln::readXMLBinarySalt", "cation charge problem");
    }
    int jSpecies = speciesIndex(jName);
    if (jSpecies < 0) {
      return;
    }
    string jspName = speciesName(jSpecies);
    if (charge[jSpecies] >= 0) {
      throw CanteraError("HMWSoln::readXMLBinarySalt", "anion charge problem");
    }
	
    int n = iSpecies * m_kk + jSpecies;
    int counter = m_CounterIJ[n];
    int num = BinSalt.nChildren();
    for (int iChild = 0; iChild < num; iChild++) {
      XML_Node &xmlChild = BinSalt.child(iChild);
      stemp = xmlChild.name();
      string nodeName = lowercase(stemp);
      /*
       * Process the binary salt child elements
       */
      if (nodeName == "beta0") {
	/*
	 * Get the string containing all of the values
	 */
	getFloatArray(xmlChild, vParams, false, "", "beta0");
	nParamsFound = vParams.size();
	if (m_formPitzerTemp == PITZER_TEMP_CONSTANT) {
	  if (nParamsFound != 1) {
	    throw CanteraError("HMWSoln::readXMLBinarySalt::beta0 for " + ispName 
			       + "::" + jspName,
			       "wrong number of params found");
	  }
	  m_Beta0MX_ij[counter] = vParams[0];
	  m_Beta0MX_ij_coeff(0,counter) = m_Beta0MX_ij[counter];
	} else  if (m_formPitzerTemp == PITZER_TEMP_LINEAR) {
	  if (nParamsFound != 2) {
	    throw CanteraError("HMWSoln::readXMLBinarySalt::beta0 for " + ispName 
			       + "::" + jspName,
			       "wrong number of params found");
	  }
	  m_Beta0MX_ij_coeff(0,counter) = vParams[0];
	  m_Beta0MX_ij_coeff(1,counter) = vParams[1];
	  m_Beta0MX_ij[counter] = vParams[0];
	} else  if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
	  if (nParamsFound != 5) {
	    throw CanteraError("HMWSoln::readXMLBinarySalt::beta0 for " + ispName 
			       + "::" + jspName,
			       "wrong number of params found");
	  }
	  for (i = 0; i < nParamsFound; i++) {
	    m_Beta0MX_ij_coeff(i, counter) = vParams[i];
	  }
	  m_Beta0MX_ij[counter] = vParams[0];
	}
      }
      if (nodeName == "beta1") {

	/*
	 * Get the string containing all of the values
	 */
	getFloatArray(xmlChild, vParams, false, "", "beta1");
	nParamsFound = vParams.size();
	if (m_formPitzerTemp == PITZER_TEMP_CONSTANT) {
	  if (nParamsFound != 1) {
	    throw CanteraError("HMWSoln::readXMLBinarySalt::beta1 for " + ispName 
			       + "::" + jspName,
			       "wrong number of params found");
	  }
	  m_Beta1MX_ij[counter] = vParams[0];
	  m_Beta1MX_ij_coeff(0,counter) = m_Beta1MX_ij[counter];
	} else  if (m_formPitzerTemp == PITZER_TEMP_LINEAR) {
	  if (nParamsFound != 2) {
	    throw CanteraError("HMWSoln::readXMLBinarySalt::beta1 for " + ispName 
			       + "::" + jspName,
			       "wrong number of params found");
	  }
	  m_Beta1MX_ij_coeff(0,counter) = vParams[0];
	  m_Beta1MX_ij_coeff(1,counter) = vParams[1];
	  m_Beta1MX_ij[counter] = vParams[0];
	} else  if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
	  if (nParamsFound < 3) {
	    throw CanteraError("HMWSoln::readXMLBinarySalt::beta1 for " + ispName 
			       + "::" + jspName,
			       "wrong number of params found");
	  }
	  for (i = 0; i < nParamsFound; i++) {
	    m_Beta1MX_ij_coeff(i, counter) = vParams[i];
	  }
	  m_Beta1MX_ij[counter] = vParams[0];
	}

      }
      if (nodeName == "beta2") {
	stemp = xmlChild.value();
	m_Beta2MX_ij[counter] = atofCheck(stemp.c_str());
      }
      if (nodeName == "cphi") {
	/*
	 * Get the string containing all of the values
	 */
	getFloatArray(xmlChild, vParams, false, "", "Cphi");
	nParamsFound = vParams.size();
	if (m_formPitzerTemp == PITZER_TEMP_CONSTANT) {
	  if (nParamsFound != 1) {
	    throw CanteraError("HMWSoln::readXMLBinarySalt::Cphi for " + ispName 
			       + "::" + jspName,
			       "wrong number of params found");
	  }
	  m_CphiMX_ij[counter] = vParams[0];
	  m_CphiMX_ij_coeff(0,counter) = m_CphiMX_ij[counter];
	} else  if (m_formPitzerTemp == PITZER_TEMP_LINEAR) {
	  if (nParamsFound != 2) {
	    throw CanteraError("HMWSoln::readXMLBinarySalt::Cphi for " + ispName 
			       + "::" + jspName,
			       "wrong number of params found");
	  }
	  m_CphiMX_ij_coeff(0,counter) = vParams[0];
	  m_CphiMX_ij_coeff(1,counter) = vParams[1];
	  m_CphiMX_ij[counter] = vParams[0];
	} else  if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
	  if (nParamsFound != 5) {
	    throw CanteraError("HMWSoln::readXMLBinarySalt::Cphi for " + ispName 
			       + "::" + jspName,
			       "wrong number of params found");
	  }
	  for (i = 0; i < nParamsFound; i++) {
	    m_CphiMX_ij_coeff(i, counter) = vParams[i];
	  }
	  m_CphiMX_ij[counter] = vParams[0];
	}
      }

      if (nodeName == "alpha1") {
	stemp = xmlChild.value();
	m_Alpha1MX_ij[counter] = atofCheck(stemp.c_str());
      }
    }
  }

  /**
   * Process an XML node called "ThetaAnion". 
   * This node contains all of the parameters necessary to describe
   * the binary interactions between two anions.
   */
  void HMWSoln::readXMLThetaAnion(XML_Node &BinSalt) {
    double *charge = DATA_PTR(m_speciesCharge);
    string stemp;
    string iName = BinSalt.attrib("anion1");
    if (iName == "") {
      throw CanteraError("HMWSoln::readXMLThetaAnion", "no anion1 attrib");
    }
    string jName = BinSalt.attrib("anion2");
    if (jName == "") {
      throw CanteraError("HMWSoln::readXMLThetaAnion", "no anion2 attrib");
    }
    /*
     * Find the index of the species in the current phase. It's not
     * an error to not find the species
     */
    int iSpecies = speciesIndex(iName);
    if (iSpecies < 0) {
      return;
    }
    if (charge[iSpecies] >= 0) {
      throw CanteraError("HMWSoln::readXMLThetaAnion", "anion1 charge problem");
    }
    int jSpecies = speciesIndex(jName);
    if (jSpecies < 0) {
      return;
    }
    if (charge[jSpecies] >= 0) {
      throw CanteraError("HMWSoln::readXMLThetaAnion", "anion2 charge problem");
    }
	
    int n = iSpecies * m_kk + jSpecies;
    int counter = m_CounterIJ[n];
    int num = BinSalt.nChildren();
    for (int i = 0; i < num; i++) {
      XML_Node &xmlChild = BinSalt.child(i);
      stemp = xmlChild.name();
      string nodeName = lowercase(stemp);
      if (nodeName == "theta") {
	stemp = xmlChild.value();
	double old = m_Theta_ij[counter];
	m_Theta_ij[counter] = atofCheck(stemp.c_str());
	if (old != 0.0) {
	  if (old != m_Theta_ij[counter]) {
	    throw CanteraError("HMWSoln::readXMLThetaAnion", "conflicting values");
	  }
	}
      }
    }
  } 

  /**
   * Process an XML node called "ThetaCation". 
   * This node contains all of the parameters necessary to describe
   * the binary interactions between two cation.
   */
  void HMWSoln::readXMLThetaCation(XML_Node &BinSalt) {
    double *charge = DATA_PTR(m_speciesCharge);
    string stemp;
    string iName = BinSalt.attrib("cation1");
    if (iName == "") {
      throw CanteraError("HMWSoln::readXMLThetaCation", "no cation1 attrib");
    }
    string jName = BinSalt.attrib("cation2");
    if (jName == "") {
      throw CanteraError("HMWSoln::readXMLThetaCation", "no cation2 attrib");
    }
    /*
     * Find the index of the species in the current phase. It's not
     * an error to not find the species
     */
    int iSpecies = speciesIndex(iName);
    if (iSpecies < 0) {
      return;
    }
    if (charge[iSpecies] <= 0) {
      throw CanteraError("HMWSoln::readXMLThetaCation", "cation1 charge problem");
    }
    int jSpecies = speciesIndex(jName);
    if (jSpecies < 0) {
      return;
    }
    if (charge[jSpecies] <= 0) {
      throw CanteraError("HMWSoln::readXMLThetaCation", "cation2 charge problem");
    }
	
    int n = iSpecies * m_kk + jSpecies;
    int counter = m_CounterIJ[n];
    int num = BinSalt.nChildren();
    for (int i = 0; i < num; i++) {
      XML_Node &xmlChild = BinSalt.child(i);
      stemp = xmlChild.name();
      string nodeName = lowercase(stemp);
      if (nodeName == "theta") {
	stemp = xmlChild.value();
	double old = m_Theta_ij[counter];
	m_Theta_ij[counter] = atofCheck(stemp.c_str());
	if (old != 0.0) {
	  if (old != m_Theta_ij[counter]) {
	    throw CanteraError("HMWSoln::readXMLThetaCation", "conflicting values");
	  }
	}
      }
    }
  }

  /**
   * Process an XML node called "readXMLPsiCommonCation". 
   * This node contains all of the parameters necessary to describe
   * the binary interactions between two anions and one common cation.
   */
  void HMWSoln::readXMLPsiCommonCation(XML_Node &BinSalt) {
    double *charge = DATA_PTR(m_speciesCharge);
    string stemp;
    string kName = BinSalt.attrib("cation");
    if (kName == "") {
      throw CanteraError("HMWSoln::readXMLPsiCommonCation", "no cation attrib");
    }
    string iName = BinSalt.attrib("anion1");
    if (iName == "") {
      throw CanteraError("HMWSoln::readXMLPsiCommonCation", "no anion1 attrib");
    }
    string jName = BinSalt.attrib("anion2");
    if (jName == "") {
      throw CanteraError("HMWSoln::readXMLPsiCommonCation", "no anion2 attrib");
    }
    /*
     * Find the index of the species in the current phase. It's not
     * an error to not find the species
     */
    int kSpecies = speciesIndex(kName);
    if (kSpecies < 0) {
      return;
    }
    if (charge[kSpecies] <= 0) {
      throw CanteraError("HMWSoln::readXMLPsiCommonCation", 
			 "cation charge problem");
    }
    int iSpecies = speciesIndex(iName);
    if (iSpecies < 0) {
      return;
    }
    if (charge[iSpecies] >= 0) {
      throw CanteraError("HMWSoln::readXMLPsiCommonCation", 
			 "anion1 charge problem");
    }
    int jSpecies = speciesIndex(jName);
    if (jSpecies < 0) {
      return;
    }
    if (charge[jSpecies] >= 0) {
      throw CanteraError("HMWSoln::readXMLPsiCommonCation", 
			 "anion2 charge problem");
    }
	
    int n = iSpecies * m_kk + jSpecies;
    int counter = m_CounterIJ[n];
    int num = BinSalt.nChildren();
    for (int i = 0; i < num; i++) {
      XML_Node &xmlChild = BinSalt.child(i);
      stemp = xmlChild.name();
      string nodeName = lowercase(stemp);
      if (nodeName == "theta") {
	stemp = xmlChild.value();
	double old = m_Theta_ij[counter];
	m_Theta_ij[counter] = atofCheck(stemp.c_str());
	if (old != 0.0) {
	  if (old != m_Theta_ij[counter]) {
	    throw CanteraError("HMWSoln::readXMLPsiCommonCation",
			       "conflicting values");
	  }
	}
      }
      if (nodeName == "psi") {
	stemp = xmlChild.value();
	double param = atofCheck(stemp.c_str());
	n = iSpecies * m_kk *m_kk + jSpecies * m_kk + kSpecies ;
	m_Psi_ijk[n] = param;
	n = iSpecies * m_kk *m_kk + kSpecies * m_kk + jSpecies ;
	m_Psi_ijk[n] = param;
	n = jSpecies * m_kk *m_kk + iSpecies * m_kk + kSpecies ;
	m_Psi_ijk[n] = param;
	n = jSpecies * m_kk *m_kk + kSpecies * m_kk + iSpecies ;
	m_Psi_ijk[n] = param;
	n = kSpecies * m_kk *m_kk + jSpecies * m_kk + iSpecies ;
	m_Psi_ijk[n] = param;
	n = kSpecies * m_kk *m_kk + iSpecies * m_kk + jSpecies ;
	m_Psi_ijk[n] = param;
      }
    }
  }


   
  /**
   * Process an XML node called "PsiCommonAnion". 
   * This node contains all of the parameters necessary to describe
   * the binary interactions between two cations and one common anion.
   */
  void HMWSoln::readXMLPsiCommonAnion(XML_Node &BinSalt) {
    double *charge = DATA_PTR(m_speciesCharge);
    string stemp;
    string kName = BinSalt.attrib("anion");
    if (kName == "") {
      throw CanteraError("HMWSoln::readXMLPsiCommonAnion", "no anion attrib");
    }
    string iName = BinSalt.attrib("cation1");
    if (iName == "") {
      throw CanteraError("HMWSoln::readXMLPsiCommonAnion", "no cation1 attrib");
    }
    string jName = BinSalt.attrib("cation2");
    if (jName == "") {
      throw CanteraError("HMWSoln::readXMLPsiCommonAnion", "no cation2 attrib");
    }
    /*
     * Find the index of the species in the current phase. It's not
     * an error to not find the species
     */
    int kSpecies = speciesIndex(kName);
    if (kSpecies < 0) {
      return;
    }
    if (charge[kSpecies] >= 0) {
      throw CanteraError("HMWSoln::readXMLPsiCommonAnion", "anion charge problem");
    }
    int iSpecies = speciesIndex(iName);
    if (iSpecies < 0) {
      return;
    }
    if (charge[iSpecies] <= 0) {
      throw CanteraError("HMWSoln::readXMLPsiCommonAnion", 
			 "cation1 charge problem");
    }
    int jSpecies = speciesIndex(jName);
    if (jSpecies < 0) {
      return;
    }
    if (charge[jSpecies] <= 0) {
      throw CanteraError("HMWSoln::readXMLPsiCommonAnion", 
			 "cation2 charge problem");
    }
	
    int n = iSpecies * m_kk + jSpecies;
    int counter = m_CounterIJ[n];
    int num = BinSalt.nChildren();
    for (int i = 0; i < num; i++) {
      XML_Node &xmlChild = BinSalt.child(i);
      stemp = xmlChild.name();
      string nodeName = lowercase(stemp);
      if (nodeName == "theta") {
	stemp = xmlChild.value();
	double old = m_Theta_ij[counter];
	m_Theta_ij[counter] = atofCheck(stemp.c_str());
	if (old != 0.0) {
	  if (old != m_Theta_ij[counter]) {
	    throw CanteraError("HMWSoln::readXMLPsiCommonAnion",
			       "conflicting values");
	  }
	}
      }
      if (nodeName == "psi") {
	stemp = xmlChild.value();
	double param = atofCheck(stemp.c_str());
	n = iSpecies * m_kk *m_kk + jSpecies * m_kk + kSpecies ;
	m_Psi_ijk[n] = param;
	n = iSpecies * m_kk *m_kk + kSpecies * m_kk + jSpecies ;
	m_Psi_ijk[n] = param;
	n = jSpecies * m_kk *m_kk + iSpecies * m_kk + kSpecies ;
	m_Psi_ijk[n] = param;
	n = jSpecies * m_kk *m_kk + kSpecies * m_kk + iSpecies ;
	m_Psi_ijk[n] = param;
	n = kSpecies * m_kk *m_kk + jSpecies * m_kk + iSpecies ;
	m_Psi_ijk[n] = param;
	n = kSpecies * m_kk *m_kk + iSpecies * m_kk + jSpecies ;
	m_Psi_ijk[n] = param;
      }
    }
  }
   
  /**
   * Process an XML node called "LambdaNeutral". 
   * This node contains all of the parameters necessary to describe
   * the binary interactions between one neutral species and
   * any other species (neutral or otherwise) in the mechanism.
   */
  void HMWSoln::readXMLLambdaNeutral(XML_Node &BinSalt) {
    double *charge = DATA_PTR(m_speciesCharge);
    string stemp;
    string iName = BinSalt.attrib("neutral");
    if (iName == "") {
      throw CanteraError("HMWSoln::readXMLLambdaNeutral", "no neutral attrib");
    }
    string jName = BinSalt.attrib("speciesj");
    if (jName == "") {
      throw CanteraError("HMWSoln::readXMLLambdaNeutral", "no speciesj attrib");
    }
    /*
     * Find the index of the species in the current phase. It's not
     * an error to not find the species
     */
    int iSpecies = speciesIndex(iName);
    if (iSpecies < 0) {
      return;
    }
    if (charge[iSpecies] != 0) {
      throw CanteraError("HMWSoln::readXMLLambdaNeutral", "neutral charge problem");
    }
    int jSpecies = speciesIndex(jName);
    if (jSpecies < 0) {
      return;
    }
	
    int num = BinSalt.nChildren();
    for (int i = 0; i < num; i++) {
      XML_Node &xmlChild = BinSalt.child(i);
      stemp = xmlChild.name();
      string nodeName = lowercase(stemp);
      if (nodeName == "lambda") {
	stemp = xmlChild.value();
	double old = m_Lambda_ij(iSpecies,jSpecies);
	m_Lambda_ij(iSpecies,jSpecies) = atofCheck(stemp.c_str());
	if (old != 0.0) {
	  if (old != m_Lambda_ij(iSpecies,jSpecies)) {
	    throw CanteraError("HMWSoln::readXMLLambdaNeutral", "conflicting values");
	  }
	}
      }
    }
  }

  /**
   *  Initialization routine for a HMWSoln phase.
   *
   * This is a virtual routine. This routine will call initThermo()
   * for the parent class as well.
   */
  void HMWSoln::initThermo() {
    MolalityVPSSTP::initThermo();
    initLengths();
  }

  /**
   *   Import, construct, and initialize a HMWSoln phase 
   *   specification from an XML tree into the current object.
   *
   * This routine is a precursor to constructPhaseXML(XML_Node*)
   * routine, which does most of the work.
   *
   * @param infile XML file containing the description of the
   *        phase
   *
   * @param id  Optional parameter identifying the name of the
   *            phase. If none is given, the first XML
   *            phase element will be used.
   */
  void HMWSoln::constructPhaseFile(string inputFile, string id) {

    if (inputFile.size() == 0) {
      throw CanteraError("HMWSoln:constructPhaseFile",
			 "input file is null");
    }
    string path = findInputFile(inputFile);
    ifstream fin(path.c_str());
    if (!fin) {
      throw CanteraError("HMWSoln:constructPhaseFile","could not open "
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
      throw CanteraError("HMWSoln:constructPhaseFile",
			 "ERROR: Can not find phase named " +
			 id + " in file named " + inputFile);
    }
    fxml_phase->copy(&phaseNode_XML);	
    constructPhaseXML(*fxml_phase, id);
    delete fxml;
  }
  
  /**
   *   Import, construct, and initialize a HMWSoln phase 
   *   specification from an XML tree into the current object.
   *
   *   Most of the work is carried out by the cantera base
   *   routine, importPhase(). That routine imports all of the
   *   species and element data, including the standard states
   *   of the species.
   *
   *   Then, In this routine, we read the information 
   *   particular to the specification of the activity 
   *   coefficient model for the Pitzer parameterization.
   *
   *   We also read information about the molar volumes of the
   *   standard states if present in the XML file.
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
  void HMWSoln::constructPhaseXML(XML_Node& phaseNode, string id) {
    string stemp;
    if (id.size() > 0) {
      string idp = phaseNode.id();
      if (idp != id) {
	throw CanteraError("HMWSoln::constructPhaseXML", 
			   "phasenode and Id are incompatible");
      }
    }

    /*
     * Find the Thermo XML node 
     */
    if (!phaseNode.hasChild("thermo")) {
      throw CanteraError("HMWSoln::constructPhaseXML",
			 "no thermo XML node");
    }
    XML_Node& thermoNode = phaseNode.child("thermo");

    /*
     * Possibly change the form of the standard concentrations
     */
    if (thermoNode.hasChild("standardConc")) {
      XML_Node& scNode = thermoNode.child("standardConc");
      m_formGC = 2;
      stemp = scNode.attrib("model");
      string formString = lowercase(stemp);
      if (formString != "") {
	if (formString == "unity") {
	  m_formGC = 0;
	  printf("exit standardConc = unity not done\n");
	  exit(-1);
	} else if (formString == "molar_volume") {
	  m_formGC = 1;
	  printf("exit standardConc = molar_volume not done\n");
	  exit(-1);
	} else if (formString == "solvent_volume") {
	  m_formGC = 2;
	} else {
	  throw CanteraError("HMWSoln::constructPhaseXML",
			     "Unknown standardConc model: " + formString);
	}
      }
    }
    /*
     * Get the Name of the Solvent:
     *      <solvent> solventName </solvent>
     */
    string solventName = "";
    if (thermoNode.hasChild("solvent")) {
      XML_Node& scNode = thermoNode.child("solvent");
      vector<string> nameSolventa;
      getStringArray(scNode, nameSolventa);
      int nsp = static_cast<int>(nameSolventa.size());
      if (nsp != 1) {
	throw CanteraError("HMWSoln::constructPhaseXML",
			   "badly formed solvent XML node");
      }
      solventName = nameSolventa[0];
    }

    /*
     * Determine the form of the Pitzer model,
     *   We will use this information to size arrays below.
     */
    if (thermoNode.hasChild("activityCoefficients")) {
      XML_Node& scNode = thermoNode.child("activityCoefficients");
      m_formPitzer = m_formPitzer;
      stemp = scNode.attrib("model");
      string formString = lowercase(stemp);
      if (formString != "") {
	if        (formString == "pitzer" || formString == "default") {
	  m_formPitzer = PITZERFORM_BASE;
	} else if (formString == "base") {
	  m_formPitzer = PITZERFORM_BASE;
	} else {
	  throw CanteraError("HMWSoln::constructPhaseXML",
			     "Unknown Pitzer ActivityCoeff model: "
			     + formString);
	}
      }

      /*
       * Determine the form of the temperature dependence
       * of the Pitzer activity coefficient model.
       */
      stemp = scNode.attrib("TempModel");
      formString = lowercase(stemp);
      if (formString != "") {
	if (formString == "constant" || formString == "default") {
	  m_formPitzerTemp = PITZER_TEMP_CONSTANT;
	} else if (formString == "linear") {
	  m_formPitzerTemp = PITZER_TEMP_LINEAR;
	} else if (formString == "complex" || formString == "complex1") {
	  m_formPitzerTemp = PITZER_TEMP_COMPLEX1;
	} else {
	  throw CanteraError("HMWSoln::constructPhaseXML",
			     "Unknown Pitzer ActivityCoeff Temp model: "
			     + formString);
	}
      }

      /*
       * Determine the reference temperature
       * of the Pitzer activity coefficient model's temperature
       * dependence formulation: defaults to 25C
       */
      stemp = scNode.attrib("TempReference");
      formString = lowercase(stemp);
      if (formString != "") {
	m_TempPitzerRef = atofCheck(formString.c_str());
      } else {
	m_TempPitzerRef = 273.15 + 25;
      }
	
    }

    /*
     * Call the Cantera importPhase() function. This will import
     * all of the species into the phase. This will also handle
     * all of the solvent and solute standard states
     */
    bool m_ok = importPhase(phaseNode, this);
    if (!m_ok) {
      throw CanteraError("HMWSoln::constructPhaseXML","importPhase failed "); 
    }
    
  }

  /**
   * Process the XML file after species are set up.
   *
   *  This gets called from importPhase(). It processes the XML file
   *  after the species are set up. This is the main routine for
   *  reading in activity coefficient parameters.
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
  void HMWSoln::
  initThermoXML(XML_Node& phaseNode, string id) {
    int k;
    string stemp;
    /*
     * Find the Thermo XML node 
     */
    if (!phaseNode.hasChild("thermo")) {
      throw CanteraError("HMWSoln::initThermoXML",
			 "no thermo XML node");
    }
    XML_Node& thermoNode = phaseNode.child("thermo");

    /*
     * Get the Name of the Solvent:
     *      <solvent> solventName </solvent>
     */
    string solventName = "";
    if (thermoNode.hasChild("solvent")) {
      XML_Node& scNode = thermoNode.child("solvent");
      vector<string> nameSolventa;
      getStringArray(scNode, nameSolventa);
      int nsp = static_cast<int>(nameSolventa.size());
      if (nsp != 1) {
	throw CanteraError("HMWSoln::initThermoXML",
			   "badly formed solvent XML node");
      }
      solventName = nameSolventa[0];
    }

    /*
     * Initialize all of the lengths of arrays in the object
     * now that we know what species are in the phase.
     */
    initLengths();

    /*
     * Reconcile the solvent name and index.
     */
    for (k = 0; k < m_kk; k++) {
      string sname = speciesName(k);
      if (solventName == sname) {
	setSolvent(k);
	if (k != 0) {
	  throw CanteraError("HMWSoln::initThermoXML",
			     "Solvent must be species 0 atm");
	}
	m_indexSolvent = k;
	break;
      }
    }
    if (m_indexSolvent == -1) {
      cout << "HMWSoln::initThermo: Solvent Name not found" 
	   << endl;
      throw CanteraError("HMWSoln::initThermoXML",
			 "Solvent name not found");
    }
    if (m_indexSolvent != 0) {
      throw CanteraError("HMWSoln::initThermoXML",
			 "Solvent " + solventName +
			 " should be first species");
    }

    /*
     * Now go get the specification of the standard states for
     * species in the solution. This includes the molar volumes
     * data blocks for incompressible species.
     */
    XML_Node& speciesList = phaseNode.child("speciesArray");
    XML_Node* speciesDB =
      get_XML_NameID("speciesData", speciesList["datasrc"],
		     &phaseNode.root());
    const vector<string>&sss = speciesNames();

    for (k = 0; k < m_kk; k++) {
      XML_Node* s =  speciesDB->findByAttr("name", sss[k]);
      if (!s) {
	throw CanteraError("HMWSoln::initThermoXML",
			 "Species Data Base " + sss[k] + " not found");
      }
      XML_Node *ss = s->findByName("standardState");
      if (!ss) {
	throw CanteraError("HMWSoln::initThermoXML",
			 "Species " + sss[k] + 
			   " standardState XML block  not found");
      }
      string modelStringa = ss->attrib("model");
      if (modelStringa == "") {
	throw CanteraError("HMWSoln::initThermoXML",
			   "Species " + sss[k] + 
			   " standardState XML block model attribute not found");
      }
      string modelString = lowercase(modelStringa);
      if (k == 0) {
	if (modelString == "wateriapws" || modelString == "real_water" ||
	    modelString == "waterpdss") {
	  /*
	   * Initialize the water standard state model
	   */
	  if (m_waterSS) delete m_waterSS;
	  m_waterSS = new WaterPDSS(this, 0);
	  /*
	   * Fill in the molar volume of water (m3/kmol)
	   * at standard conditions to fill in the m_speciesSize entry
	   * with something reasonable.
	   */
	  m_waterSS->setState_TP(300., OneAtm);
	  double dens = m_waterSS->density();
	  double mw = m_waterSS->molecularWeight();
	  m_speciesSize[0] = mw / dens;
#ifdef DEBUG_HKM_NOT
	  cout << "Solvent species " << sss[k] << " has volume " <<  
	    m_speciesSize[k] << endl;
#endif
	} else {
	  throw CanteraError("HMWSoln::initThermoXML",
			     "Solvent SS Model \"" + modelStringa + 
			     "\" is not allowed");
	}
      } else {
	if (modelString != "constant_incompressible") {
	  throw CanteraError("HMWSoln::initThermoXML",
			     "Solute SS Model \"" + modelStringa + 
			     "\" is not known");
	}
	m_speciesSize[k] = getFloat(*ss, "molarVolume", "-");
#ifdef DEBUG_HKM_NOT
	cout << "species " << sss[k] << " has volume " <<  
   	  m_speciesSize[k] << endl;
#endif
      }
    }

    /*
     * Initialize the water property calculator. It will share
     * the internal eos water calculator.
     */
    m_waterProps = new WaterProps(m_waterSS);

    /*
     * Go get all of the coefficients and factors in the
     * activityCoefficients XML block
     */
    XML_Node *acNodePtr = 0;
    if (thermoNode.hasChild("activityCoefficients")) {
      XML_Node& acNode = thermoNode.child("activityCoefficients");
      acNodePtr = &acNode;
      /*
       * Look for parameters for A_Debye
       */
      if (acNode.hasChild("A_Debye")) {
	XML_Node &ADebye = acNode.child("A_Debye");
	m_form_A_Debye = A_DEBYE_CONST;
	stemp = "model";
	if (ADebye.hasAttrib(stemp)) {
	  string atemp = ADebye.attrib(stemp);
	  stemp = lowercase(atemp);
	  if (stemp == "water") {
	    m_form_A_Debye = A_DEBYE_WATER;
	  }
	}
	if (m_form_A_Debye == A_DEBYE_CONST) {
	  m_A_Debye = getFloat(acNode, "A_Debye");
	}
#ifdef DEBUG_HKM_NOT
	cout << "A_Debye = " << m_A_Debye << endl;
#endif
      }

      /*
       * Look for Parameters for the Maximum Ionic Strength
       */
      if (acNode.hasChild("maxIonicStrength")) {
	m_maxIionicStrength = getFloat(acNode, "maxIonicStrength");
#ifdef DEBUG_HKM_NOT
	cout << "m_maxIionicStrength = "
	     <<m_maxIionicStrength << endl;
#endif
      }


      /*
       * Look for parameters for the Ionic radius
       */
      if (acNode.hasChild("ionicRadius")) {
	XML_Node& irNode = acNode.child("ionicRadius");

	string Aunits = "";
	double Afactor = 1.0;
	if (irNode.hasAttrib("units")) {
	  string Aunits = irNode.attrib("units");
	  Afactor = toSI(Aunits); 
	}

	if (irNode.hasAttrib("default")) {
	  string ads = irNode.attrib("default");
	  double ad = fpValue(ads);
	  for (int k = 0; k < m_kk; k++) {
	    m_Aionic[k] = ad * Afactor;
	  }
	}

      }


      /*
       * Fill in parameters for the calculation of the 
       * stoichiometric Ionic Strength
       *
       * The default is that stoich charge is the same as the
       * regular charge.
       */
      m_speciesCharge_Stoich.resize(m_kk, 0.0);
      for (k = 0; k < m_kk; k++) {
	m_speciesCharge_Stoich[k] = m_speciesCharge[k];
      }
      /*
       * First look at the species database.
       *  -> Look for the subelement "stoichIsMods"
       *     in each of the species SS databases.
       */
      const XML_Node *phaseSpecies = speciesData();
      if (phaseSpecies) {
	string kname, jname;
	vector<XML_Node*> xspecies;
	phaseSpecies->getChildren("species", xspecies);
	int jj = xspecies.size();
	for (k = 0; k < m_kk; k++) {
	  int jmap = -1;
	  kname = speciesName(k);
	  for (int j = 0; j < jj; j++) {
	    const XML_Node& sp = *xspecies[j];
	    jname = sp["name"];
	    if (jname == kname) {
	      jmap = j;
	      break;
	    }
	  }
	  if (jmap > -1) {
	    const XML_Node& sp = *xspecies[jmap];
	    if (sp.hasChild("stoichIsMods")) {
	      double val = getFloat(sp, "stoichIsMods");
	      m_speciesCharge_Stoich[k] = val;
	    }
	  }
	}
      }
      /*
       * Now look at the activity coefficient database
       */
      if (acNodePtr) {
	if (acNodePtr->hasChild("stoichIsMods")) {
	  XML_Node& sIsNode = acNodePtr->child("stoichIsMods");

	  map<string, string> msIs;
	  getMap(sIsNode, msIs);
	  map<string,string>::const_iterator _b = msIs.begin();
	  for (; _b != msIs.end(); ++_b) {
	    int kk = speciesIndex(_b->first);
	    if (kk < 0) {
	      //throw CanteraError(
	      //   "HMWSoln::initThermo error",
	      //   "no species match was found"
	      //   );
	    } else {
	      double val = fpValue(_b->second);
	      m_speciesCharge_Stoich[kk] = val;
	    }
	  }
	}
      }

      /*
       * Loop through the children getting multiple instances of
       * parameters
       */
      if (acNodePtr) {
	int n = acNodePtr->nChildren();
	for (int i = 0; i < n; i++) {
	  XML_Node &xmlACChild = acNodePtr->child(i);
	  stemp = xmlACChild.name();
	  string nodeName = lowercase(stemp);
	  /*
	   * Process a binary salt field, or any of the other XML fields
	   * that make up the Pitzer Database. Entries will be ignored
	   * if any of the species in the entry isn't in the solution.
	   */
	  if (nodeName == "binarysaltparameters") {
	    readXMLBinarySalt(xmlACChild);
	  } else if (nodeName == "thetaanion") {
	    readXMLThetaAnion(xmlACChild);
	  } else if (nodeName == "thetacation") {
	    readXMLThetaCation(xmlACChild);
	  } else if (nodeName == "psicommonanion") { 
	    readXMLPsiCommonAnion(xmlACChild);
	  } else if (nodeName == "psicommoncation") { 
	    readXMLPsiCommonCation(xmlACChild);
	  } else if (nodeName == "lambdaneutral") {
	    readXMLLambdaNeutral(xmlACChild);
	  }
	}
      }


    }

    /*
     * Fill in the vector specifying the electrolyte species
     * type
     *
     *   First fill in default values. Everthing is either
     *   a charge species, a nonpolar neutral, or the solvent.
     */
    for (k = 0; k < m_kk; k++) {
      if (fabs(m_speciesCharge[k]) > 0.0001) {
	m_electrolyteSpeciesType[k] = cEST_chargedSpecies;
	if (fabs(m_speciesCharge_Stoich[k] - m_speciesCharge[k])
	    > 0.0001) {
	  m_electrolyteSpeciesType[k] = cEST_weakAcidAssociated;
	}
      } else if (fabs(m_speciesCharge_Stoich[k]) > 0.0001) {
	m_electrolyteSpeciesType[k] = cEST_weakAcidAssociated;
      } else {
	m_electrolyteSpeciesType[k] = cEST_nonpolarNeutral;
      }
    }
    m_electrolyteSpeciesType[m_indexSolvent] = cEST_solvent;
    /*
     * First look at the species database.
     *  -> Look for the subelement "stoichIsMods"
     *     in each of the species SS databases.
     */
    const XML_Node *phaseSpecies = speciesData();
    const XML_Node *spPtr = 0;
    if (phaseSpecies) {
      string kname;
      for (k = 0; k < m_kk; k++) {
	kname = speciesName(k);
	spPtr = speciesXML_Node(kname, phaseSpecies);
	if (!spPtr) {
	  if (spPtr->hasChild("electrolyteSpeciesType")) {
	    string est = getString(*spPtr, "electrolyteSpeciesType");
	    if ((m_electrolyteSpeciesType[k] = interp_est(est)) == -1) {
	      throw CanteraError("HMWSoln::initThermoXML",
				 "Bad electrolyte type: " + est);
	    }
	  }
	}
      }
    }
    /*
     * Then look at the phase thermo specification
     */
    if (acNodePtr) {
      if (acNodePtr->hasChild("electrolyteSpeciesType")) {
	XML_Node& ESTNode = acNodePtr->child("electrolyteSpeciesType");
	map<string, string> msEST;
	getMap(ESTNode, msEST);
	map<string,string>::const_iterator _b = msEST.begin();
	for (; _b != msEST.end(); ++_b) {
	  int kk = speciesIndex(_b->first);
	  if (kk < 0) {
	  } else {
	    string est = _b->second;
	    if ((m_electrolyteSpeciesType[kk] = interp_est(est))  == -1) {
	      throw CanteraError("HMWSoln::initThermoXML",
				 "Bad electrolyte type: " + est);
	    }
	  }
	}
      }
    }

    /*
     * Lastly set the state
     */
    if (phaseNode.hasChild("state")) {
      XML_Node& stateNode = phaseNode.child("state");
      setStateFromXML(stateNode);
    }

  }
}
