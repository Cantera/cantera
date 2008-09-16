/**
 *  @file VPSSMgrFactory.cpp
 *    Definitions for factory to build instances of classes that manage the
 *    calculation of standard state properties for all the species in a phase
 *    (see \ref spthermo and class 
 *    \link Cantera::VPSSMgrFactory VPSSMgrFactory\endlink);
 */
/*
 * $Id$
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#ifdef WIN32
#pragma warning(disable:4786)
#endif


#include "SpeciesThermo.h"


#include "VPSSMgr.h"
#include "VPSSMgrFactory.h"

#include "VPStandardStateTP.h"

#include "VPSSMgr_IdealGas.h"
#include "VPSSMgr_ConstVol.h"
#include "VPSSMgr_Water_ConstVol.h"
#include "VPSSMgr_Water_HKFT.h"
#include "VPSSMgr_General.h"

#include "VPSSMgr_types.h"

#include "SpeciesThermoMgr.h"
#include "speciesThermoTypes.h"
#include "SpeciesThermo.h"
#include "SpeciesThermoFactory.h"
#include "GeneralSpeciesThermo.h"

#include "mix_defs.h"

#include "xml.h"
#include "ctml.h"

using namespace ctml;
using namespace std;


namespace Cantera {

  VPSSMgrFactory* VPSSMgrFactory::s_factory = 0;

#if defined(THREAD_SAFE_CANTERA)
  // Defn of the static mutex variable that locks the 
  // %VPSSMgr factory singelton
  boost::mutex VPSSMgrFactory::vpss_species_thermo_mutex;
#endif
 
  /*
   * Examine the types of species thermo parameterizations,
   * and return a flag indicating the type of parameterization
   * needed by the species.
   * 
   *  @param spData_node Species Data XML node. This node contains a list
   *                     of species XML nodes underneath it.
   * 
   * @todo Make sure that spDadta_node is species Data XML node by checking its name is speciesData
   */
  static void getVPSSMgrTypes(XML_Node* spData_node, 
			      int& has_nasa, int& has_shomate, int& has_simple,
			      int &has_water,
			      int &has_tpx,
			      int &has_hptx,
			      int &has_other) {

    const XML_Node& sparray = *spData_node;
    std::vector<XML_Node*> sp;

    // get all of the species nodes
    sparray.getChildren("species",sp);
    size_t n, ns = sp.size();
    for (n = 0; n < ns; n++) {
      bool ifound = false;
      XML_Node* spNode = sp[n];
      if (spNode->hasChild("standardState")) {
	const XML_Node& ssN = sp[n]->child("standardState");
	string mm = ssN["model"];
	if (mm == "waterIAPWS" || mm == "waterPDSS") {
	  has_water++;
	  ifound = true;
	}
	if (mm == "HKFT") {
	  has_hptx++;
	  ifound = true;
	}
      } 
      if (!ifound) {
	if (spNode->hasChild("thermo")) {
	  const XML_Node& th = sp[n]->child("thermo");
	  if (th.hasChild("NASA")) {
	    has_nasa++;
	    ifound = true;
	  }
	  if (th.hasChild("Shomate")) {
	    has_shomate++;
	    ifound = true;
	  }
	  if (th.hasChild("const_cp")){
	    has_simple = 1;
	    ifound = true;
	  }
	  if (th.hasChild("poly")) {
	    if (th.child("poly")["order"] == "1") {
	      has_simple = 1;
	      ifound = true;
	    } else throw CanteraError("newSpeciesThermo",
				      "poly with order > 1 not yet supported");
	  }
	  if (th.hasChild("Mu0")) {
	    has_other++;
	    ifound = true;
	  }
	  if (th.hasChild("NASA9")) {
	    has_other++;
	    ifound = true;
	  }
	  if (th.hasChild("NASA9MULTITEMP")) {
	    has_other++;
	    ifound = true;
	  }
	  if (th.hasChild("adsorbate")) {
	    has_other++;
	    ifound = true;
	  }
	  if (th.hasChild("HKFT")) {
	    has_hptx++;
	    ifound = true;
	  }
	} else {
	  throw UnknownVPSSMgrModel("getVPSSMgrTypes:",
				    spNode->attrib("name"));
	}
      }
    }
  }

  VPSSMgr_enumType
  VPSSMgrFactory::VPSSMgr_StringConversion(std::string ssModel) const {
    VPSSMgr_enumType type;
    if (ssModel == "IdealGas") {
      type = cVPSSMGR_IDEALGAS;
    } else if (ssModel == "ConstVol") {
      type = cVPSSMGR_CONSTVOL;
    } else if (ssModel == "PureFuild") {
      type = cVPSSMGR_PUREFLUID;
    } else if (ssModel == "Water_ConstVol") {
      type = cVPSSMGR_WATER_CONSTVOL;
    } else if (ssModel == "Water_HKFT") {
      type = cVPSSMGR_WATER_HKFT;
    } else if (ssModel == "General") {
      type = cVPSSMGR_GENERAL;
    } else {
      type = cVPSSMGR_UNDEF;
    }
    return type;
  }

  // Stub out of new capabilities.
    
  VPSSMgr* 
  VPSSMgrFactory::newVPSSMgr(VPStandardStateTP *vp_ptr, 
			     XML_Node* phaseNode_ptr,
			     XML_Node* spData_node) {
    std::string ssModel="";
    VPSSMgr *vpss = 0;
    // First look for any explicit instructions within the XML Data
    if (phaseNode_ptr) {
      if (phaseNode_ptr->hasChild("thermo")) {
	const XML_Node& thermoNode = phaseNode_ptr->child("thermo");
	if (thermoNode.hasChild("standardState")) {
	  const XML_Node& ssNode = thermoNode.child("standardState");
	  ssModel = ssNode["model"];
	}
      }
    }
  

    // first get the reference state handler
    SpeciesThermo *spth = newSpeciesThermoMgr(spData_node);
    vp_ptr->setSpeciesThermo(spth);

    if (ssModel != "") {
      VPSSMgr_enumType type = VPSSMgr_StringConversion(ssModel);
      vpss = newVPSSMgr(type, vp_ptr);
      return vpss;
    }

    // If it comes back as general, then there may be some unknown 
    // parameterizations to the SpeciesThermo factory routine.
    bool haveSomeUnknowns = true;
    GeneralSpeciesThermo *ttmp = dynamic_cast<GeneralSpeciesThermo *>(spth);
    if (ttmp == 0) {
      haveSomeUnknowns = false;
    }

  
    if (vp_ptr->eosType() == cVPSS_IdealGas) {
      vpss = new VPSSMgr_IdealGas(vp_ptr, spth);

    }

    if (vp_ptr->eosType() == cVPSS_ConstVol) {
      vpss = new VPSSMgr_ConstVol(vp_ptr, spth);
    }

    int inasa = 0, ishomate = 0, isimple = 0, iwater = 0, itpx = 0, iother = 0;
    int ihptx = 0;
  
    try {
      getVPSSMgrTypes(spData_node, inasa, ishomate, isimple, iwater, 
		      itpx, ihptx, iother);
    } catch (UnknownSpeciesThermoModel) {
      iother = 1;
      popError();
    }
  
    if (iwater == 1) {
      if (ihptx == 0) {
	vpss = new VPSSMgr_Water_ConstVol(vp_ptr, spth);
      } else {
	vpss = new VPSSMgr_Water_HKFT(vp_ptr, spth);
      }
    }
    // The default here is to fall back to use the completely
    // general representation.
    if (vpss == 0) {
      vpss = new VPSSMgr_General(vp_ptr, spth);
    }
    return vpss;
  }

  VPSSMgr* 
  VPSSMgrFactory::newVPSSMgr(VPStandardStateTP *vp_ptr, 
			     XML_Node* phaseNode_ptr,
			     std::vector<XML_Node*> spData_nodes) {

    std::string ssModel="";
    VPSSMgr *vpss = 0;
    // First look for any explicit instructions within the XML Data
    if (phaseNode_ptr) {
      if (phaseNode_ptr->hasChild("thermo")) {
	const XML_Node& thermoNode = phaseNode_ptr->child("thermo");
	if (thermoNode.hasChild("standardState")) {
	  const XML_Node& ssNode = thermoNode.child("standardState");
	  ssModel = ssNode["model"];
	}
      }
    }

    // first get the reference state handler
    SpeciesThermo *spth = newSpeciesThermoMgr(spData_nodes);
    vp_ptr->setSpeciesThermo(spth);

    if (ssModel != "") {
      VPSSMgr_enumType type = VPSSMgr_StringConversion(ssModel);
      vpss = newVPSSMgr(type, vp_ptr);
      return vpss;
    }

    // If it comes back as general, then there may be some unknown 
    // parameterizations to the SpeciesThermo factory routine.
    bool haveSomeUnknowns = true;
    GeneralSpeciesThermo *ttmp = dynamic_cast<GeneralSpeciesThermo *>(spth);
    if (ttmp == 0) {
      haveSomeUnknowns = false;
    }

    if (vp_ptr->eosType() == cIdealSolnGasVPSS) {
      vpss = new VPSSMgr_IdealGas(vp_ptr, spth);
    }

    if (vp_ptr->eosType() == cIdealSolnGasVPSS_iscv) {
      vpss = new VPSSMgr_ConstVol(vp_ptr, spth);
    }

    int n = static_cast<int>(spData_nodes.size());
    int inasa = 0, ishomate = 0, isimple = 0, iwater = 0, itpx = 0, iother = 0;
    int ihptx = 0;
    for (int j = 0; j < n; j++) {
      try {
	getVPSSMgrTypes(spData_nodes[j], inasa, ishomate, isimple, iwater, 
				  itpx, ihptx, iother);
      } catch (UnknownSpeciesThermoModel) {
	iother = 1;
	popError();
      }
    }
    if (iwater == 1) {
      if (ihptx == 0) {
	vpss = new VPSSMgr_Water_ConstVol(vp_ptr, spth);
      } else {
	vpss = new VPSSMgr_Water_HKFT(vp_ptr, spth);
      }
    }
    if (vpss == 0) {
      vpss = new VPSSMgr_General(vp_ptr, spth);
    }
    return vpss;
  }

  

  // I don't think this is currently used. However, this is a virtual 
  // function where additional capabilities may be added.
  VPSSMgr* 
  VPSSMgrFactory::newVPSSMgr(VPSSMgr_enumType type, VPStandardStateTP *vp_ptr) {
    SpeciesThermo &spthermoRef = vp_ptr->speciesThermo();
    switch (type) {
    case cVPSSMGR_IDEALGAS:
      return new VPSSMgr_IdealGas(vp_ptr, &spthermoRef);
      break;
    case cVPSSMGR_CONSTVOL:
      return new VPSSMgr_ConstVol(vp_ptr, &spthermoRef);
      break;
    case cVPSSMGR_PUREFLUID:
      throw CanteraError("VPSSMgrFactory::newVPSSMgr",
			 "unimplemented");
    case cVPSSMGR_WATER_CONSTVOL:
      return new VPSSMgr_Water_ConstVol(vp_ptr, &spthermoRef);
      break;
    case cVPSSMGR_WATER_HKFT:
      return new VPSSMgr_Water_HKFT(vp_ptr, &spthermoRef);
      break;
    case cVPSSMGR_GENERAL:
      return new VPSSMgr_General(vp_ptr, &spthermoRef);
      break;
    case cVPSSMGR_UNDEF:
    default:
      throw UnknownVPSSMgrModel("VPSSMgrFactory::newVPSSMgr", int2str(type));
      return 0; 
    }
  }

  // I don't think this is currently used
  VPSSMgr* newVPSSMgr(VPSSMgr_enumType type, VPStandardStateTP *vp_ptr,
		      Cantera::VPSSMgrFactory* f) {
    if (f == 0) {
      f = VPSSMgrFactory::factory();
    }
    VPSSMgr* vpsssptherm = f->newVPSSMgr(type, vp_ptr);
    return vpsssptherm;
  }

  VPSSMgr* newVPSSMgr(VPStandardStateTP *tp_ptr,
		      XML_Node* phaseNode_ptr,
		      XML_Node* spData_node, 
		      VPSSMgrFactory* f) {
    if (f == 0) {
      f = VPSSMgrFactory::factory();
    }
    VPSSMgr* vpsssptherm = f->newVPSSMgr(tp_ptr, phaseNode_ptr, spData_node);
    return vpsssptherm;
  }

  VPSSMgr* newVPSSMgr(VPStandardStateTP *tp_ptr,
		      XML_Node* phaseNode_ptr,
		      std::vector<XML_Node*> spData_nodes, 
		      VPSSMgrFactory* f) {
    if (f == 0) {
      f = VPSSMgrFactory::factory();
    }
    VPSSMgr* vpsssptherm = f->newVPSSMgr(tp_ptr, phaseNode_ptr, spData_nodes);
    return vpsssptherm;
  }


}
