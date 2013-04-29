/**
 *  @file VPSSMgrFactory.cpp
 *    Definitions for factory to build instances of classes that manage the
 *    calculation of standard state properties for all the species in a phase
 *    (see \ref spthermo and class
 *    \link Cantera::VPSSMgrFactory VPSSMgrFactory\endlink);
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/thermo/SpeciesThermo.h"

#include "cantera/thermo/VPSSMgr.h"
#include "VPSSMgrFactory.h"

#include "cantera/thermo/VPStandardStateTP.h"

#include "cantera/thermo/VPSSMgr_IdealGas.h"
#include "cantera/thermo/VPSSMgr_ConstVol.h"
#include "cantera/thermo/VPSSMgr_Water_ConstVol.h"
#include "cantera/thermo/VPSSMgr_Water_HKFT.h"
#include "cantera/thermo/VPSSMgr_General.h"

#include "cantera/thermo/SpeciesThermoMgr.h"
#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/thermo/SpeciesThermo.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/thermo/GeneralSpeciesThermo.h"

#include "cantera/thermo/mix_defs.h"

#include "cantera/base/xml.h"
#include "cantera/base/ctml.h"

using namespace ctml;
using namespace std;

namespace Cantera
{

VPSSMgrFactory* VPSSMgrFactory::s_factory = 0;

// Defn of the static mutex variable that locks the VPSSMgr factory singleton
mutex_t VPSSMgrFactory::vpss_species_thermo_mutex;

//! Examine the types of species thermo parameterizations, and return a flag indicating the type of parameterization
//! needed by the species.
/*!
 *  @param spDataNodeList            Species Data XML node. This node contains a list
 *                                     of species XML nodes underneath it.
 *  @param has_nasa_idealGas         Boolean indicating that one species has a NASA ideal gas standard state
 *  @param has_nasa_constVol         Boolean indicating that one species has a NASA ideal solution standard state
 *  @param has_shomate_idealGas      Boolean indicating that one species has a shomate ideal gas standard state
 *  @param has_shomate_constVol      Boolean indicating that one species has a shomate ideal solution standard state
 *  @param has_simple_idealGas       Boolean indicating that one species has a simple ideal gas standard state
 *  @param has_simple_constVol       Boolean indicating that one species has a simple ideal solution standard state
 *  @param has_water                 Boolean indicating that one species has a water standard state
 *  @param has_tpx                   Boolean indicating that one species has a tpx standard state
 *  @param has_hptx                  Boolean indicating that one species has a htpx standard state
 *  @param has_other                 Boolean indicating that one species has different standard state than the ones listed above
 *
 * @todo Make sure that spDadta_node is species Data XML node by checking
 *      its name is speciesData
 */
static void getVPSSMgrTypes(std::vector<XML_Node*> & spDataNodeList,
                            int& has_nasa_idealGas,
                            int& has_nasa_constVol,
                            int& has_shomate_idealGas,
                            int& has_shomate_constVol,
                            int& has_simple_idealGas,
                            int& has_simple_constVol,
                            int& has_water,
                            int& has_tpx,
                            int& has_hptx,
                            int& has_other)
{

    XML_Node* ss_ptr = 0;
    string ssModel = "idealGas";
    size_t ns = spDataNodeList.size();
    for (size_t n = 0; n < ns; n++) {
        bool ifound = false;
        XML_Node* spNode = spDataNodeList[n];
        if (spNode->hasChild("standardState")) {
            const XML_Node& ssN = spNode->child("standardState");
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
                const XML_Node& th = spNode->child("thermo");
                if (spNode->hasChild("standardState")) {
                    ss_ptr = &(spNode->child("standardState"));
                    ssModel = ss_ptr->attrib("model");
                }
                if (th.hasChild("NASA")) {
                    if (ssModel == "idealGas") {
                        has_nasa_idealGas++;
                    } else if (ssModel == "constant_incompressible" ||
                               ssModel == "constantVolume") {
                        has_nasa_constVol++;
                    } else if (ssModel == "temperature_polynomial" ||
                               ssModel == "density_temperature_polynomial"  ||
                               ssModel == "constant") {
                        has_other++;
                    } else {
                        throw UnknownVPSSMgrModel("getVPSSMgrTypes:",
                                                  spNode->attrib("name"));
                    }
                    ifound = true;
                }
                if (th.hasChild("Shomate")) {
                    if (ssModel == "idealGas") {
                        has_shomate_idealGas++;
                    } else if (ssModel == "constant_incompressible" ||
                               ssModel == "constantVolume") {
                        has_shomate_constVol++;
                    } else if (ssModel == "temperature_polynomial" ||
                               ssModel == "density_temperature_polynomial"  ||
                               ssModel == "constant") {
                        has_other++;
                    } else {
                        throw UnknownVPSSMgrModel("getVPSSMgrTypes:",
                                                  spNode->attrib("name"));
                    }
                    ifound = true;
                }
                if (th.hasChild("const_cp")) {
                    if (ssModel == "idealGas") {
                        has_simple_idealGas++;
                    } else if (ssModel == "constant_incompressible" ||
                               ssModel == "constantVolume") {
                        has_simple_constVol++;
                    } else if (ssModel == "temperature_polynomial" ||
                               ssModel == "density_temperature_polynomial"  ||
                               ssModel == "constant") {
                        has_other++;
                    } else {
                        throw UnknownVPSSMgrModel("getVPSSMgrTypes:",
                                                  spNode->attrib("name"));
                    }
                    ifound = true;
                }
                if (th.hasChild("poly")) {
                    if (th.child("poly")["order"] == "1") {
                        has_simple_constVol = 1;
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

void VPSSMgrFactory::deleteFactory()
{
    ScopedLock lock(vpss_species_thermo_mutex);
    if (s_factory) {
        delete s_factory;
        s_factory = 0;
    }
}

VPSSMgr_enumType
VPSSMgrFactory::VPSSMgr_StringConversion(const std::string& ssModel) const
{
    std::string lssModel = lowercase(ssModel);
    VPSSMgr_enumType type;
    if (lssModel == "idealgas") {
        type = cVPSSMGR_IDEALGAS;
    } else if (lssModel == "constvol") {
        type = cVPSSMGR_CONSTVOL;
    } else if (lssModel == "purefuild") {
        type = cVPSSMGR_PUREFLUID;
    } else if (lssModel == "water_constvol") {
        type = cVPSSMGR_WATER_CONSTVOL;
    } else if (lssModel == "water_hkft") {
        type = cVPSSMGR_WATER_HKFT;
    } else if (lssModel == "general") {
        type = cVPSSMGR_GENERAL;
    } else {
        type = cVPSSMGR_UNDEF;
    }
    return type;
}

VPSSMgr*
VPSSMgrFactory::newVPSSMgr(VPStandardStateTP* vp_ptr,
                           XML_Node* phaseNode_ptr,
                           std::vector<XML_Node*> & spDataNodeList)
{

    std::string ssManager="";
    std::string vpssManager="";

    // First look for any explicit instructions within the XML Database
    // for the standard state manager and the variable pressure
    // standard state manager
    if (phaseNode_ptr) {
        if (phaseNode_ptr->hasChild("thermo")) {
            const XML_Node& thermoNode = phaseNode_ptr->child("thermo");
            if (thermoNode.hasChild("standardStateManager")) {
                const XML_Node& ssNode = thermoNode.child("standardStateManager");
                ssManager = ssNode["model"];
            }
            if (thermoNode.hasChild("variablePressureStandardStateManager")) {
                const XML_Node& vpssNode = thermoNode.child("variablePressureStandardStateManager");
                vpssManager = vpssNode["model"];
            }
        }
    }

    // first get the reference state handler. If we have explicit instructions,
    // use them to spawn the object.
    SpeciesThermo* spth = 0;
    if (ssManager != "") {
        spth = newSpeciesThermoMgr(ssManager);
    } else {
        spth = newSpeciesThermoMgr(spDataNodeList);
    }
    vp_ptr->setSpeciesThermo(spth);

    // Next, if we have specific directions, use them to get the VPSSSMgr object
    // and return immediately
    if (vpssManager != "") {
        VPSSMgr_enumType type = VPSSMgr_StringConversion(vpssManager);
        return newVPSSMgr(type, vp_ptr);
    }

    // Handle special cases based on the VPStandardState types
    if (vp_ptr->eosType() == cVPSS_IdealGas) {
        return new VPSSMgr_IdealGas(vp_ptr, spth);
    } else if (vp_ptr->eosType() == cVPSS_ConstVol) {
        return new VPSSMgr_ConstVol(vp_ptr, spth);
    }


    int inasaIG = 0, inasaCV = 0, ishomateIG = 0, ishomateCV = 0,
        isimpleIG = 0, isimpleCV = 0,
        iwater = 0, itpx = 0, iother = 0;
    int ihptx = 0;

    try {
        getVPSSMgrTypes(spDataNodeList, inasaIG, inasaCV, ishomateIG, ishomateCV,
                        isimpleIG, isimpleCV, iwater, itpx, ihptx, iother);
    } catch (UnknownSpeciesThermoModel) {
        iother = 1;
        popError();
    }

    if (iwater == 1) {
        if (ihptx == 0) {
            if (inasaIG ||  ishomateIG || isimpleIG) {
                throw CanteraError("newVPSSMgr", "Ideal gas with liquid water");
            } else {
                return new VPSSMgr_Water_ConstVol(vp_ptr, spth);
            }
        } else {
            if (inasaIG ||  ishomateIG || isimpleIG) {
                throw CanteraError("newVPSSMgr", "Ideal gas with liquid water");
            } else if (inasaCV || ishomateCV ||  isimpleCV) {
                return new VPSSMgr_General(vp_ptr, spth);
            } else {
                return new VPSSMgr_Water_HKFT(vp_ptr, spth);
            }
        }
    }
    if (inasaCV || ishomateCV || isimpleCV) {
        if (!inasaIG && !ishomateIG && !isimpleIG && !itpx && !ihptx && !iother) {
            return new VPSSMgr_ConstVol(vp_ptr, spth);
        }
    }

    return new VPSSMgr_General(vp_ptr, spth);
}

// I don't think this is currently used. However, this is a virtual
// function where additional capabilities may be added.
VPSSMgr*
VPSSMgrFactory::newVPSSMgr(VPSSMgr_enumType type, VPStandardStateTP* vp_ptr)
{
    SpeciesThermo& spthermoRef = vp_ptr->speciesThermo();
    switch (type) {
    case cVPSSMGR_IDEALGAS:
        return new VPSSMgr_IdealGas(vp_ptr, &spthermoRef);
    case cVPSSMGR_CONSTVOL:
        return new VPSSMgr_ConstVol(vp_ptr, &spthermoRef);
    case cVPSSMGR_PUREFLUID:
        throw CanteraError("VPSSMgrFactory::newVPSSMgr",
                           "unimplemented");
    case cVPSSMGR_WATER_CONSTVOL:
        return new VPSSMgr_Water_ConstVol(vp_ptr, &spthermoRef);
    case cVPSSMGR_WATER_HKFT:
        return new VPSSMgr_Water_HKFT(vp_ptr, &spthermoRef);
    case cVPSSMGR_GENERAL:
        return new VPSSMgr_General(vp_ptr, &spthermoRef);
    case cVPSSMGR_UNDEF:
    default:
        throw UnknownVPSSMgrModel("VPSSMgrFactory::newVPSSMgr", int2str(type));
        return 0;
    }
}

// I don't think this is currently used
VPSSMgr* newVPSSMgr(VPSSMgr_enumType type, VPStandardStateTP* vp_ptr,
                    Cantera::VPSSMgrFactory* f)
{
    if (f == 0) {
        f = VPSSMgrFactory::factory();
    }
    return f->newVPSSMgr(type, vp_ptr);
}


VPSSMgr* newVPSSMgr(VPStandardStateTP* tp_ptr,
                    XML_Node* phaseNode_ptr,
                    std::vector<XML_Node*> & spDataNodeList,
                    VPSSMgrFactory* f)
{
    if (f == 0) {
        f = VPSSMgrFactory::factory();
    }
    return f->newVPSSMgr(tp_ptr, phaseNode_ptr, spDataNodeList);
}

}
