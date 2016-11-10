/**
 *  @file VPSSMgrFactory.cpp
 *    Definitions for factory to build instances of classes that manage the
 *    calculation of standard state properties for all the species in a phase
 *    (see \ref spthermo and class
 *    \link Cantera::VPSSMgrFactory VPSSMgrFactory\endlink);
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "VPSSMgrFactory.h"
#include "cantera/thermo/VPStandardStateTP.h"
#include "cantera/thermo/VPSSMgr_IdealGas.h"
#include "cantera/thermo/VPSSMgr_ConstVol.h"
#include "cantera/thermo/VPSSMgr_Water_ConstVol.h"
#include "cantera/thermo/VPSSMgr_Water_HKFT.h"
#include "cantera/thermo/VPSSMgr_General.h"

#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/thermo/MultiSpeciesThermo.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"

using namespace std;

namespace Cantera
{

VPSSMgrFactory* VPSSMgrFactory::s_factory = 0;

// Defn of the static mutex variable that locks the VPSSMgr factory singleton
std::mutex VPSSMgrFactory::vpss_species_thermo_mutex;

//! Examine the types of species thermo parameterizations, and return a flag
//! indicating the type of parameterization needed by the species.
/*!
 *  @param spDataNodeList            Species Data XML node. This node contains a list
 *                                     of species XML nodes underneath it.
 *  @param has_nasa_idealGas         Boolean indicating that one species has a
 *      NASA ideal gas standard state
 *  @param has_nasa_constVol         Boolean indicating that one species has a
 *      NASA ideal solution standard state
 *  @param has_shomate_idealGas      Boolean indicating that one species has a
 *      Shomate ideal gas standard state
 *  @param has_shomate_constVol      Boolean indicating that one species has a
 *      Shomate ideal solution standard state
 *  @param has_simple_idealGas       Boolean indicating that one species has a
 *      simple ideal gas standard state
 *  @param has_simple_constVol       Boolean indicating that one species has a
 *      simple ideal solution standard state
 *  @param has_water                 Boolean indicating that one species has a
 *      water standard state
 *  @param has_tpx                   Boolean indicating that one species has a
 *      tpx standard state
 *  @param has_hptx                  Boolean indicating that one species has a
 *      htpx standard state
 *  @param has_other                 Boolean indicating that one species has
 *      different standard state than the ones listed above
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
    string ssModel = "idealGas";
    for (size_t n = 0; n < spDataNodeList.size(); n++) {
        bool ifound = false;
        XML_Node* spNode = spDataNodeList[n];
        if (spNode->hasChild("standardState")) {
            string mm = spNode->child("standardState")["model"];
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
                    ssModel = spNode->child("standardState")["model"];
                }
                if (th.hasChild("NASA")) {
                    if (ssModel == "idealGas") {
                        has_nasa_idealGas++;
                    } else if (ssModel == "constant_incompressible" ||
                               ssModel == "constantVolume") {
                        has_nasa_constVol++;
                    } else if (ssModel == "temperature_polynomial" ||
                               ssModel == "density_temperature_polynomial" ||
                               ssModel == "constant") {
                        has_other++;
                    } else {
                        throw CanteraError("getVPSSMgrTypes",
                            "Specified VPSSMgr model {} does not match any known type.",
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
                               ssModel == "density_temperature_polynomial" ||
                               ssModel == "constant") {
                        has_other++;
                    } else {
                        throw CanteraError("getVPSSMgrTypes",
                            "Specified VPSSMgr model {} does not match any known type.",
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
                               ssModel == "density_temperature_polynomial" ||
                               ssModel == "constant") {
                        has_other++;
                    } else {
                        throw CanteraError("getVPSSMgrTypes",
                            "Specified VPSSMgr model {} does not match any known type.",
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
                throw CanteraError("getVPSSMgrTypes",
                    "Specified VPSSMgr model {} does not match any known type.",
                    spNode->attrib("name"));
            }
        }
    }
}

VPSSMgrFactory::VPSSMgrFactory()
{
    reg("idealgas",
        [] (VPStandardStateTP* tp, MultiSpeciesThermo* st) {
            return new VPSSMgr_IdealGas(tp, st); });
    reg("constvol",
        [] (VPStandardStateTP* tp, MultiSpeciesThermo* st) {
            return new VPSSMgr_ConstVol(tp, st); });
    reg("water_constvol",
        [] (VPStandardStateTP* tp, MultiSpeciesThermo* st) {
            return new VPSSMgr_Water_ConstVol(tp, st); });
    reg("water_hkft",
        [] (VPStandardStateTP* tp, MultiSpeciesThermo* st) {
            return new VPSSMgr_Water_HKFT(tp, st); });
    reg("general",
        [] (VPStandardStateTP* tp, MultiSpeciesThermo* st) {
            return new VPSSMgr_General(tp, st); });
}

void VPSSMgrFactory::deleteFactory()
{
    std::unique_lock<std::mutex> lock(vpss_species_thermo_mutex);
    delete s_factory;
    s_factory = 0;
}

VPSSMgr_enumType
VPSSMgrFactory::VPSSMgr_StringConversion(const std::string& ssModel) const
{
    warn_deprecated("VPSSMgrFactory::VPSSMgr_StringConversion",
        "To be removed after Cantera 2.3.");
    std::string lssModel = ba::to_lower_copy(ssModel);
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

VPSSMgr* VPSSMgrFactory::newVPSSMgr(VPStandardStateTP* vp_ptr,
                                    XML_Node* phaseNode_ptr,
                                    std::vector<XML_Node*> & spDataNodeList)
{
    std::string ssManager;
    std::string vpssManager;

    // First look for any explicit instructions within the XML Database for the
    // standard state manager and the variable pressure standard state manager
    if (phaseNode_ptr && phaseNode_ptr->hasChild("thermo")) {
        const XML_Node& thermoNode = phaseNode_ptr->child("thermo");
        if (thermoNode.hasChild("standardStateManager")) {
            const XML_Node& ssNode = thermoNode.child("standardStateManager");
            ssManager = ssNode["model"];
        }
        if (thermoNode.hasChild("variablePressureStandardStateManager")) {
            const XML_Node& vpssNode = thermoNode.child("variablePressureStandardStateManager");
            vpssManager = ba::to_lower_copy(vpssNode["model"]);
        }
    }

    // first get the reference state handler.
    MultiSpeciesThermo* spth = &vp_ptr->speciesThermo();

    // Next, if we have specific directions, use them to get the VPSSSMgr object
    // and return immediately
    if (vpssManager != "") {
        return create(vpssManager, vp_ptr, spth);
    }

    int inasaIG = 0, inasaCV = 0, ishomateIG = 0, ishomateCV = 0,
        isimpleIG = 0, isimpleCV = 0, iwater = 0, itpx = 0, iother = 0;
    int ihptx = 0;

    try {
        getVPSSMgrTypes(spDataNodeList, inasaIG, inasaCV, ishomateIG, ishomateCV,
                        isimpleIG, isimpleCV, iwater, itpx, ihptx, iother);
    } catch (CanteraError) {
        iother = 1;
    }

    if (iwater == 1) {
        if (ihptx == 0) {
            if (inasaIG || ishomateIG || isimpleIG) {
                throw CanteraError("newVPSSMgr", "Ideal gas with liquid water");
            } else {
                return new VPSSMgr_Water_ConstVol(vp_ptr, spth);
            }
        } else {
            if (inasaIG || ishomateIG || isimpleIG) {
                throw CanteraError("newVPSSMgr", "Ideal gas with liquid water");
            } else if (inasaCV || ishomateCV || isimpleCV) {
                return new VPSSMgr_General(vp_ptr, spth);
            } else {
                return new VPSSMgr_Water_HKFT(vp_ptr, spth);
            }
        }
    }
    if ((inasaCV || ishomateCV || isimpleCV) &&
        !inasaIG && !ishomateIG && !isimpleIG && !itpx && !ihptx && !iother) {
        return new VPSSMgr_ConstVol(vp_ptr, spth);
    }
    return new VPSSMgr_General(vp_ptr, spth);
}

// I don't think this is currently used. However, this is a virtual
// function where additional capabilities may be added.
VPSSMgr* VPSSMgrFactory::newVPSSMgr(VPSSMgr_enumType type,
                                    VPStandardStateTP* vp_ptr)
{
    warn_deprecated("VPSSMgrFactory::newVPSSMgr(VPSSMgr_enumType, VPStandardStateTP*)",
        "To be removed after Cantera 2.3.");
    static unordered_map<int, std::string> types {
        {cVPSSMGR_IDEALGAS, "idealgas"},
        {cVPSSMGR_CONSTVOL, "constvol"},
        {cVPSSMGR_WATER_CONSTVOL, "water_constvol"},
        {cVPSSMGR_WATER_HKFT, "water_hkft"},
        {cVPSSMGR_GENERAL, "general"}
    };
    MultiSpeciesThermo& spthermoRef = vp_ptr->speciesThermo();
    return create(types.at(type), vp_ptr, &spthermoRef);
}

// I don't think this is currently used
VPSSMgr* newVPSSMgr(VPSSMgr_enumType type, VPStandardStateTP* vp_ptr,
                    Cantera::VPSSMgrFactory* f)
{
    warn_deprecated("newVPSSMgr(VPSSMgr_enumType, VPStandardStateTP*, VPSSMgrFactory*)",
        "To be removed after Cantera 2.3.");
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
    } else {
        warn_deprecated("newVPSSMgr(VPStandardStateTP*, XML_Node*, vector<XML_Node*>, VPSSMgrFactory*)",
            "The `VPSSMgrFactory*` argument to this function is deprecated and"
            " will be removed after Cantera 2.3.");
    }
    return f->newVPSSMgr(tp_ptr, phaseNode_ptr, spDataNodeList);
}

}
