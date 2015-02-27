/**
 *  @file SpeciesThermoFactory.cpp
 *    Definitions for factory to build instances of classes that manage the
 *    standard-state thermodynamic properties of a set of species
 *    (see \ref spthermo and class \link Cantera::SpeciesThermoFactory SpeciesThermoFactory\endlink);
 */
// Copyright 2001  California Institute of Technology

#include "cantera/thermo/SpeciesThermoFactory.h"

#include "cantera/thermo/SpeciesThermo.h"
#include "NasaThermo.h"
#include "ShomateThermo.h"
#include "cantera/thermo/SimpleThermo.h"
#include "cantera/thermo/GeneralSpeciesThermo.h"
#include "cantera/thermo/Mu0Poly.h"
#include "cantera/thermo/Nasa9PolyMultiTempRegion.h"
#include "cantera/thermo/Nasa9Poly1.h"
#include "cantera/thermo/StatMech.h"
#include "cantera/thermo/NasaPoly2.h"
#include "cantera/thermo/ConstCpPoly.h"
#include "cantera/thermo/AdsorbateThermo.h"
#include "cantera/thermo/SpeciesThermoMgr.h"
#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/thermo/VPSSMgr.h"
#include "cantera/thermo/VPStandardStateTP.h"

#include "cantera/base/ctml.h"

using namespace std;
using namespace ctml;

namespace Cantera
{

SpeciesThermoFactory* SpeciesThermoFactory::s_factory = 0;
mutex_t SpeciesThermoFactory::species_thermo_mutex;

//! Examine the types of species thermo parameterizations,
//! and return a flag indicating the type of reference state thermo manager
//! that will be needed in order to evaluate them all.
/*!
 *
 *  @param spDataNodeList    This vector contains a list
 *                           of species XML nodes that will be in the phase
 *  @param has_nasa          Return int that indicates whether the phase has a NASA polynomial form for one of its species
 *  @param has_shomate       Return int that indicates whether the phase has a SHOMATE polynomial form for one of its species
 *  @param has_simple        Return int that indicates whether the phase has a SIMPLE polynomial form for one of its species
 *  @param has_other         Return int that indicates whether the phase has a form for one of its species that is not one of the ones listed above.
 *
 * @todo Make sure that spDadta_node is species Data XML node by checking its name is speciesData
 * @deprecated
 */
static void getSpeciesThermoTypes(std::vector<XML_Node*> & spDataNodeList,
                                  int& has_nasa, int& has_shomate, int& has_simple,
                                  int& has_other)
{
    for (size_t n = 0; n < spDataNodeList.size(); n++) {
        XML_Node* spNode = spDataNodeList[n];
        if (spNode->hasChild("standardState")) {
            string mname = spNode->child("standardState")["model"];
            if (mname == "water" || mname == "waterIAPWS") {
                has_other = 1;
                continue;
            }
        }
        if (spNode->hasChild("thermo")) {
            const XML_Node& th = spNode->child("thermo");
            if (th.hasChild("NASA")) {
                has_nasa = 1;
            } else if (th.hasChild("Shomate")) {
                has_shomate = 1;
            } else if (th.hasChild("MinEQ3")) {
                has_shomate = 1;
            } else if (th.hasChild("const_cp")) {
                has_simple = 1;
            } else if (th.hasChild("poly")) {
                if (th.child("poly")["order"] == "1") {
                    has_simple = 1;
                } else throw CanteraError("newSpeciesThermo",
                                              "poly with order > 1 not yet supported");
            } else if (th.hasChild("Mu0")) {
                has_other = 1;
            } else if (th.hasChild("NASA9")) {
                has_other = 1;
            } else if (th.hasChild("NASA9MULTITEMP")) {
                has_other = 1;
            } else if (th.hasChild("adsorbate")) {
                has_other = 1;
            } else {
                has_other = 1;
            }
        } else {
            throw CanteraError("getSpeciesThermoTypes:",
                               spNode->attrib("name") + " is missing the thermo XML node");
        }
    }
}

SpeciesThermoFactory* SpeciesThermoFactory::factory()
{
    warn_deprecated("class SpeciesThermoFactory",
                    "To be removed after Cantera 2.2.");
    ScopedLock lock(species_thermo_mutex);
    if (!s_factory) {
        s_factory = new SpeciesThermoFactory;
    }
    return s_factory;
}

void SpeciesThermoFactory::deleteFactory()
{
    ScopedLock lock(species_thermo_mutex);
    delete s_factory;
    s_factory = 0;
}

SpeciesThermo* SpeciesThermoFactory::newSpeciesThermo(std::vector<XML_Node*> & spDataNodeList) const
{
    warn_deprecated("SpeciesThermoFactory::newSpeciesThermo",
        "To be removed after Cantera 2.2. Use class GeneralSpeciesThermo directly.");
    int inasa = 0, ishomate = 0, isimple = 0, iother = 0;
    try {
        getSpeciesThermoTypes(spDataNodeList, inasa, ishomate, isimple, iother);
    } catch (UnknownSpeciesThermoModel) {
        iother = 1;
        popError();
    }
    if (iother) {
        return new GeneralSpeciesThermo();
    }
    return newSpeciesThermo(NASA*inasa
                            + SHOMATE*ishomate + SIMPLE*isimple);
}

SpeciesThermo* SpeciesThermoFactory::newSpeciesThermo(int type) const
{
    warn_deprecated("SpeciesThermoFactory::newSpeciesThermo",
        "To be removed after Cantera 2.2. Use class GeneralSpeciesThermo directly.");
    switch (type) {
    case NASA:
        return new NasaThermo;
    case SHOMATE:
        return new ShomateThermo;
    case SIMPLE:
        return new SimpleThermo;
    case NASA + SHOMATE:
        return new SpeciesThermoDuo<NasaThermo, ShomateThermo>;
    case NASA + SIMPLE:
        return new SpeciesThermoDuo<NasaThermo, SimpleThermo>;
    case SHOMATE + SIMPLE:
        return new SpeciesThermoDuo<ShomateThermo, SimpleThermo>;
    default:
        throw UnknownSpeciesThermo("SpeciesThermoFactory::newSpeciesThermo",
                                   type);
        return 0;
    }
}

SpeciesThermo* SpeciesThermoFactory::newSpeciesThermoManager(const std::string& stype) const
{
    warn_deprecated("SpeciesThermoFactory::newSpeciesThermo",
        "To be removed after Cantera 2.2. Use class GeneralSpeciesThermo directly.");
    std::string ltype = lowercase(stype);
    if (ltype == "nasa") {
        return new NasaThermo;
    } else if (ltype == "shomate") {
        return new ShomateThermo;
    } else if (ltype ==  "simple" || ltype == "constant_cp") {
        return new SimpleThermo;
    } else if (ltype ==  "nasa_shomate_duo") {
        return new SpeciesThermoDuo<NasaThermo, ShomateThermo>;
    } else if (ltype ==  "nasa_simple_duo") {
        return new SpeciesThermoDuo<NasaThermo, SimpleThermo>;
    } else if (ltype ==  "shomate_simple_duo") {
        return new SpeciesThermoDuo<ShomateThermo, SimpleThermo>;
    } else if (ltype ==   "general") {
        return new GeneralSpeciesThermo();
    } else if (ltype ==  "") {
        return (SpeciesThermo*) 0;
    } else {
        throw UnknownSpeciesThermo("SpeciesThermoFactory::newSpeciesThermoManager",
                                   stype);
    }
    return (SpeciesThermo*) 0;
}

SpeciesThermoInterpType* newSpeciesThermoInterpType(int type, double tlow,
    double thigh, double pref, const double* coeffs)
{
    switch (type) {
    case NASA1:
        return new NasaPoly1(tlow, thigh, pref, coeffs);
    case SHOMATE1:
        return new ShomatePoly(tlow, thigh, pref, coeffs);
    case CONSTANT_CP:
    case SIMPLE:
        return new ConstCpPoly(tlow, thigh, pref, coeffs);
    case MU0_INTERP:
        return new Mu0Poly(tlow, thigh, pref, coeffs);
    case SHOMATE2:
        return new ShomatePoly2(tlow, thigh, pref, coeffs);
    case NASA2:
        return new NasaPoly2(tlow, thigh, pref, coeffs);
    case ADSORBATE:
        return new Adsorbate(tlow, thigh, pref, coeffs);
    default:
        throw CanteraError("newSpeciesThermoInterpType",
                           "Unknown species thermo type: " + int2str(type) + ".");
    }
}

SpeciesThermoInterpType* newSpeciesThermoInterpType(const std::string& stype,
    double tlow, double thigh, double pref, const double* coeffs)
{
    int itype = -1;
    std::string type = lowercase(stype);
    if (type == "nasa2" || type == "nasa") {
        itype = NASA2; // two-region 7-coefficient NASA polynomials
    } else if (type == "const_cp" || type == "simple") {
        itype = CONSTANT_CP;
    } else if (type == "shomate" || type == "shomate1") {
        itype = SHOMATE1; // single-region Shomate polynomial
    } else if (type == "shomate2") {
        itype = SHOMATE2; // two-region Shomate polynomials
    } else if (type == "nasa1") {
        itype = NASA1; // single-region, 7-coefficient NASA polynomial
    } else if (type == "nasa9") {
        itype = NASA9; // single-region, 9-coefficient NASA polynomial
    } else if (type == "nasa9multi") {
        itype = NASA9MULTITEMP; // multi-region, 9-coefficient NASA polynomials
    } else if (type == "mu0") {
        itype = MU0_INTERP;
    } else if (type == "adsorbate") {
        itype = ADSORBATE;
    } else {
        throw CanteraError("newSpeciesThermoInterpType",
                           "Unknown species thermo type: '" + stype + "'.");
    }
    return newSpeciesThermoInterpType(itype, tlow, thigh, pref, coeffs);
}

//! Create a NASA polynomial thermodynamic property parameterization for a
//! species from a set ! of XML nodes
/*!
 * This is called if a 'NASA' node is found in the XML input.
 *
 *  @param speciesName  name of the species
 *  @param nodes        vector of 1 or 2 'NASA' XML_Nodes, each defining the
 *      coefficients for a temperature range
 */
static SpeciesThermoInterpType* newNasaThermoFromXML(
    const std::string& speciesName, vector<XML_Node*> nodes)
{
    const XML_Node& f0 = *nodes[0];
    bool dualRange = (nodes.size() > 1);

    double tmin0 = fpValue(f0["Tmin"]);
    double tmax0 = fpValue(f0["Tmax"]);

    doublereal p0 = OneAtm;
    if (f0.hasAttrib("P0")) {
        p0 = fpValue(f0["P0"]);
    }
    if (f0.hasAttrib("Pref")) {
        p0 = fpValue(f0["Pref"]);
    }
    p0 = OneAtm;

    double tmin1 = tmax0;
    double tmax1 = tmin1 + 0.0001;
    if (dualRange) {
        tmin1 = fpValue(nodes[1]->attrib("Tmin"));
        tmax1 = fpValue(nodes[1]->attrib("Tmax"));
    }

    vector_fp c0, c1;
    doublereal tmin, tmid, tmax;
    if (fabs(tmax0 - tmin1) < 0.01) {
        // f0 has the lower T data, and f1 the higher T data
        tmin = tmin0;
        tmid = tmax0;
        tmax = tmax1;
        getFloatArray(f0.child("floatArray"), c0, false);
        if (dualRange) {
            getFloatArray(nodes[1]->child("floatArray"), c1, false);
        } else {
            // if there is no higher range data, then copy c0 to c1.
            c1.resize(7,0.0);
            copy(c0.begin(), c0.end(), c1.begin());
        }
    } else if (fabs(tmax1 - tmin0) < 0.01) {
        // f1 has the lower T data, and f0 the higher T data
        tmin = tmin1;
        tmid = tmax1;
        tmax = tmax0;
        getFloatArray(nodes[1]->child("floatArray"), c0, false);
        getFloatArray(f0.child("floatArray"), c1, false);
    } else {
        throw CanteraError("installNasaThermo",
                           "non-continuous temperature ranges.");
    }

    vector_fp c(15);
    c[0] = tmid;
    copy(c1.begin(), c1.begin()+7, c.begin() + 1); // high-T coefficients
    copy(c0.begin(), c0.begin()+7, c.begin() + 8); // low-T coefficients
    return newSpeciesThermoInterpType(NASA, tmin, tmax, p0, &c[0]);
}

//! Create a Shomate polynomial from an XML node giving the 'EQ3' coefficients
/*!
 *  This is called if a 'MinEQ3' node is found in the XML input.
 *
 *  @param name         name of the species
 *  @param MinEQ3node   The XML_Node containing the MinEQ3 parameterization
 */
SpeciesThermoInterpType* newShomateForMineralEQ3(const std::string& name,
                                                 const XML_Node& MinEQ3node)
{
    doublereal tmin0 = strSItoDbl(MinEQ3node["Tmin"]);
    doublereal tmax0 = strSItoDbl(MinEQ3node["Tmax"]);
    doublereal p0 = strSItoDbl(MinEQ3node["Pref"]);

    doublereal deltaG_formation_pr_tr =
        getFloatDefaultUnits(MinEQ3node, "DG0_f_Pr_Tr", "cal/gmol", "actEnergy");
    doublereal deltaH_formation_pr_tr =
        getFloatDefaultUnits(MinEQ3node, "DH0_f_Pr_Tr", "cal/gmol", "actEnergy");
    doublereal Entrop_pr_tr = getFloatDefaultUnits(MinEQ3node, "S0_Pr_Tr", "cal/gmol/K");
    doublereal a = getFloatDefaultUnits(MinEQ3node, "a", "cal/gmol/K");
    doublereal b = getFloatDefaultUnits(MinEQ3node, "b", "cal/gmol/K2");
    doublereal c = getFloatDefaultUnits(MinEQ3node, "c", "cal-K/gmol");
    doublereal dg = deltaG_formation_pr_tr * 4.184 * 1.0E3;
    doublereal DHjmol = deltaH_formation_pr_tr * 1.0E3 * 4.184;
    doublereal fac = DHjmol - dg - 298.15 * Entrop_pr_tr * 1.0E3 * 4.184;
    doublereal Mu0_tr_pr = fac + dg;
    doublereal e = Entrop_pr_tr * 1.0E3 * 4.184;
    doublereal Hcalc = Mu0_tr_pr + 298.15 * e;

    /*
     * Now calculate the shomate polynomials
     *
     * Cp first
     *
     *  Shomate: (Joules / gmol / K)
     *    Cp = As + Bs * t + Cs * t*t + Ds * t*t*t + Es / (t*t)
     *     where
     *          t = temperature(Kelvin) / 1000
     */
    double As = a * 4.184;
    double Bs = b * 4.184 * 1000.;
    double Cs = 0.0;
    double Ds = 0.0;
    double Es = c * 4.184 / (1.0E6);

    double t = 298.15 / 1000.;
    double H298smFs = As * t + Bs * t * t / 2.0 - Es / t;

    double HcalcS = Hcalc / 1.0E6;
    double Fs = HcalcS - H298smFs;

    double S298smGs = As * log(t) + Bs * t - Es/(2.0*t*t);
    double ScalcS = e / 1.0E3;
    double Gs = ScalcS - S298smGs;

    double c0[7] = {As, Bs, Cs, Ds, Es, Fs, Gs};

    return newSpeciesThermoInterpType(SHOMATE1, tmin0, tmax0, p0, c0);
}

//! Create a Shomate polynomial thermodynamic property parameterization for a
//! species
/*!
 *  This is called if a 'Shomate' node is found in the XML input.
 *
 *  @param speciesName  name of the species
 *  @param nodes        vector of 1 or 2 'Shomate' XML_Nodes, each defining the
 *      coefficients for a temperature range
 */
static SpeciesThermoInterpType* newShomateThermoFromXML(
    const std::string& speciesName, vector<XML_Node*>& nodes)
{
    bool dualRange = false;
    if (nodes.size() == 2) {
        dualRange = true;
    }
    double tmin0 = fpValue(nodes[0]->attrib("Tmin"));
    double tmax0 = fpValue(nodes[0]->attrib("Tmax"));

    doublereal p0 = OneAtm;
    if (nodes[0]->hasAttrib("P0")) {
        p0 = fpValue(nodes[0]->attrib("P0"));
    }
    if (nodes[0]->hasAttrib("Pref")) {
        p0 = fpValue(nodes[0]->attrib("Pref"));
    }
    p0 = OneAtm;

    double tmin1 = tmax0;
    double tmax1 = tmin1 + 0.0001;
    if (dualRange) {
        tmin1 = fpValue(nodes[1]->attrib("Tmin"));
        tmax1 = fpValue(nodes[1]->attrib("Tmax"));
    }

    vector_fp c0, c1;
    doublereal tmin, tmid, tmax;
    if (fabs(tmax0 - tmin1) < 0.01) {
        tmin = tmin0;
        tmid = tmax0;
        tmax = tmax1;
        getFloatArray(nodes[0]->child("floatArray"), c0, false);
        if (dualRange) {
            getFloatArray(nodes[1]->child("floatArray"), c1, false);
        } else {
            if(c0.size() != 7)
            {
              throw CanteraError("installShomateThermoFromXML",
                                 "Shomate thermo requires 7 coefficients in float array.");
            }
            c1.resize(7,0.0);
            copy(c0.begin(), c0.begin()+7, c1.begin());
        }
    } else if (fabs(tmax1 - tmin0) < 0.01) {
        tmin = tmin1;
        tmid = tmax1;
        tmax = tmax0;
        getFloatArray(nodes[1]->child("floatArray"), c0, false);
        getFloatArray(nodes[0]->child("floatArray"), c1, false);
    } else {
        throw CanteraError("installShomateThermoFromXML",
                           "non-continuous temperature ranges.");
    }
    if(c0.size() != 7 || c1.size() != 7)
    {
      throw CanteraError("installShomateThermoFromXML",
                         "Shomate thermo requires 7 coefficients in float array.");
    }
    vector_fp c(15);
    c[0] = tmid;
    copy(c0.begin(), c0.begin()+7, c.begin() + 1);
    copy(c1.begin(), c1.begin()+7, c.begin() + 8);
    return newSpeciesThermoInterpType(SHOMATE, tmin, tmax, p0, &c[0]);
}

//! Create a "simple" constant heat capacity thermodynamic property
//! parameterization for a ! species
/*!
 * This is called if a 'const_cp' XML node is found
 *
 *  @param speciesName  name of the species
 *  @param f            'const_cp' XML node
 */
static SpeciesThermoInterpType* newConstCpThermoFromXML(
    const std::string& speciesName, XML_Node& f)
{
    double tmin = fpValue(f["Tmin"]);
    double tmax = fpValue(f["Tmax"]);
    if (tmax == 0.0) {
        tmax = 1.0e30;
    }

    vector_fp c(4);
    c[0] = getFloat(f, "t0", "toSI");
    c[1] = getFloat(f, "h0", "toSI");
    c[2] = getFloat(f, "s0", "toSI");
    c[3] = getFloat(f, "cp0", "toSI");
    doublereal p0 = OneAtm;
    return newSpeciesThermoInterpType(CONSTANT_CP, tmin, tmax, p0, &c[0]);
}


//! Create a NASA9 polynomial thermodynamic property parameterization for a
//! species
/*!
 *  This is called if a 'NASA9' Node is found in the XML input.
 *
 *  @param speciesName  name of the species
 *  @param tp           Vector of XML Nodes that make up the parameterization
 */
static SpeciesThermoInterpType* newNasa9ThermoFromXML(
    const std::string& speciesName, const std::vector<XML_Node*>& tp)
{
    int nRegions = 0;
    vector_fp cPoly;
    std::vector<Nasa9Poly1*> regionPtrs;
    doublereal pref = OneAtm;
    // Loop over all of the possible temperature regions
    for (size_t i = 0; i < tp.size(); i++) {
        const XML_Node& fptr = *tp[i];
        if (fptr.name() == "NASA9" && fptr.hasChild("floatArray")) {
            double tmin = fpValue(fptr["Tmin"]);
            double tmax = fpValue(fptr["Tmax"]);
            if (fptr.hasAttrib("P0")) {
                pref = fpValue(fptr["P0"]);
            }
            if (fptr.hasAttrib("Pref")) {
                pref = fpValue(fptr["Pref"]);
            }

            getFloatArray(fptr.child("floatArray"), cPoly, false);
            if (cPoly.size() != 9) {
                throw CanteraError("installNasa9ThermoFromXML",
                                   "Expected 9 coeff polynomial");
            }
            regionPtrs.push_back(new Nasa9Poly1(tmin, tmax, pref, &cPoly[0]));
            nRegions++;
        }
    }
    if (nRegions == 0) {
        throw UnknownSpeciesThermoModel("installThermoForSpecies",
                                        speciesName, "  ");
    } else if (nRegions == 1)  {
        return regionPtrs[0];
    } else {
        return new Nasa9PolyMultiTempRegion(regionPtrs);
    }
}

/**
 * Create a stat mech based property solver for a species
 * @deprecated
 */
static StatMech* newStatMechThermoFromXML(const std::string& speciesName, XML_Node& f)
{
    doublereal tmin = fpValue(f["Tmin"]);
    doublereal tmax = fpValue(f["Tmax"]);
    doublereal pref = OneAtm;
    if (f.hasAttrib("P0")) {
        pref = fpValue(f["P0"]);
    }
    if (f.hasAttrib("Pref")) {
        pref = fpValue(f["Pref"]);
    }

    // set properties
    tmin = 0.1;
    vector_fp coeffs(1);
    coeffs[0] = 0.0;
    return new StatMech(0, tmin, tmax, pref, &coeffs[0], speciesName);
}

//! Create an Adsorbate polynomial thermodynamic property parameterization for a
//! species
/*!
 * This is called if a 'Adsorbate' node is found in the XML input.
 *
 *  @param speciesName  name of the species
 *  @param f            XML Node that contains the parameterization
 */
static SpeciesThermoInterpType* newAdsorbateThermoFromXML(
    const std::string& speciesName, const XML_Node& f)
{
    vector_fp freqs;
    doublereal pref = OneAtm;
    double tmin = fpValue(f["Tmin"]);
    double tmax = fpValue(f["Tmax"]);
    if (f.hasAttrib("P0")) {
        pref = fpValue(f["P0"]);
    }
    if (f.hasAttrib("Pref")) {
        pref = fpValue(f["Pref"]);
    }
    if (tmax == 0.0) {
        tmax = 1.0e30;
    }

    if (f.hasChild("floatArray")) {
        getFloatArray(f.child("floatArray"), freqs, false);
    }
    for (size_t n = 0; n < freqs.size(); n++) {
        freqs[n] *= 3.0e10;
    }
    vector_fp coeffs(freqs.size() + 2);
    coeffs[0] = static_cast<double>(freqs.size());
    coeffs[1] = getFloat(f, "binding_energy", "toSI");
    copy(freqs.begin(), freqs.end(), coeffs.begin() + 2);
    return new Adsorbate(0, tmin, tmax, pref, &coeffs[0]);
}

void SpeciesThermoFactory::installThermoForSpecies
(size_t k, const XML_Node& speciesNode, ThermoPhase* th_ptr,
 SpeciesThermo& spthermo, const XML_Node* phaseNode_ptr) const
{
    SpeciesThermoInterpType* stit = newSpeciesThermoInterpType(speciesNode);
    stit->validate(speciesNode["name"]);
    spthermo.install_STIT(k, stit);
}

void SpeciesThermoFactory::installVPThermoForSpecies(size_t k,
    const XML_Node& speciesNode,
    VPStandardStateTP* vp_ptr,
    VPSSMgr* vpssmgr_ptr,
    SpeciesThermo* spthermo_ptr,
    const XML_Node* phaseNode_ptr) const
{
    warn_deprecated("SpeciesThermoFactory::installVPThermoForSpecies",
                    "Call VPStandardStateTP::createInstallPDSS directly.");
    // Call the VPStandardStateTP object to install the pressure dependent species
    // standard state into the object.
    //
    // We don't need to pass spthermo_ptr down, because it's already installed
    // into vp_ptr.
    //
    // We don't need to pass vpssmgr_ptr down, because it's already installed
    // into vp_ptr.
    vp_ptr->createInstallPDSS(k, speciesNode,  phaseNode_ptr);
}

SpeciesThermoInterpType* newSpeciesThermoInterpType(const XML_Node& speciesNode)
{
    /*
     * Check to see that the species block has a thermo block
     * before processing. Throw an error if not there.
     */
    if (!(speciesNode.hasChild("thermo"))) {
        throw UnknownSpeciesThermoModel("installThermoForSpecies",
                                        speciesNode["name"], "<nonexistent>");
    }
    const XML_Node& thermo = speciesNode.child("thermo");

    // Get the children of the thermo XML node. In the next bit of code we take out the comments that
    // may have been children of the thermo XML node by doing a selective copy.
    // These shouldn't interfere with the algorithm at any point.
    const std::vector<XML_Node*>& tpWC = thermo.children();
    std::vector<XML_Node*> tp;
    for (size_t i = 0; i < tpWC.size(); i++) {
        if (!(tpWC[i])->isComment()) {
            tp.push_back(tpWC[i]);
        }
    }

    std::string thermoType = lowercase(tp[0]->name());
    std::string specName = speciesNode["name"];

    for (size_t i = 1; i < tp.size(); i++) {
        if (lowercase(tp[i]->name()) != thermoType) {
            throw CanteraError("newSpeciesThermoInterpType",
                "Encounter unsupported mixed species thermo parameterizations "
                "for species '" + specName + "'.");
        }
    }
    if ((tp.size() > 2 && thermoType != "nasa9") ||
        (tp.size() > 1 && (thermoType == "const_cp" ||
                           thermoType == "mu0" ||
                           thermoType == "adsorbate"))) {
        throw CanteraError("newSpeciesThermoInterpType",
            "Too many regions in thermo parameterization for species '" +
            specName + "'.");
    }

    if (thermo["model"] == "MineralEQ3") {
        if (thermoType != "mineq3") {
            throw CanteraError("SpeciesThermoFactory::installThermoForSpecies",
                               "confused: expected MinEQ3");
        }
        return newShomateForMineralEQ3(specName, *tp[0]);
    } else if (thermoType == "shomate") {
        return newShomateThermoFromXML(specName, tp);
    } else if (thermoType == "const_cp") {
        return newConstCpThermoFromXML(specName, *tp[0]);
    } else if (thermoType == "nasa") {
        return newNasaThermoFromXML(specName, tp);
    } else if (thermoType == "mu0") {
        return newMu0ThermoFromXML(specName, *tp[0]);
    } else if (thermoType == "nasa9") {
        return newNasa9ThermoFromXML(specName, tp);
    } else if (thermoType == "adsorbate") {
        return newAdsorbateThermoFromXML(specName, *tp[0]);
    } else if (thermoType == "statmech") {
        return newStatMechThermoFromXML(specName, *tp[0]);
    } else {
        throw UnknownSpeciesThermoModel("installThermoForSpecies",
                                        specName, thermoType);
    }
}

SpeciesThermo* newSpeciesThermoMgr(int type, SpeciesThermoFactory* f)
{
    warn_deprecated("newSpeciesThermoMgr", "To be removed after Cantera 2.2. "
        "Use class GeneralSpeciesThermo directly.");
    if (f == 0) {
        f = SpeciesThermoFactory::factory();
    }
    return f->newSpeciesThermo(type);
}

SpeciesThermo* newSpeciesThermoMgr(const std::string& stype,
                                   SpeciesThermoFactory* f)
{
    warn_deprecated("newSpeciesThermoMgr", "To be removed after Cantera 2.2. "
        "Use class GeneralSpeciesThermo directly.");
    if (f == 0) {
        f = SpeciesThermoFactory::factory();
    }
    return f->newSpeciesThermoManager(stype);
}

SpeciesThermo* newSpeciesThermoMgr(std::vector<XML_Node*> spData_nodes,
                                   SpeciesThermoFactory* f)
{
    warn_deprecated("newSpeciesThermoMgr", "To be removed after Cantera 2.2. "
        "Use class GeneralSpeciesThermo directly.");
    if (f == 0) {
        f = SpeciesThermoFactory::factory();
    }
    return f->newSpeciesThermo(spData_nodes);
}

}
