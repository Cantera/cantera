/**
 *  @file SpeciesThermoFactory.cpp
 *    Definitions for factory functions to build instances of classes that
 *    manage the standard-state thermodynamic properties of a set of species
 *    (see \ref spthermo);
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/thermo/MultiSpeciesThermo.h"
#include "cantera/thermo/Mu0Poly.h"
#include "cantera/thermo/Nasa9PolyMultiTempRegion.h"
#include "cantera/thermo/Nasa9Poly1.h"
#include "cantera/thermo/NasaPoly2.h"
#include "cantera/thermo/ShomatePoly.h"
#include "cantera/thermo/ConstCpPoly.h"
#include "cantera/thermo/AdsorbateThermo.h"
#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/thermo/VPStandardStateTP.h"
#include "cantera/base/ctml.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{

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
                           "Unknown species thermo type: {}.", type);
    }
}

SpeciesThermoInterpType* newSpeciesThermoInterpType(const std::string& stype,
    double tlow, double thigh, double pref, const double* coeffs)
{
    int itype = -1;
    std::string type = toLowerCopy(stype);
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
 *  @param nodes        vector of 1 or 2 'NASA' XML_Nodes, each defining the
 *      coefficients for a temperature range
 */
static SpeciesThermoInterpType* newNasaThermoFromXML(vector<XML_Node*> nodes)
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
            c1 = c0;
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
 *  @param MinEQ3node   The XML_Node containing the MinEQ3 parameterization
 */
SpeciesThermoInterpType* newShomateForMineralEQ3(const XML_Node& MinEQ3node)
{
    doublereal tmin0 = strSItoDbl(MinEQ3node["Tmin"]);
    doublereal tmax0 = strSItoDbl(MinEQ3node["Tmax"]);
    doublereal p0 = strSItoDbl(MinEQ3node["Pref"]);

    doublereal deltaG_formation_pr_tr =
        getFloat(MinEQ3node, "DG0_f_Pr_Tr", "actEnergy") / actEnergyToSI("cal/gmol");
    doublereal deltaH_formation_pr_tr =
        getFloat(MinEQ3node, "DH0_f_Pr_Tr", "actEnergy") / actEnergyToSI("cal/gmol");
    doublereal Entrop_pr_tr = getFloat(MinEQ3node, "S0_Pr_Tr", "toSI") / toSI("cal/gmol/K");
    doublereal a = getFloat(MinEQ3node, "a", "toSI") / toSI("cal/gmol/K");
    doublereal b = getFloat(MinEQ3node, "b", "toSI") / toSI("cal/gmol/K2");
    doublereal c = getFloat(MinEQ3node, "c", "toSI") / toSI("cal-K/gmol");
    doublereal dg = deltaG_formation_pr_tr * toSI("cal/gmol");
    doublereal DHjmol = deltaH_formation_pr_tr * toSI("cal/gmol");
    doublereal fac = DHjmol - dg - 298.15 * Entrop_pr_tr * toSI("cal/gmol");
    doublereal Mu0_tr_pr = fac + dg;
    doublereal e = Entrop_pr_tr * toSI("cal/gmol");
    doublereal Hcalc = Mu0_tr_pr + 298.15 * e;

    // Now calculate the shomate polynomials
    //
    // Cp first
    //
    //  Shomate: (Joules / gmol / K)
    //    Cp = As + Bs * t + Cs * t*t + Ds * t*t*t + Es / (t*t)
    //     where
    //          t = temperature(Kelvin) / 1000
    double As = a * toSI("cal");
    double Bs = b * toSI("cal") * 1000.;
    double Cs = 0.0;
    double Ds = 0.0;
    double Es = c * toSI("cal") / (1.0E6);

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
 *  @param nodes        vector of 1 or 2 'Shomate' XML_Nodes, each defining the
 *      coefficients for a temperature range
 */
static SpeciesThermoInterpType* newShomateThermoFromXML(
    vector<XML_Node*>& nodes)
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
 *  @param f            'const_cp' XML node
 */
static SpeciesThermoInterpType* newConstCpThermoFromXML(XML_Node& f)
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
 *  @param tp           Vector of XML Nodes that make up the parameterization
 */
static SpeciesThermoInterpType* newNasa9ThermoFromXML(
    const std::vector<XML_Node*>& tp)
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
        throw CanteraError("newNasa9ThermoFromXML", "zero regions found");
    } else if (nRegions == 1) {
        return regionPtrs[0];
    } else {
        return new Nasa9PolyMultiTempRegion(regionPtrs);
    }
}

//! Create an Adsorbate polynomial thermodynamic property parameterization for a
//! species
/*!
 * This is called if a 'Adsorbate' node is found in the XML input.
 *
 *  @param f            XML Node that contains the parameterization
 */
static SpeciesThermoInterpType* newAdsorbateThermoFromXML(const XML_Node& f)
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
    return new Adsorbate(tmin, tmax, pref, &coeffs[0]);
}

SpeciesThermoInterpType* newSpeciesThermoInterpType(const XML_Node& thermo)
{
    std::string model = toLowerCopy(thermo["model"]);
    if (model == "hkft" || model == "ionfromneutral") {
        // Some PDSS species use the 'thermo' node, but don't specify a
        // SpeciesThermoInterpType parameterization. This function needs to
        // just ignore this data.
        return 0;
    }

    // Get the children of the thermo XML node. In the next bit of code we take
    // out the comments that may have been children of the thermo XML node by
    // doing a selective copy. These shouldn't interfere with the algorithm at
    // any point.
    const std::vector<XML_Node*>& tpWC = thermo.children();
    std::vector<XML_Node*> tp;
    for (size_t i = 0; i < tpWC.size(); i++) {
        if (!tpWC[i]->isComment()) {
            tp.push_back(tpWC[i]);
        }
    }

    std::string thermoType = toLowerCopy(tp[0]->name());

    for (size_t i = 1; i < tp.size(); i++) {
        if (!caseInsensitiveEquals(tp[i]->name(), thermoType)) {
            throw CanteraError("newSpeciesThermoInterpType",
                "Encountered unsupported mixed species thermo "
                "parameterizations, '{}' and '{}'", tp[i]->name(), thermoType);
        }
    }
    if ((tp.size() > 2 && thermoType != "nasa9") ||
        (tp.size() > 1 && (thermoType == "const_cp" ||
                           thermoType == "mu0" ||
                           thermoType == "adsorbate"))) {
        throw CanteraError("newSpeciesThermoInterpType",
            "Too many regions in thermo parameterization.");
    }

    if (model == "mineraleq3") {
        if (thermoType != "mineq3") {
            throw CanteraError("newSpeciesThermoInterpType",
                               "confused: expected MinEQ3");
        }
        return newShomateForMineralEQ3(*tp[0]);
    } else if (thermoType == "shomate") {
        return newShomateThermoFromXML(tp);
    } else if (thermoType == "const_cp") {
        return newConstCpThermoFromXML(*tp[0]);
    } else if (thermoType == "nasa") {
        return newNasaThermoFromXML(tp);
    } else if (thermoType == "mu0") {
        return newMu0ThermoFromXML(*tp[0]);
    } else if (thermoType == "nasa9") {
        return newNasa9ThermoFromXML(tp);
    } else if (thermoType == "adsorbate") {
        return newAdsorbateThermoFromXML(*tp[0]);
    } else {
        throw CanteraError("newSpeciesThermoInterpType",
            "Unknown species thermo model '" + thermoType + "'.");
    }
}

}
