/**
 *  @file SpeciesThermoFactory.cpp
 *    Definitions for factory functions to build instances of classes that
 *    manage the standard-state thermodynamic properties of a set of species
 *    (see \ref spthermo);
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/thermo/MultiSpeciesThermo.h"
#include "cantera/thermo/Mu0Poly.h"
#include "cantera/thermo/Nasa9PolyMultiTempRegion.h"
#include "cantera/thermo/Nasa9Poly1.h"
#include "cantera/thermo/NasaPoly2.h"
#include "cantera/thermo/ShomatePoly.h"
#include "cantera/thermo/ConstCpPoly.h"
#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/thermo/VPStandardStateTP.h"
#include "cantera/base/ctml.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/Units.h"

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
    case NASA9MULTITEMP:
        return new Nasa9PolyMultiTempRegion(tlow, thigh, pref, coeffs);
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
        throw CanteraError("newNasaThermoFromXML",
                           "non-continuous temperature ranges.");
    }

    vector_fp c(15);
    c[0] = tmid;
    copy(c1.begin(), c1.begin()+7, c.begin() + 1); // high-T coefficients
    copy(c0.begin(), c0.begin()+7, c.begin() + 8); // low-T coefficients
    return newSpeciesThermoInterpType(NASA, tmin, tmax, p0, &c[0]);
}

void setupSpeciesThermo(SpeciesThermoInterpType& thermo,
                        const AnyMap& node)
{
    double Pref = node.convert("reference-pressure", "Pa", OneAtm);
    thermo.setRefPressure(Pref);
}

void setupNasaPoly(NasaPoly2& thermo, const AnyMap& node)
{
    setupSpeciesThermo(thermo, node);
    vector_fp Tranges = node.convertVector("temperature-ranges", "K", 2, 3);
    const auto& data = node["data"].asVector<vector_fp>(Tranges.size()-1);
    for (const auto& poly : data) {
        if (poly.size() != 7) {
            throw CanteraError("setupNasaPoly", "Wrong number of coefficients "
                "for NASA polynomial. Expected 7, but got {}", poly.size());
        }
    }
    thermo.setMinTemp(Tranges.front());
    thermo.setMaxTemp(Tranges.back());
    if (Tranges.size() == 3) { // standard 2 temperature range polynomial
        thermo.setParameters(Tranges[1], data[0], data[1]);
    } else { // Repeat data for single temperature range for both ranges
        thermo.setParameters(Tranges[1], data[0], data[0]);
    }
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
              throw CanteraError("newShomateThermoFromXML",
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
        throw CanteraError("newShomateThermoFromXML",
                           "non-continuous temperature ranges.");
    }
    if(c0.size() != 7 || c1.size() != 7)
    {
      throw CanteraError("newShomateThermoFromXML",
                         "Shomate thermo requires 7 coefficients in float array.");
    }
    vector_fp c(15);
    c[0] = tmid;
    copy(c0.begin(), c0.begin()+7, c.begin() + 1);
    copy(c1.begin(), c1.begin()+7, c.begin() + 8);
    return newSpeciesThermoInterpType(SHOMATE, tmin, tmax, p0, &c[0]);
}


void setupShomatePoly(ShomatePoly2& thermo, const AnyMap& node)
{
    setupSpeciesThermo(thermo, node);
    vector_fp Tranges = node.convertVector("temperature-ranges", "K", 2, 3);
    const auto& data = node["data"].asVector<vector_fp>(Tranges.size()-1);
    for (const auto& poly : data) {
        if (poly.size() != 7) {
            throw CanteraError("setupShomatePoly", "Wrong number of coefficients "
                "for Shomate polynomial. Expected 7, but got {}", poly.size());
        }
    }
    thermo.setMinTemp(Tranges.front());
    thermo.setMaxTemp(Tranges.back());
    if (Tranges.size() == 3) { // standard 2 temperature range polynomial
        thermo.setParameters(Tranges[1], data[0], data[1]);
    } else { // Repeat data for single temperature range for both ranges
        thermo.setParameters(Tranges[1], data[0], data[0]);
    }
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

void setupConstCp(ConstCpPoly& thermo, const AnyMap& node)
{
    setupSpeciesThermo(thermo, node);
    double T0 = node.convert("T0", "K", 298.15);
    double h0 = node.convert("h0", "J/kmol", 0.0);
    double s0 = node.convert("s0", "J/kmol/K", 0.0);
    double cp0 = node.convert("cp0", "J/kmol/K", 0.0);
    thermo.setParameters(T0, h0, s0, cp0);
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
                throw CanteraError("newNasa9ThermoFromXML",
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


void setupNasa9Poly(Nasa9PolyMultiTempRegion& thermo, const AnyMap& node)
{
    setupSpeciesThermo(thermo, node);
    vector_fp Tranges = node.convertVector("temperature-ranges", "K", 2, 999);
    const auto& data = node["data"].asVector<vector_fp>(Tranges.size()-1);
    map<double, vector_fp> regions;
    for (size_t i = 0; i < data.size(); i++) {
        if (data[i].size() != 9) {
            throw CanteraError("setupNasa9Poly", "Wrong number of coefficients "
                "for NASA9 polynomial. Expected 9, but got {}", data[i].size());
        }
        regions[Tranges[i]] = data[i];
    }
    thermo.setMinTemp(Tranges.front());
    thermo.setMaxTemp(Tranges.back());
    thermo.setParameters(regions);
}


void setupMu0(Mu0Poly& thermo, const AnyMap& node)
{
    setupSpeciesThermo(thermo, node);
    bool dimensionless = node.getBool("dimensionless", false);
    double h0 = node.convert("h0", "J/kmol", 0.0);
    map<double, double> T_mu;
    for (const auto& item : node["data"]) {
        double T = node.units().convert(fpValueCheck(item.first), "K");
        if (dimensionless) {
            T_mu[T] = item.second.asDouble() * GasConstant * T;
        } else {
            T_mu[T] = node.units().convert(item.second, "J/kmol");
        }
    }
    thermo.setParameters(h0, T_mu);
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
                           thermoType == "mu0"))) {
        throw CanteraError("newSpeciesThermoInterpType",
            "Too many regions in thermo parameterization.");
    }

    if (thermoType == "shomate") {
        return newShomateThermoFromXML(tp);
    } else if (thermoType == "const_cp") {
        return newConstCpThermoFromXML(*tp[0]);
    } else if (thermoType == "nasa") {
        return newNasaThermoFromXML(tp);
    } else if (thermoType == "mu0") {
        return newMu0ThermoFromXML(*tp[0]);
    } else if (thermoType == "nasa9") {
        return newNasa9ThermoFromXML(tp);
    } else {
        throw CanteraError("newSpeciesThermoInterpType",
            "Unknown species thermo model '" + thermoType + "'.");
    }
}


unique_ptr<SpeciesThermoInterpType> newSpeciesThermo(const AnyMap& node)
{
    std::string model = node["model"].asString();
    if (model == "NASA7") {
        unique_ptr<NasaPoly2> thermo(new NasaPoly2());
        setupNasaPoly(*thermo, node);
        return unique_ptr<SpeciesThermoInterpType>(move(thermo));
    } else if (model == "Shomate") {
        unique_ptr<ShomatePoly2> thermo(new ShomatePoly2());
        setupShomatePoly(*thermo, node);
        return unique_ptr<SpeciesThermoInterpType>(move(thermo));
    } else if (model == "NASA9") {
        unique_ptr<Nasa9PolyMultiTempRegion> thermo(new Nasa9PolyMultiTempRegion());
        setupNasa9Poly(*thermo, node);
        return unique_ptr<SpeciesThermoInterpType>(move(thermo));
    } else if (model == "constant-cp") {
        unique_ptr<ConstCpPoly> thermo(new ConstCpPoly());
        setupConstCp(*thermo, node);
        return unique_ptr<SpeciesThermoInterpType>(move(thermo));
    } else if (model == "piecewise-Gibbs") {
        unique_ptr<Mu0Poly> thermo(new Mu0Poly());
        setupMu0(*thermo, node);
        return unique_ptr<SpeciesThermoInterpType>(move(thermo));
    } else {
        throw CanteraError("newSpeciesThermo",
            "Unknown thermo model '{}'", model);
    }
}

}
