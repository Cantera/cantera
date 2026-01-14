/**
 *  @file SpeciesThermoFactory.cpp
 *    Definitions for factory functions to build instances of classes that
 *    manage the standard-state thermodynamic properties of a set of species
 *    (see @ref spthermo);
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
#include "cantera/base/stringUtils.h"
#include "cantera/base/Units.h"

namespace Cantera
{

SpeciesThermoInterpType* newSpeciesThermoInterpType(int type, double tlow,
    double thigh, double pref, span<const double> coeffs)
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

SpeciesThermoInterpType* newSpeciesThermoInterpType(const string& stype,
    double tlow, double thigh, double pref, span<const double> coeffs)
{
    int itype = -1;
    string type = toLowerCopy(stype);
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

void setupSpeciesThermo(SpeciesThermoInterpType& thermo,
                        const AnyMap& node)
{
    double Pref = node.convert("reference-pressure", "Pa", OneAtm);
    thermo.setRefPressure(Pref);
    thermo.input() = node;
}

void setupNasaPoly(NasaPoly2& thermo, const AnyMap& node)
{
    setupSpeciesThermo(thermo, node);
    vector<double> Tranges = node.convertVector("temperature-ranges", "K", 2, 3);
    const auto& data = node["data"].asVector<vector<double>>(Tranges.size()-1);
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

void setupShomatePoly(ShomatePoly2& thermo, const AnyMap& node)
{
    setupSpeciesThermo(thermo, node);
    vector<double> Tranges = node.convertVector("temperature-ranges", "K", 2, 3);
    const auto& data = node["data"].asVector<vector<double>>(Tranges.size()-1);
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

void setupConstCp(ConstCpPoly& thermo, const AnyMap& node)
{
    setupSpeciesThermo(thermo, node);
    if (node.hasKey("T-min")) {
        thermo.setMinTemp(node.convert("T-min", "K"));
    }
    if (node.hasKey("T-max")) {
        thermo.setMaxTemp(node.convert("T-max", "K"));
    }
    double T0 = node.convert("T0", "K", 298.15);
    double h0 = node.convert("h0", "J/kmol", 0.0);
    double s0 = node.convert("s0", "J/kmol/K", 0.0);
    double cp0 = node.convert("cp0", "J/kmol/K", 0.0);
    thermo.setParameters(T0, h0, s0, cp0);
}

void setupNasa9Poly(Nasa9PolyMultiTempRegion& thermo, const AnyMap& node)
{
    setupSpeciesThermo(thermo, node);
    vector<double> Tranges = node.convertVector("temperature-ranges", "K", 2, 999);
    const auto& data = node["data"].asVector<vector<double>>(Tranges.size()-1);
    map<double, vector<double>> regions;
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
    if (node.hasKey("T-min")) {
        thermo.setMinTemp(node.convert("T-min", "K"));
    }
    if (node.hasKey("T-max")) {
        thermo.setMaxTemp(node.convert("T-max", "K"));
    }
    bool dimensionless = node.getBool("dimensionless", false);
    double h0 = node.convert("h0", "J/kmol", 0.0);
    map<double, double> T_mu;
    for (const auto& [T_str, mu] : node["data"]) {
        double T = node.units().convertTo(fpValueCheck(T_str), "K");
        if (dimensionless) {
            T_mu[T] = mu.asDouble() * GasConstant * T;
        } else {
            T_mu[T] = node.units().convert(mu, "J/kmol");
        }
    }
    thermo.setParameters(h0, T_mu);
}

unique_ptr<SpeciesThermoInterpType> newSpeciesThermo(const AnyMap& node)
{
    string model = node["model"].asString();
    if (model == "NASA7") {
        auto thermo = make_unique<NasaPoly2>();
        setupNasaPoly(*thermo, node);
        return thermo;
    } else if (model == "Shomate") {
        auto thermo = make_unique<ShomatePoly2>();
        setupShomatePoly(*thermo, node);
        return thermo;
    } else if (model == "NASA9") {
        auto thermo = make_unique<Nasa9PolyMultiTempRegion>();
        setupNasa9Poly(*thermo, node);
        return thermo;
    } else if (model == "constant-cp") {
        auto thermo = make_unique<ConstCpPoly>();
        setupConstCp(*thermo, node);
        return thermo;
    } else if (model == "piecewise-Gibbs") {
        auto thermo = make_unique<Mu0Poly>();
        setupMu0(*thermo, node);
        return thermo;
    } else {
        throw CanteraError("newSpeciesThermo",
            "Unknown thermo model '{}'", model);
    }
}

}
