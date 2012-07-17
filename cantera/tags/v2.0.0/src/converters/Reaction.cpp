/**
 *  @file Reaction.cpp
 *
 */

// Copyright 2001  California Institute of Technology

#include "Reaction.h"
#include <iostream>
#include <stdio.h>

using namespace std;

namespace ckr
{

Reaction forwardReaction(const Reaction& rxn)
{
    Reaction r(rxn);
    r.isReversible = false;
    r.krev = RateCoeff();
    return r;
}

Reaction reverseReaction(const Reaction& rxn)
{
    Reaction r(rxn);
    if (rxn.isReversible && (r.krev.A > 0.0)) {
        r.isReversible = false;
        r.products = rxn.reactants;
        r.reactants = rxn.products;
        r.kf = rxn.krev;
        r.krev = RateCoeff();
    } else {
        r.reactants.clear();
        r.products.clear();
    }
    return r;
}

Reaction& Reaction::operator=(const Reaction& b)
{
    if (this == &b) {
        return *this;
    }
    type = b.type;
    reactants = b.reactants;
    fwdOrder = b.fwdOrder;
    products = b.products;
    thirdBody = b.thirdBody;
    e3b = b.e3b;
    kf = b.kf;
    kf_aux = b.kf_aux;
    krev = b.krev;
    duplicate = b.duplicate;
    falloffType = b.falloffType;
    falloffParameters = b.falloffParameters;
    otherAuxData = b.otherAuxData;
    isFalloffRxn = b.isFalloffRxn;
    isChemActRxn = b.isChemActRxn;
    isThreeBodyRxn = b.isThreeBodyRxn;
    isDuplicate = b.isDuplicate;
    lines = b.lines;
    number = b.number;
    comment = b.comment;
    return *this;
}

void Reaction::write(std::ostream& s) const
{
    int nl = static_cast<int>(lines.size());
    for (int nn = 0; nn < nl; nn++) {
        s << lines[nn] << std::endl;
    }
    //         int nr = reactants.size();
    //         int np = products.size();
    //         int n;
    //         double nu;
    //         for (n = 0; n < nr; n++) {
    //             nu = reactants[n].number;
    //             if (nu > 1) s << nu;
    //             s << reactants[n].name;
    //             if (n < nr - 1)
    //                 s << " + ";
    //             else if (isThreeBodyRxn)
    //                 s << " + " << thirdBody;
    //             else if (isFalloffRxn || isChemActRxn)
    //                 s << " (+ " << thirdBody << ")";
    //         }

    //         if (isReversible) s << " = ";
    //         else s << "=>";
    //         for (n = 0; n < np; n++) {
    //             nu = products[n].number;
    //             if (nu > 1) s << nu;
    //             s << products[n].name;
    //             if (n < np - 1)
    //                 s << " + ";
    //             else if (isThreeBodyRxn)
    //                 s << " + " << thirdBody;
    //             else if (isFalloffRxn || isChemActRxn)
    //                 s << " (+ " << thirdBody << ")";
    //         }
    //         char kfstr[100];
    //         sprintf(kfstr, "  %14.5g   %4.2g   %f", kf.A, kf.n, kf.E);
    //         s << kfstr << endl;
}


/**
 *  stoichiometric coefficient of species s in the reaction.  Negative
 *  for reactants, positive for products, and zero if the species does
 *  not participate in the reaction.
 */

double Reaction::stoichCoefficient(const std::string& s) const
{
    int k;
    int nr = static_cast<int>(reactants.size());
    for (k = 0; k < nr; k++)
        if (reactants[k].name == s) {
            return -reactants[k].number;
        }
    int np = static_cast<int>(products.size());
    for (k = 0; k < np; k++)
        if (products[k].name == s) {
            return products[k].number;
        }
    return 0.0;
}

/**
 *  used to find undeclared duplicate reactions.
 *  @todo could be made faster
 */
bool Reaction::operator==(const Reaction& r) const
{
    int nr = static_cast<int>(reactants.size());
    int np = static_cast<int>(products.size());
    if (int(r.reactants.size()) != nr ||
            int(r.products.size()) != np || r.thirdBody != thirdBody) {
        return false;
    }

    std::string nm;
    std::map<std::string, double> coeffs;
    for (int ir = 0; ir < nr; ir++) {
        coeffs[reactants[ir].name] = -reactants[ir].number;
    }
    for (int ip = 0; ip < np; ip++) {
        coeffs[products[ip].name] = products[ip].number;
    }
    for (int jr = 0; jr < nr; jr++) {
        nm = r.reactants[jr].name;
        if (coeffs[nm] == 0.0) {
            return false;
        }
        coeffs[nm] /= -r.reactants[jr].number;
    }
    for (int jp = 0; jp < np; jp++) {
        nm = r.products[jp].name;
        if (coeffs[nm] == 0.0) {
            return false;
        }
        coeffs[nm] /= products[jp].number;
    }
    int nc = static_cast<int>(coeffs.size());
    std::vector<double> ratios;
    getMapValues(coeffs, ratios);

    if (!isReversible && ratios[0] < 0.0) {
        return false;
    }

    for (int ic = 0; ic < nc; ic++) {
        if (ratios[ic] != ratios[0]) {
            return false;
        }
    }
    return true;
}

}













