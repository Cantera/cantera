/**
 *  @file importKinetics.cpp
 *     Declarations of global routines for the importing
 *     of kinetics data from XML files (see \ref inputfiles).
 *
 *     This file contains routines which are global routines, i.e.,
 *     not part of any object. These routine take as input, ctml
 *     pointers to data, and pointers to %Cantera objects. The purpose
 *     of these routines is to initialize the %Cantera objects with data
 *     from the ctml tree structures.
 */
// Copyright 2002  California Institute of Technology

#include "cantera/kinetics/importKinetics.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/kinetics/ReactionData.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"

#include <cstring>

using namespace std;

namespace Cantera
{

ReactionRules::ReactionRules() :
    skipUndeclaredSpecies(false),
    skipUndeclaredThirdBodies(false),
    allowNegativeA(false)
{
}

void checkRxnElementBalance(Kinetics& kin,
                            const ReactionData& rdata, doublereal errorTolerance)
{
    warn_deprecated("checkRxnElementBalance", "Now handled by "
        "Kinetics::checkReactionBalance. To be removed after Cantera 2.2.");
    doublereal kstoich;

    map<string, double> bal, balr, balp;
    bal.clear();
    balp.clear();
    balr.clear();
    size_t np = rdata.products.size();

    // iterate over the products
    for (size_t index = 0; index < np; index++) {
        size_t kp = rdata.products[index];     // index of the product in 'kin'
        size_t n = kin.speciesPhaseIndex(kp);  // phase this product belongs to
        size_t klocal = kp - kin.kineticsSpeciesIndex(0,n); // index within this phase
        kstoich = rdata.pstoich[index]; // product stoichiometric coeff
        const ThermoPhase& ph = kin.speciesPhase(kp);
        for (size_t m = 0; m < ph.nElements(); m++) {
            bal[ph.elementName(m)] += kstoich*ph.nAtoms(klocal,m);
            balp[ph.elementName(m)] += kstoich*ph.nAtoms(klocal,m);
        }
    }
    for (size_t index = 0; index < rdata.reactants.size(); index++) {
        size_t kr = rdata.reactants[index];
        size_t n = kin.speciesPhaseIndex(kr);
        size_t klocal = kr - kin.kineticsSpeciesIndex(0,n);
        kstoich = rdata.rstoich[index];
        const ThermoPhase& ph = kin.speciesPhase(kr);
        for (size_t m = 0; m < ph.nElements(); m++) {
            bal[ph.elementName(m)] -= kstoich*ph.nAtoms(klocal,m);
            balr[ph.elementName(m)] += kstoich*ph.nAtoms(klocal,m);
        }
    }

    map<string, double>::iterator b = bal.begin();
    string msg = "\n\tElement    Reactants    Products";
    bool ok = true;
    doublereal err, elemsum;
    for (; b != bal.end(); ++b) {
        elemsum = fabs(balr[b->first]) + fabs(balp[b->first]);
        if (elemsum > 0.0) {
            err = fabs(b->second/elemsum);
            if (err > errorTolerance) {
                ok = false;
                msg += "\n\t"+b->first+"           "+ fp2str(balr[b->first])
                       +"           "+ fp2str(balp[b->first]);
            }
        }
    }
    if (!ok) {
        msg = "The following reaction is unbalanced:\n\t"
              + rdata.equation + "\n" + msg + "\n";
        throw CanteraError("checkRxnElementBalance",msg);
    }
}

bool getReagents(const XML_Node& rxn, Kinetics& kin, int rp,
                 std::string default_phase, std::vector<size_t>& spnum,
                 vector_fp& stoich, vector_fp& order,
                 const ReactionRules& rules)
{
    warn_deprecated("getReagents", "Now handled through newReaction() and its "
        "support functions. To be removed after Cantera 2.2.");
    string rptype;

    /*
     * The id of reactants and products are kept in child elements
     * of reaction, named "reactants" and "products". We search
     * the XML tree for these children based on the value of rp,
     * and store the XML element pointer here.
     */
    if (rp == 1) {
        rptype = "reactants";
    } else {
        rptype = "products";
    }
    const XML_Node& rg = rxn.child(rptype);

    /*
     * The species and stoichiometric coefficient for the species
     * are stored as a colon separated pair. Get all of these
     * pairs in the reactions/products object.
     */
    std::vector<string> key, val;
    getPairs(rg, key, val);

    /*
     * Loop over each of the pairs and process them
     */
    doublereal ord, stch;
    string ph, spName;
    map<string, size_t> speciesMap;
    for (size_t n = 0; n < key.size(); n++) {
        spName = key[n]; // sp is the string name for species
        ph = "";
        /*
         * Search for the species in the kinetics object using the
         * member function kineticsSpeciesIndex(). We will search
         * for the species in all phases defined in the kinetics operator.
         */
        size_t isp = kin.kineticsSpeciesIndex(spName);
        if (isp == npos) {
            if (rules.skipUndeclaredSpecies) {
                return false;
            } else {
                throw CanteraError("getReagents",
                                   "Undeclared reactant or product species " + spName);
                return false;
            }
        }

        /*
         * For each reagent, we store the the species number, isp
         * the stoichiometric coefficient, val[n], and the order
         * species in the reaction rate expression. We assume mass
         * action kinetics here, but will modify this below for
         * specified species.
         */
        spnum.push_back(isp);
        stch = fpValue(val[n]);
        stoich.push_back(stch);
        ord = doublereal(stch);
        order.push_back(ord);

        /*
         * Needed to process reaction orders below.
         */
        speciesMap[spName] = order.size();
    }

    /*
     * Check to see if reaction orders have been specified.
     */

    if (rp == 1 && rxn.hasChild("order")) {
        std::vector<XML_Node*> ord = rxn.getChildren("order");
        doublereal forder;
        for (size_t nn = 0; nn < ord.size(); nn++) {
            const XML_Node& oo = *ord[nn];
            string sp = oo["species"];
            size_t loc = speciesMap[sp];
            if (loc == 0)
                throw CanteraError("getReagents",
                                   "reaction order specified for non-reactant: "
                                   +sp);
            forder = oo.fp_value();
            if (forder < 0.0) {
                throw CanteraError("getReagents",
                                   "reaction order must be non-negative");
            }
            // replace the stoichiometric coefficient
            // stored above in 'order' with the specified
            // reaction order
            order[loc-1] = forder;
        }
    }
    return true;
}

/**
 * getArrhenius() parses the XML element called Arrhenius.
 * The Arrhenius expression is
 * \f[        k =  A T^(b) exp (-E_a / RT). \f]
 * @deprecated to be removed after Cantera 2.2.
 */
static void getArrhenius(const XML_Node& node, int& labeled,
                         doublereal& A, doublereal& b, doublereal& E)
{
    if (node["name"] == "k0") {
        labeled = -1;
    } else if (node["name"] == "kHigh") {
        labeled = 1;
    } else {
        labeled = 0;
    }
    /*
     * We parse the children for the A, b, and E components.
     */
    A = getFloat(node, "A", "toSI");
    b = getFloat(node, "b");
    E = getFloat(node, "E", "actEnergy");
    E /= GasConstant;
}

/**
 * getStick() processes the XML element called Stick that specifies
 * the sticking coefficient reaction. This routine will
 * translate the sticking coefficient value into a "normal"
 * rate constant for the surface reaction.
 *
 *  Output
 * -----------
 * Output is the normal Arrhenius expressions for a surface
 * reaction rate constant.
 *
 *   A - units such that rate of rxn has kmol/m^2/s when
 *       A is multiplied by activity concentrations of
 *       reactants in the normal manner.
 *   n - unitless
 *   E - Units 1/Kelvin
 * @deprecated to be removed after Cantera 2.2.
 */
static void getStick(const XML_Node& node, Kinetics& kin,
                     ReactionData& r, doublereal& A, doublereal& b, doublereal& E)
{
    size_t nr = r.reactants.size();
    size_t k, klocal, not_surf = 0;
    size_t np = 0;
    doublereal f = 1.0;
    doublereal order;
    /*
     * species is the name of the special reactant whose surface
     * flux rate will be calculated.
     *      isp = species # in the local phase
     *      ispKinetics = species # in the kinetics object
     *      ispPhaseIndex = phase # of the special species
     */
    string spname = node["species"];
    ThermoPhase& th = kin.speciesPhase(spname);
    size_t isp = th.speciesIndex(spname);
    size_t ispKinetics = kin.kineticsSpeciesIndex(spname);
    size_t ispPhaseIndex = kin.speciesPhaseIndex(ispKinetics);

    doublereal ispMW = th.molecularWeights()[isp];
    doublereal sc;

    // loop over the reactants
    for (size_t n = 0; n < nr; n++) {
        k = r.reactants[n];
        order = r.rorder[n];    // stoich coeff

        // get the phase species k belongs to
        np = kin.speciesPhaseIndex(k);
        const ThermoPhase& p = kin.thermo(np);

        // get the local index of species k in this phase
        klocal = p.speciesIndex(kin.kineticsSpeciesName(k));

        // if it is a surface species, divide f by the standard
        // concentration for this species, in order to convert
        // from concentration units used in the law of mass action
        // to coverages used in the sticking probability expression
        if (p.eosType() == cSurf || p.eosType() == cEdge) {
            sc = p.standardConcentration(klocal);
            f /= pow(sc, order);
        }
        // Otherwise:
        else {
            // We only allow one species to be in the phase containing the
            // special sticking coefficient species.
            if (ispPhaseIndex == np) {
                not_surf++;
            }
            // Other bulk phase species on the other side of ther interface are
            // treated like surface species.
            else {
                sc = p.standardConcentration(klocal);
                f /= pow(sc, order);
            }
        }
    }
    if (not_surf != 1) {
        throw CanteraError("getStick",
                           "reaction probabilities can only be used in "
                           "reactions with exactly 1 gas/liquid species.");
    }

    doublereal cbar = sqrt(8.0*GasConstant/(Pi*ispMW));
    A = 0.25 * getFloat(node, "A", "toSI") * cbar * f;
    b = getFloat(node, "b") + 0.5;
    E = getFloat(node, "E", "actEnergy");
    E /= GasConstant;
}

//!  Read the XML data concerning the coverage dependence of an interfacial reaction
/*!
 *     @param node       XML node with name reaction containing the reaction information
 *     @param surfphase  Surface phase
 *     @param rdata      Reaction data for the reaction.
 *
 *  Example:
 * @verbatim
        <coverage species="CH3*">
           <a> 1.0E-5 </a>
           <m> 0.0 </m>
           <actEnergy> 0.0 </actEnergy>
        </coverage>
@endverbatim
 * @deprecated to be removed after Cantera 2.2.
 */
static void getCoverageDependence(const XML_Node& node,
                                  thermo_t& surfphase, ReactionData& rdata)
{
    vector<XML_Node*> cov = node.getChildren("coverage");
    size_t k, nc = cov.size();
    doublereal e;
    string spname;
    if (nc > 0) {
        for (size_t n = 0; n < nc; n++) {
            const XML_Node& cnode = *cov[n];
            spname = cnode["species"];
            k = surfphase.speciesIndex(spname);
            rdata.cov.push_back(doublereal(k));
            rdata.cov.push_back(getFloat(cnode, "a"));
            rdata.cov.push_back(getFloat(cnode, "m"));
            e = getFloat(cnode, "e", "actEnergy");
            rdata.cov.push_back(e/GasConstant);
        }
    }
}

//! Get falloff parameters for a reaction.
/*!
 *  This routine reads the falloff XML node and extracts parameters into a
 *  vector of doubles
 *
 * @verbatim
 <falloff type="Troe"> 0.5 73.2 5000. 9999. </falloff>
 @endverbatim
 * @deprecated to be removed after Cantera 2.2.
*/
static void getFalloff(const XML_Node& f, ReactionData& rdata)
{
    string type = f["type"];
    vector<string> p;
    getStringArray(f,p);
    vector_fp c;
    size_t np = p.size();
    for (size_t n = 0; n < np; n++) {
        c.push_back(fpValue(p[n]));
    }
    if (type == "Troe") {
        if (np == 3 || np == 4) {
            rdata.falloffType = TROE_FALLOFF;
        } else {
            throw CanteraError("getFalloff()", "Troe parameterization is specified by number of parameters, "
                               + int2str(np) + ", is not equal to 3 or 4");
        }
    } else if (type == "SRI") {
        if (np == 3 || np == 5) {
            rdata.falloffType = SRI_FALLOFF;
        } else {
            throw CanteraError("getFalloff()", "SRI parameterization is specified by number of parameters, "
                               + int2str(np) + ", is not equal to 3 or 5");
        }
    }
    rdata.falloffParameters = c;
}

/**
 * Get the enhanced collision efficiencies. It is assumed that the
 * reaction mechanism is homogeneous, so that all species belong
 * to phase(0) of 'kin'.
 * @deprecated to be removed after Cantera 2.2.
 */
static void getEfficiencies(const XML_Node& eff, Kinetics& kin,
                            ReactionData& rdata, const ReactionRules& rules)
{
    // set the default collision efficiency
    rdata.default_3b_eff = fpValue(eff["default"]);

    vector<string> key, val;
    getPairs(eff, key, val);
    string nm;
    string phse = kin.thermo(0).id();
    for (size_t n = 0; n < key.size(); n++) {
        nm = key[n];
        size_t k = kin.kineticsSpeciesIndex(nm, phse);
        if (k != npos) {
            rdata.thirdBodyEfficiencies[k] = fpValue(val[n]);
        } else if (!rules.skipUndeclaredThirdBodies) {
            throw CanteraError("getEfficiencies", "Encountered third-body "
                               "efficiency for undefined species \"" + nm + "\"\n"
                               "while adding reaction " + int2str(rdata.number+1) + ".");
        }
    }
}

void getRateCoefficient(const XML_Node& kf, Kinetics& kin,
                        ReactionData& rdata, const ReactionRules& rules)
{
    warn_deprecated("getRateCoefficent", "Now handled through newReaction() "
        "and its support functions. To be removed after Cantera 2.2.");
    if (rdata.reactionType == PLOG_RXN) {
        rdata.rateCoeffType = PLOG_REACTION_RATECOEFF_TYPE;
        for (size_t m = 0; m < kf.nChildren(); m++) {
            const XML_Node& node = kf.child(m);
            double p = getFloat(node, "P", "toSI");
            vector_fp& rate = rdata.plogParameters.insert(
                                  std::make_pair(p, vector_fp()))->second;
            rate.resize(3);
            rate[0] = getFloat(node, "A", "toSI");
            rate[1] = getFloat(node, "b");
            rate[2] = getFloat(node, "E", "actEnergy") / GasConstant;
        }

    } else if (rdata.reactionType == CHEBYSHEV_RXN) {
        rdata.rateCoeffType = CHEBYSHEV_REACTION_RATECOEFF_TYPE;
        rdata.chebTmin = getFloat(kf, "Tmin", "toSI");
        rdata.chebTmax = getFloat(kf, "Tmax", "toSI");
        rdata.chebPmin = getFloat(kf, "Pmin", "toSI");
        rdata.chebPmax = getFloat(kf, "Pmax", "toSI");
        const XML_Node& coeffs = kf.child("floatArray");
        rdata.chebDegreeP = atoi(coeffs["degreeP"].c_str());
        rdata.chebDegreeT = atoi(coeffs["degreeT"].c_str());
        getFloatArray(kf, rdata.chebCoeffs, false);

    } else {

        string type = kf.attrib("type");
        if (type == "") {
            type = "Arrhenius";
            rdata.rateCoeffType = ARRHENIUS_REACTION_RATECOEFF_TYPE;
        }
        if (type == "ExchangeCurrentDensity") {
            rdata.rateCoeffType = EXCHANGE_CURRENT_REACTION_RATECOEFF_TYPE;
        } else if (type == "Arrhenius") {

        } else {
            throw CanteraError("getRateCoefficient", "Unknown type: " + type);
        }

        vector_fp c_alt(3,0.0), c_base(3,0.0);
        for (size_t m = 0; m < kf.nChildren(); m++) {
            const XML_Node& c = kf.child(m);
            string nm = c.name();
            int labeled=0;

            if (nm == "Arrhenius") {
                vector_fp coeff(3);
                if (c["type"] == "stick") {
                    getStick(c, kin, rdata, coeff[0], coeff[1], coeff[2]);
                    c_base = coeff;
                } else {
                    getArrhenius(c, labeled, coeff[0], coeff[1], coeff[2]);
                    if (labeled == 0 || rdata.reactionType == THREE_BODY_RXN
                            || rdata.reactionType == ELEMENTARY_RXN) {
                        c_base = coeff;
                    } else {
                        c_alt = coeff;
                    }
                }
                if (rdata.reactionType == SURFACE_RXN || rdata.reactionType == EDGE_RXN) {
                    getCoverageDependence(c,
                                          kin.thermo(kin.surfacePhaseIndex()), rdata);
                }

                if (coeff[0] < 0.0 && !rules.allowNegativeA) {
                    throw CanteraError("getRateCoefficient",
                                       "negative A coefficient for reaction "+int2str(rdata.number));
                }
            } else if (nm == "Arrhenius_ExchangeCurrentDensity") {
                vector_fp coeff(3);
                getArrhenius(c, labeled, coeff[0], coeff[1], coeff[2]);
                c_base = coeff;
                rdata.rateCoeffType = EXCHANGE_CURRENT_REACTION_RATECOEFF_TYPE;
            } else if (nm == "falloff") {
                getFalloff(c, rdata);
            } else if (nm == "efficiencies") {
                getEfficiencies(c, kin, rdata, rules);
            } else if (nm == "electrochem") {
                rdata.beta = fpValue(c["beta"]);
            }
        }
        /*
         * Store the coefficients in the ReactionData object for return
         * from this function.
         */
        if (rdata.reactionType == FALLOFF_RXN) {
            rdata.rateCoeffParameters = c_base;
            rdata.auxRateCoeffParameters = c_alt;
        } else if (rdata.reactionType == CHEMACT_RXN) {
            rdata.rateCoeffParameters = c_alt;
            rdata.auxRateCoeffParameters = c_base;
        } else {
            rdata.rateCoeffParameters = c_base;
        }

    }
}

doublereal isDuplicateReaction(std::map<int, doublereal>& r1,
                               std::map<int, doublereal>& r2)
{
    warn_deprecated("isDuplicateReaction", "Now handled by "
        "Kinetics::checkDuplicateStoich. To be removed after Cantera 2.2.");
    map<int, doublereal>::const_iterator b = r1.begin(), e = r1.end();
    int k1 = b->first;
    // check for duplicate written in the same direction
    doublereal ratio = 0.0;
    if (r1[k1] && r2[k1]) {
        ratio = r2[k1]/r1[k1];
        ++b;
        bool different = false;
        for (; b != e; ++b) {
            k1 = b->first;
            if (!r1[k1] || !r2[k1] || fabs(r2[k1]/r1[k1] - ratio) > 1.e-8) {
                different = true;
                break;
            }
        }
        if (!different) {
            return ratio;
        }
    }

    // check for duplicate written in the reverse direction
    b = r1.begin();
    k1 = b->first;
    if (r1[k1] == 0.0 || r2[-k1] == 0.0) {
        return 0.0;
    }
    ratio = r2[-k1]/r1[k1];
    ++b;
    for (; b != e; ++b) {
        k1 = b->first;
        if (!r1[k1] || !r2[-k1] || fabs(r2[-k1]/r1[k1] - ratio) > 1.e-8) {
            return 0.0;
        }
    }
    return ratio;
}

bool installReactionArrays(const XML_Node& p, Kinetics& kin,
                           std::string default_phase, bool check_for_duplicates)
{
    int itot = 0;
    /*
     * Search the children of the phase element for the
     * XML element named reactionArray. If we can't find it,
     * then return signaling having not found any reactions.
     * Apparently, we allow multiple reactionArray elements here
     * Each one will be processed sequentially, with the
     * end result being purely additive.
     */
    vector<XML_Node*> rarrays = p.getChildren("reactionArray");
    if (rarrays.empty()) {
        kin.finalize();
        return false;
    }
    for (size_t n = 0; n < rarrays.size(); n++) {
        /*
         * Go get a reference to the current XML element,
         * reactionArray. We will process this element now.
         */
        const XML_Node& rxns = *rarrays[n];
        /*
         * The reactionArray element has an attribute called,
         * datasrc. The value of the attribute is the XML
         * element comprising the top of the
         * tree of reactions for the phase.
         * Find this datasrc element starting with the root
         * of the current XML node.
         */
        const XML_Node* rdata = get_XML_Node(rxns["datasrc"], &rxns.root());
        /*
         * If the reactionArray element has a child element named "skip", and
         * if the attribute of skip called "species" has a value of "undeclared",
         * we will set rxnrule.skipUndeclaredSpecies to 'true'. rxnrule is
         * passed to the routine that parses each individual reaction so that
         * the parser will skip all reactions containing an undefined species
         * without throwing an error.
         *
         * Similarly, an attribute named "third_bodies" with the value of
         * "undeclared" will skip undeclared third body efficiencies (while
         * retaining the reaction and any other efficiencies).
         */
        if (rxns.hasChild("skip")) {
            const XML_Node& sk = rxns.child("skip");
            if (sk["species"] == "undeclared") {
                kin.skipUndeclaredSpecies(true);
            }
            if (sk["third_bodies"] == "undeclared") {
                kin.skipUndeclaredThirdBodies(true);
            }
        }
        /*
         * Search for child elements called include. We only include
         * a reaction if it's tagged by one of the include fields.
         * Or, we include all reactions if there are no include fields.
         */
        vector<XML_Node*> incl = rxns.getChildren("include");
        vector<XML_Node*> allrxns = rdata->getChildren("reaction");
        // if no 'include' directive, then include all reactions
        if (incl.empty()) {
            for (size_t i = 0; i < allrxns.size(); i++) {
                kin.addReaction(newReaction(*allrxns[i]));
                ++itot;
            }
        } else {
            for (size_t nii = 0; nii < incl.size(); nii++) {
                const XML_Node& ii = *incl[nii];
                string imin = ii["min"];
                string imax = ii["max"];

                string::size_type iwild = string::npos;
                if (imax == imin) {
                    iwild = imin.find("*");
                    if (iwild != string::npos) {
                        imin = imin.substr(0,iwild);
                        imax = imin;
                    }
                }

                for (size_t i = 0; i < allrxns.size(); i++) {
                    const XML_Node* r = allrxns[i];
                    string rxid;
                    if (r) {
                        rxid = r->attrib("id");
                        if (iwild != string::npos) {
                            rxid = rxid.substr(0,iwild);
                        }
                        /*
                         * To decide whether the reaction is included or not
                         * we do a lexical min max and operation. This
                         * sometimes has surprising results.
                         */
                        if ((rxid >= imin) && (rxid <= imax)) {
                            kin.addReaction(newReaction(*r));
                            ++itot;
                        }
                    }
                }
            }
        }
    }

    if (check_for_duplicates) {
        kin.checkDuplicates();
    }
    /*
     * Finalize the installation of the kinetics, now that we know
     * the true number of reactions in the mechanism, itot.
     */
    kin.finalize();

    return true;
}

bool importKinetics(const XML_Node& phase, std::vector<ThermoPhase*> th,
                    Kinetics* k)
{
    if (k == 0) {
        return false;
    }

    // This phase will be the owning phase for the kinetics operator
    // For interfaces, it is the surface phase between two volumes.
    // For homogeneous kinetics, it's the current volumetric phase.
    string owning_phase = phase["id"];

    bool check_for_duplicates = false;
    if (phase.parent()) {
        if (phase.parent()->hasChild("validate")) {
            const XML_Node& d = phase.parent()->child("validate");
            if (d["reactions"] == "yes") {
                check_for_duplicates = true;
            }
        }
    }

    // if other phases are involved in the reaction mechanism,
    // they must be listed in a 'phaseArray' child
    // element. Homogeneous mechanisms do not need to include a
    // phaseArray element.

    vector<string> phase_ids;
    if (phase.hasChild("phaseArray")) {
        const XML_Node& pa = phase.child("phaseArray");
        getStringArray(pa, phase_ids);
    }
    phase_ids.push_back(owning_phase);

    int np = static_cast<int>(phase_ids.size());
    int nt = static_cast<int>(th.size());

    // for each referenced phase, attempt to find its id among those
    // phases specified.
    bool phase_ok;

    string phase_id;
    string msg = "";
    for (int n = 0; n < np; n++) {
        phase_id = phase_ids[n];
        phase_ok = false;

        // loop over the supplied 'ThermoPhase' objects representing
        // phases, to find an object with the same id.
        for (int m = 0; m < nt; m++) {
            if (th[m]->id() == phase_id) {
                phase_ok = true;

                // if no phase with this id has been added to
                //the kinetics manager yet, then add this one
                if (k->phaseIndex(phase_id) == npos) {
                    k->addPhase(*th[m]);
                }
            }
            msg += " "+th[m]->id();
        }
        if (!phase_ok) {
            throw CanteraError("importKinetics",
                               "phase "+phase_id+" not found. Supplied phases are:"+msg);
        }
    }

    // allocates arrays, etc. Must be called after the phases have
    // been added to 'kin', so that the number of species in each
    // phase is known.
    k->init();

    // Install the reactions.
    return installReactionArrays(phase, *k, owning_phase, check_for_duplicates);
}

bool buildSolutionFromXML(XML_Node& root, const std::string& id,
                          const std::string& nm, ThermoPhase* th, Kinetics* kin)
{
    XML_Node* x;
    x = get_XML_NameID(nm, string("#")+id, &root);
    if (!x) {
        return false;
    }

    /*
     * Fill in the ThermoPhase object by querying the
     * const XML_Node tree located at x.
     */
    importPhase(*x, th);
    /*
     * Create a vector of ThermoPhase pointers of length 1
     * having the current th ThermoPhase as the entry.
     */
    std::vector<ThermoPhase*> phases(1);
    phases[0] = th;
    /*
     * Fill in the kinetics object k, by querying the
     * const XML_Node tree located by x. The source terms and
     * eventually the source term vector will be constructed
     * from the list of ThermoPhases in the vector, phases.
     */
    importKinetics(*x, phases, kin);
    return true;
}

}
