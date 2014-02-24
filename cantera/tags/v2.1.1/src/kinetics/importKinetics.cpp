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
#include "cantera/thermo/mix_defs.h"

//   Cantera includes
#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/EdgePhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/kinetics/reaction_defs.h"
#include "cantera/kinetics/ReactionData.h"
#include "cantera/base/global.h"
#include "cantera/base/stringUtils.h"

#include "cantera/base/xml.h"
#include "cantera/base/ctml.h"

using namespace ctml;
using namespace std;

namespace Cantera
{

ReactionRules::ReactionRules() :
    skipUndeclaredSpecies(false),
    skipUndeclaredThirdBodies(false),
    allowNegativeA(false)
{
}

//! these are all used to check for duplicate reactions
class rxninfo
{
public:
    std::vector<ReactionData*> m_rdata;

    //! Map of (key indicating participating species) to reaction numbers
    //! Used to speed up duplicate reaction checks.
    std::map<unsigned long int, std::vector<size_t> > m_participants;

    /**
     *  Install an individual reaction into a kinetics manager. The
     *  data for the reaction is in the xml_node r. In other words, r
     *  points directly to a ctml element named "reaction". i refers
     *  to the number id of the reaction in the kinetics object.
     *
     * @param i Reaction number.
     * @param r XML_Node containing reaction data.
     * @param kin Kinetics manager to which reaction will be added.
     * @param default_phase Default phase for locating a species
     * @param rules Rule for handling reactions with missing species
     *             (skip or flag as error)
     * @param validate_rxn If true, check that this reaction is not a
     *                     duplicate of one already entered, and check that
     *                     the reaction balances.
     *
     * @ingroup kineticsmgr
     */
    bool installReaction(int i, const XML_Node& r, Kinetics& kin,
                         std::string default_phase, ReactionRules& rules,
                         bool validate_rxn) ;

    ~rxninfo() {
        for (size_t i = 0; i < m_rdata.size(); i++) {
            delete m_rdata[i];
        }
    }
};

void checkRxnElementBalance(Kinetics& kin,
                            const ReactionData& rdata, doublereal errorTolerance)
{
    doublereal kstoich;

    map<string, double> bal, balr, balp;
    bal.clear();
    balp.clear();
    balr.clear();
    //cout << "checking " << rdata.equation << endl;
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
            //cout << "product species " << ph.speciesName(klocal) << " has " << ph.nAtoms(klocal,m)
            //     << " atoms of " << ph.elementName(m) << " and kstoich = " << kstoich << endl;
        }
    }
    for (size_t index = 0; index < rdata.reactants.size(); index++) {
        size_t kr = rdata.reactants[index];
        size_t n = kin.speciesPhaseIndex(kr);
        //klocal = kr - kin.start(n);
        size_t klocal = kr - kin.kineticsSpeciesIndex(0,n);
        kstoich = rdata.rstoich[index];
        const ThermoPhase& ph = kin.speciesPhase(kr);
        for (size_t m = 0; m < ph.nElements(); m++) {
            bal[ph.elementName(m)] -= kstoich*ph.nAtoms(klocal,m);
            balr[ph.elementName(m)] += kstoich*ph.nAtoms(klocal,m);
            //cout << "reactant species " << ph.speciesName(klocal) << " has " << ph.nAtoms(klocal,m)
            //     << " atoms of " << ph.elementName(m) << " and kstoich = " << kstoich << endl;
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

    string rptype;

    /*
     * The id of reactants and products are kept in child elements
     * of reaction, named "reactants" and "products". We search
     * the xml tree for these children based on the value of rp,
     * and store the xml element pointer here.
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
    vector<string> key, val;
    getPairs(rg, key, val);

    /*
     * Loop over each of the pairs and process them
     */
    doublereal ord, stch;
    string ph, sp;
    map<string, size_t> speciesMap;
    for (size_t n = 0; n < key.size(); n++) {
        sp = key[n]; // sp is the string name for species
        ph = "";
        /*
         * Search for the species in the kinetics object using the
         * member function kineticsSpeciesIndex(). We will search
         * for the species in all phases defined in the kinetics operator.
         */
        size_t isp = kin.kineticsSpeciesIndex(sp);
        if (isp == npos) {
            if (rules.skipUndeclaredSpecies) {
                return false;
            } else {
                throw CanteraError("getReagents",
                                   "Undeclared reactant or product species "+sp);
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
        //cout << key[n] << " " << isp << " " << stch << endl;

        /*
         * Needed to process reaction orders below.
         */
        speciesMap[sp] = order.size();
    }

    /*
     * Check to see if reaction orders have been specified.
     */
    if (rp == 1 && rxn.hasChild("order")) {
        vector<XML_Node*> ord;
        rxn.getChildren("order",ord);
        doublereal forder;
        for (size_t nn = 0; nn < ord.size(); nn++) {
            const XML_Node& oo = *ord[nn];
            string sp = oo["species"];
            size_t loc = speciesMap[sp];
            if (loc == 0)
                throw CanteraError("getReagents",
                                   "reaction order specified for non-reactant: "
                                   +sp);
            forder = fpValue(oo());
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
 * getArrhenius() parses the xml element called Arrhenius.
 * The Arrhenius expression is
 * \f[        k =  A T^(b) exp (-E_a / RT). \f]
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
        // to coverages used in the sticking probability
        // expression
        if (p.eosType() == cSurf || p.eosType() == cEdge) {
            sc = p.standardConcentration(klocal);
            f /= pow(sc, order);
        }
        // Otherwise:
        else {
            // We only allow one species to be in the phase
            // containing the special sticking coefficient
            // species.
            if (ispPhaseIndex == np) {
                not_surf++;
            }
            // Other bulk phase species on the other side
            // of ther interface are treated like surface
            // species.
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

static void getCoverageDependence(const XML_Node& node,
                                  thermo_t& surfphase, ReactionData& rdata)
{
    vector<XML_Node*> cov;
    node.getChildren("coverage", cov);
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
*/
static void getFalloff(const XML_Node& f, ReactionData& rdata)
{
    string type = f["type"];
    vector<string> p;
    getStringArray(f,p);
    vector_fp c;
    int np = static_cast<int>(p.size());
    for (int n = 0; n < np; n++) {
        c.push_back(fpValue(p[n]));
    }
    if (type == "Troe") {
        if (np == 4) {
            rdata.falloffType = TROE4_FALLOFF;
        } else if (np == 3) {
            rdata.falloffType = TROE3_FALLOFF;
        } else {
            throw CanteraError("getFalloff()", "Troe parameterization is specified by number of parameters, "
                               + int2str(np) + ", is not equal to 3 or 4");
        }
    } else if (type == "SRI") {
        if (np == 5) {
            rdata.falloffType = SRI5_FALLOFF;
            if (c[2] < 0.0) {
                throw CanteraError("getFalloff()", "SRI5 m_c parameter is less than zero: " + fp2str(c[2]));
            }
            if (c[3] < 0.0) {
                throw CanteraError("getFalloff()", "SRI5 m_d parameter is less than zero: " + fp2str(c[3]));
            }
        } else if (np == 3) {
            rdata.falloffType = SRI3_FALLOFF;
            if (c[2] < 0.0) {
                throw CanteraError("getFalloff()", "SRI3 m_c parameter is less than zero: " + fp2str(c[2]));
            }
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
    for (size_t n = 0; n < key.size(); n++) { // ; bb != ee; ++bb) {
        nm = key[n];// bb->first;
        size_t k = kin.kineticsSpeciesIndex(nm, phse);
        if (k != npos) {
            rdata.thirdBodyEfficiencies[k] = fpValue(val[n]); // bb->second;
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

    map<int, doublereal>::const_iterator b = r1.begin(), e = r1.end();
    int k1 = b->first;
    doublereal ratio = 0.0;
    if (r1[k1] == 0.0 || r2[k1] == 0.0) {
        goto next;
    }
    ratio = r2[k1]/r1[k1];
    ++b;
    for (; b != e; ++b) {
        k1 = b->first;
        if (r1[k1] == 0.0 || r2[k1] == 0.0) {
            goto next;
        }
        if (fabs(r2[k1]/r1[k1] - ratio) > 1.e-8) {
            goto next;
        }
    }
    return ratio;
next:
    ratio = 0.0;
    b = r1.begin();
    k1 = b->first;
    if (r1[k1] == 0.0 || r2[-k1] == 0.0) {
        return 0.0;
    }
    ratio = r2[-k1]/r1[k1];
    ++b;
    for (; b != e; ++b) {
        k1 = b->first;
        if (r1[k1] == 0.0 || r2[-k1] == 0.0) {
            return 0.0;
        }
        if (fabs(r2[-k1]/r1[k1] - ratio) > 1.e-8) {
            return 0.0;
        }
    }
    return ratio;
}

bool rxninfo::installReaction(int iRxn, const XML_Node& r, Kinetics& kin,
                              string default_phase, ReactionRules& rules,
                              bool validate_rxn)
{
    // Check to see that we are in fact at a reaction node
    if (r.name() != "reaction") {
        throw CanteraError(" rxninfo::installReaction",
                           " expected xml node reaction, got " + r.name());
    }

    // We use the ReactionData object to store initial values read in from the
    // xml data. Then, when we have collected everything we add the reaction to
    // the kinetics object, kin, at the end of the routine.
    ReactionData& rdata = **m_rdata.insert(m_rdata.end(), new ReactionData());
    rdata.validate = validate_rxn;

    // Check to see if the reaction is specified to be a duplicate of another
    // reaction. It's an error if the reaction is a duplicate and this is not
    // set.
    rdata.duplicate = (r.hasAttrib("duplicate")) ? 1 : 0;

    // Check to see if the reaction rate constant can be negative. It's an
    // error if a negative rate constant is found and this is not set.
    rules.allowNegativeA = (r.hasAttrib("negative_A")) ? 1 : 0;

    // Use the contents of the "equation" child element as the reaction's
    // string representation. Post-process to convert "[" and "]" characters
    // back into "<" and ">" which cannot easily be stored in an XML file. This
    // reaction string is used only for display purposes. It is not parsed for
    //  the identities of reactants or products.
    rdata.equation = (r.hasChild("equation")) ? r("equation") : "<no equation>";
    for (size_t nn = 0; nn < rdata.equation.size(); nn++) {
        if (rdata.equation[nn] == '[') {
            rdata.equation[nn] = '<';
        } else if (rdata.equation[nn] == ']') {
            rdata.equation[nn] = '>';
        }
    }

    // get the reactants
    bool ok = getReagents(r, kin, 1, default_phase, rdata.reactants,
                          rdata.rstoich, rdata.rorder, rules);

    // Get the products. We store the id of products in rdata.products
    ok = ok && getReagents(r, kin, -1, default_phase, rdata.products,
                           rdata.pstoich, rdata.porder, rules);

    // if there was a problem getting either the reactants or the products,
    // then abort.
    if (!ok) {
        return false;
    }

    // check whether the reaction is specified to be
    // reversible. Default is irreversible.
    string isrev = r["reversible"];
    rdata.reversible = (isrev == "yes" || isrev == "true");

    // If reaction orders are specified, then this reaction does not follow
    // mass-action kinetics, and is not an elementary reaction. So check that
    // it is not reversible, since computing the reverse rate from
    // thermochemistry only works for elementary reactions. Set the type to
    // global, so that kinetics managers will know to process the reaction
    // orders.
    if (r.hasChild("order")) {
        if (rdata.reversible == true)
            throw CanteraError("installReaction",
                               "reaction orders may only be given for "
                               "irreversible reactions");
        rdata.global = true;
    }

    // Some reactions can be elementary reactions but have fractional
    // stoichiometries wrt to some products and reactants. An example of these
    // are solid reactions involving phase transformations. Species with
    // fractional stoichiometries must be from single-species phases with
    // unity activities. For these reactions set the bool isReversibleWithFrac
    // to true.
    if (rdata.reversible == true) {
        for (size_t i = 0; i < rdata.products.size(); i++) {
            doublereal po = rdata.porder[i];
            AssertTrace(po == rdata.pstoich[i]);
            doublereal chk = po - 1.0 * int(po);
            if (chk != 0.0) {
                size_t k = rdata.products[i];
                // Special case when k is a single species phase.
                if (kin.speciesPhase(k).nSpecies() == 1) {
                    rdata.porder[i] = 0.0;
                }

                rdata.isReversibleWithFrac = true;
            }
        }
        for (size_t i = 0; i < rdata.reactants.size(); i++) {
            doublereal ro = rdata.rorder[i];
            AssertTrace(ro == rdata.rstoich[i]);
            doublereal chk = ro - 1.0 * int(ro);
            if (chk != 0.0) {
                size_t k = rdata.reactants[i];
                // Special case when k is a single species phase.
                if (kin.speciesPhase(k).nSpecies() == 1) {
                    rdata.rorder[i] = 0.0;
                }
                rdata.isReversibleWithFrac = true;
            }
        }
    }

    /*
     * Search the reaction element for the attribute "type".
     * If found, then branch on the type, to fill in appropriate
     * fields in rdata.
     */
    rdata.reactionType = ELEMENTARY_RXN;
    string typ = r["type"];
    if (typ == "falloff") {
        rdata.reactionType = FALLOFF_RXN;
        rdata.falloffType = SIMPLE_FALLOFF;
    } else if (typ == "chemAct") {
        rdata.reactionType = CHEMACT_RXN;
        rdata.falloffType = SIMPLE_FALLOFF;
    } else if (typ == "threeBody") {
        rdata.reactionType = THREE_BODY_RXN;
    } else if (typ == "plog") {
        rdata.reactionType = PLOG_RXN;
    } else if (typ == "chebyshev") {
        rdata.reactionType = CHEBYSHEV_RXN;
    } else if (typ == "surface") {
        rdata.reactionType = SURFACE_RXN;
    } else if (typ == "edge") {
        rdata.reactionType = EDGE_RXN;
    } else if (typ != "") {
        throw CanteraError("installReaction", "Unknown reaction type: " + typ);
    }

    rdata.number = iRxn;
    rdata.rxn_number = iRxn;

    // Read the rate coefficient data from the XML file. Trigger an
    // exception for negative A unless specifically authorized.
    getRateCoefficient(r.child("rateCoeff"), kin, rdata, rules);

    if (validate_rxn) {
        // Look for undeclared duplicate reactions.
        unsigned long int participants = 0;
        for (size_t nn = 0; nn < rdata.reactants.size(); nn++) {
            rdata.net_stoich[-1 - int(rdata.reactants[nn])] -= rdata.rstoich[nn];
            participants += rdata.reactants[nn];
        }
        for (size_t nn = 0; nn < rdata.products.size(); nn++) {
            rdata.net_stoich[int(rdata.products[nn])+1] += rdata.pstoich[nn];
            participants += 1000000 * rdata.products[nn];
        }

        vector<size_t>& related = m_participants[participants];
        for (size_t mm = 0; mm < related.size(); mm++) {
            ReactionData& other = *m_rdata[related[mm]];
            if (rdata.reactants.size() != other.reactants.size()) {
                continue; // different numbers of reactants
            } else if (rdata.reactionType != other.reactionType) {
                continue; // different reaction types
            } else if (rdata.duplicate && other.duplicate) {
                continue; // marked duplicates
            }
            doublereal c = isDuplicateReaction(rdata.net_stoich, other.net_stoich);
            if (c == 0) {
                continue; // stoichiometries differ (not by a multiple)
            } else if (c < 0.0 && !rdata.reversible && !other.reversible) {
                continue; // irreversible reactions in opposite directions
            } else if (rdata.reactionType == FALLOFF_RXN ||
                       rdata.reactionType == THREE_BODY_RXN ||
                       rdata.reactionType == CHEMACT_RXN) {
                bool thirdBodyOk = true;
                for (size_t k = 0; k < kin.nTotalSpecies(); k++) {
                    if (rdata.efficiency(k) * other.efficiency(k) != 0.0) {
                        thirdBodyOk = false;
                    }
                }
                if (thirdBodyOk) {
                    continue; // No overlap in third body efficiencies
                }
            }
            string msg = string("Undeclared duplicate reactions detected: \n")
                         +"Reaction "+int2str(other.number+1)+": "+other.equation
                         +"\nReaction "+int2str(iRxn+1)+": "+rdata.equation+"\n";
            throw CanteraError("installReaction", msg);
        }
        m_participants[participants].push_back(m_rdata.size() - 1);

        // Check to see that the elements balance in the reaction.
        // Throw an error if they don't
        checkRxnElementBalance(kin, rdata);
    }

    // Ok we have read everything in about the reaction. Add it  to the
    // kinetics object by calling the Kinetics member function addReaction()
    kin.addReaction(rdata);
    return true;
}

bool installReactionArrays(const XML_Node& p, Kinetics& kin,
                           std::string default_phase, bool check_for_duplicates)
{
    rxninfo _rxns;

    vector<XML_Node*> rarrays;
    int itot = 0;
    /*
     * Search the children of the phase element for the
     * xml element named reactionArray. If we can't find it,
     * then return signaling having not found any reactions.
     * Apparently, we allow multiple reactionArray elements here
     * Each one will be processed sequentially, with the
     * end result being purely additive.
     */
    p.getChildren("reactionArray",rarrays);
    int na = static_cast<int>(rarrays.size());
    if (na == 0) {
        kin.finalize();
        return false;
    }
    for (int n = 0; n < na; n++) {
        /*
         * Go get a reference to the current xml element,
         * reactionArray. We will process this element now.
         */
        const XML_Node& rxns = *rarrays[n];
        /*
         * The reactionArray element has an attribute called,
         * datasrc. The value of the attribute is the xml
         * element comprising the top of the
         * tree of reactions for the phase.
         * Find this datasrc element starting with the root
         * of the current xml node.
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
        ReactionRules rxnrule;
        if (rxns.hasChild("skip")) {
            const XML_Node& sk = rxns.child("skip");
            string sskip = sk["species"];
            if (sskip == "undeclared") {
                rxnrule.skipUndeclaredSpecies = true;
            }
            if (sk["third_bodies"] == "undeclared") {
                rxnrule.skipUndeclaredThirdBodies = true;
            }
        }
        int i, nrxns = 0;
        /*
         * Search for child elements called include. We only include
         * a reaction if it's tagged by one of the include fields.
         * Or, we include all reactions if there are no include fields.
         */
        vector<XML_Node*> incl;
        rxns.getChildren("include",incl);
        int ninc = static_cast<int>(incl.size());

        vector<XML_Node*> allrxns;
        rdata->getChildren("reaction",allrxns);
        nrxns = static_cast<int>(allrxns.size());
        // if no 'include' directive, then include all reactions
        if (ninc == 0) {
            for (i = 0; i < nrxns; i++) {
                const XML_Node* r = allrxns[i];
                if (r) {
                    if (_rxns.installReaction(itot, *r, kin,
                                              default_phase, rxnrule, check_for_duplicates)) {
                        ++itot;
                    }
                }
            }
        } else {
            for (int nii = 0; nii < ninc; nii++) {
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

                for (i = 0; i < nrxns; i++) {
                    const XML_Node* r = allrxns[i];
                    string rxid;
                    if (r) {
                        rxid = (*r)["id"];
                        if (iwild != string::npos) {
                            rxid = rxid.substr(0,iwild);
                        }
                        /*
                         * To decide whether the reaction is included or not
                         * we do a lexical min max and operation. This
                         * sometimes has surprising results.
                         */
                        if ((rxid >= imin) && (rxid <= imax)) {
                            if (_rxns.installReaction(itot, *r, kin,
                                                      default_phase, rxnrule, check_for_duplicates)) {
                                ++itot;
                            }
                        }
                    }
                }
            }
        }
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

    Kinetics& kin = *k;

    // This phase will be the owning phase for the kinetics operator
    // For interfaces, it is the surface phase between two volumes.
    // For homogeneous kinetics, it's the current volumetric phase.
    string owning_phase = phase["id"];

    bool check_for_duplicates = false;
    if (phase.parent()->hasChild("validate")) {
        const XML_Node& d = phase.parent()->child("validate");
        if (d["reactions"] == "yes") {
            check_for_duplicates = true;
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
                if (kin.phaseIndex(phase_id) == npos) {
                    kin.addPhase(*th[m]);
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
    kin.init();

    // Install the reactions.
    return installReactionArrays(phase, kin, owning_phase, check_for_duplicates);
}

bool buildSolutionFromXML(XML_Node& root, const std::string& id,
                          const std::string& nm, ThermoPhase* th, Kinetics* k)
{
    XML_Node* x;
    x = get_XML_NameID(nm, string("#")+id, &root);
    //            x = get_XML_Node(string("#")+id, &root);
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
    vector<ThermoPhase*> phases(1);
    phases[0] = th;
    /*
     * Fill in the kinetics object k, by querying the
     * const XML_Node tree located by x. The source terms and
     * eventually the source term vector will be constructed
     * from the list of ThermoPhases in the vector, phases.
     */
    importKinetics(*x, phases, k);
    return true;
}

}
