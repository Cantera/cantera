//------------------------------------------------
///
///  @file ReactionStoichMgr.cpp
///
///
//------------------------------------------------

#include "cantera/kinetics/ReactionStoichMgr.h"

#include "cantera/base/ctexceptions.h"
#include "cantera/kinetics/ReactionData.h"
#include "cantera/kinetics/Reaction.h"

#include <fstream>

using namespace std;

namespace Cantera
{
ReactionStoichMgr::ReactionStoichMgr()
{
    warn_deprecated("class ReactionStoichMgr",
        "To be removed after Cantera 2.2.");
    m_dummy.resize(10,1.0);
}

ReactionStoichMgr::ReactionStoichMgr(const  ReactionStoichMgr& right) :
    m_reactants(right.m_reactants),
    m_revproducts(right.m_revproducts),
    m_irrevproducts(right.m_irrevproducts),
    m_dummy(right.m_dummy)
{
}

ReactionStoichMgr& ReactionStoichMgr::operator=(const ReactionStoichMgr& right)
{
    if (this != &right) {

        m_reactants = right.m_reactants;
        m_revproducts = right.m_revproducts;
        m_irrevproducts = right.m_irrevproducts;
        m_dummy = right.m_dummy;
    }
    return *this;
}

void ReactionStoichMgr::add(size_t rxn, const std::vector<size_t>& reactants,
                            const std::vector<size_t>& products,
                            bool reversible)
{

    m_reactants.add(rxn, reactants);

    if (reversible) {
        m_revproducts.add(rxn, products);
    } else {
        m_irrevproducts.add(rxn, products);
    }
}

// Add the reaction into the stoichiometric manager
void ReactionStoichMgr::add(size_t rxn, const ReactionData& r)
{
    size_t k;
    std::vector<size_t> rk;
    doublereal frac;
    doublereal oo, os, of;
    bool doGlobal = false;
    std::vector<size_t> extReactants = r.reactants;
    vector_fp extRStoich = r.rstoich;
    vector_fp extROrder = r.rorder;

    //
    //  If we have a complete global reaction then we need to do something more complete
    //  than the previous treatment.  Basically we will use the reactant manager to calculate the
    //  global forward reaction rate of progress.
    //
    if (r.forwardFullOrder_.size() > 0) {
	//
	// Trigger a treatment where the order of the reaction and the stoichiometry
	// are treated as different.
	//
	doGlobal = true;
	size_t nsp = r.forwardFullOrder_.size();
	//
	// Set up a signal vector to indicate whether the species has been added into 
	// the input vectors for the stoich manager
	//
	vector_int kHandled(nsp, 0);
	//
	//  Loop over the reactants which are also nonzero stoichioemtric entries
	//  making sure the forwardFullOrder_ entries take precedence over rorder entries
	//  
	for (size_t kk = 0; kk < r.reactants.size(); kk++) {
	    k = r.reactants[kk];
	    os = r.rstoich[kk];
	    oo = r.rorder[kk];
	    of = r.forwardFullOrder_[k];
	    if (of != oo) {
		extROrder[kk] = of;
	    }
	    kHandled[k] = 1;
	}
	for (k = 0; k < nsp; k++) {
	    of = r.forwardFullOrder_[k];
	    if (of != 0.0) {
		if (kHandled[k] == 0) {
		    //
		    // Add extra entries to reactant inputs. Set their reactant stoichiometric entries to zero.
		    //
		    extReactants.push_back(k);
		    extROrder.push_back(of);
		    extRStoich.push_back(0.0);
		}
	    }
	}
    }
    
    bool isfrac = false;
    for (size_t n = 0; n < r.reactants.size(); n++) {
        size_t ns = size_t(r.rstoich[n]);
        frac = r.rstoich[n] - 1.0*int(r.rstoich[n]);
        if (frac != 0.0) {
            isfrac = true;
        }
        for (size_t m = 0; m < ns; m++) {
            rk.push_back(r.reactants[n]);
        }
    }

    //
    //  If the reaction is non-mass action add it in in a general way
    //  Reactants get extra terms for the forward reaction rate of progress
    //  that may have zero stoichiometries.
    //
    if (doGlobal) {
	m_reactants.add(rxn, extReactants, extROrder, extRStoich);
    } else {
	//
	//  this is confusing. The only issue should be whether rorder is different than rstoich!
	//
	if (isfrac || r.global || rk.size() > 3) {
	    m_reactants.add(rxn, r.reactants, r.rorder, r.rstoich);
	} else {
	    m_reactants.add(rxn, rk);
	}
    }

    std::vector<size_t> pk;
    isfrac = false;
    for (size_t n = 0; n < r.products.size(); n++) {
        size_t ns = size_t(r.pstoich[n]);
        frac = r.pstoich[n] - 1.0*int(r.pstoich[n]);
        if (frac != 0.0) {
            isfrac = true;
        }
        for (size_t m = 0; m < ns; m++) {
            pk.push_back(r.products[n]);
        }
    }

    if (r.reversible) {
	//
	//  this is confusing. The only issue should be whether porder is different than pstoich!
	//
        if (pk.size() > 3 || r.isReversibleWithFrac) {
            m_revproducts.add(rxn, r.products, r.porder, r.pstoich);
        } else {
            m_revproducts.add(rxn, pk);
        }
    } else {
	//
	//  this is confusing. The only issue should be whether porder is different than pstoich!
	//
	if (isfrac || pk.size() > 3) {
	    m_irrevproducts.add(rxn, r.products, r.porder, r.pstoich);
	} else {
	    m_irrevproducts.add(rxn, pk);
	}
    }
}

void ReactionStoichMgr::getCreationRates(size_t nsp, const doublereal* ropf,
                                         const doublereal* ropr, doublereal* c)
{
    // zero out the output array
    fill(c, c + nsp, 0.0);

    // the forward direction creates product species
    m_revproducts.incrementSpecies(ropf, c);
    m_irrevproducts.incrementSpecies(ropf, c);

    // the reverse direction creates reactant species
    m_reactants.incrementSpecies(ropr, c);
}

void ReactionStoichMgr::getDestructionRates(size_t nsp, const doublereal* ropf,
                                            const doublereal* ropr,
                                            doublereal* d)
{
    fill(d, d + nsp, 0.0);
    // the reverse direction destroys products in reversible reactions
    m_revproducts.incrementSpecies(ropr, d);
    // the forward direction destroys reactants
    m_reactants.incrementSpecies(ropf, d);
}

void ReactionStoichMgr::getNetProductionRates(size_t nsp, 
                                              const doublereal* ropnet,
                                              doublereal* w)
{
    fill(w, w + nsp, 0.0);
    // products are created for positive net rate of progress
    m_revproducts.incrementSpecies(ropnet, w);
    m_irrevproducts.incrementSpecies(ropnet, w);
    // reactants are destroyed for positive net rate of progress
    m_reactants.decrementSpecies(ropnet, w);
}

void ReactionStoichMgr::getReactionDelta(size_t nr, const doublereal* g,
                                         doublereal* dg)
{
    fill(dg, dg + nr, 0.0);
    // products add
    m_revproducts.incrementReactions(g, dg);
    m_irrevproducts.incrementReactions(g, dg);
    // reactants subtract
    m_reactants.decrementReactions(g, dg);
}

void ReactionStoichMgr::getRevReactionDelta(size_t nr, const doublereal* g,
                                            doublereal* dg)
{
    fill(dg, dg + nr, 0.0);
    m_revproducts.incrementReactions(g, dg);
    m_reactants.decrementReactions(g, dg);
}

void ReactionStoichMgr::multiplyReactants(const doublereal* c, doublereal* r)
{
    m_reactants.multiply(c, r);
}

void ReactionStoichMgr::multiplyRevProducts(const doublereal* c, doublereal* r)
{
    m_revproducts.multiply(c, r);
}

void ReactionStoichMgr::write(const string& filename)
{
    ofstream f(filename.c_str());
    f << "namespace mech {" << endl;
    writeCreationRates(f);
    writeDestructionRates(f);
    writeNetProductionRates(f);
    writeMultiplyReactants(f);
    writeMultiplyRevProducts(f);
    f << "} // namespace mech" << endl;
    f.close();
}

void ReactionStoichMgr::writeCreationRates(ostream& f)
{
    f << "    void getCreationRates(const doublereal* rf, const doublereal* rb," << endl;
    f << "          doublereal* c) {" << endl;
    map<size_t, string> out;
    m_revproducts.writeIncrementSpecies("rf",out);
    m_irrevproducts.writeIncrementSpecies("rf",out);
    m_reactants.writeIncrementSpecies("rb",out);
    map<size_t, string>::iterator b;
    for (b = out.begin(); b != out.end(); ++b) {
        string rhs = wrapString(b->second);
        rhs[1] = '=';
        f << "     c[" << b->first << "] " << rhs << ";" << endl;
    }
    f << "    }" << endl << endl << endl;
}

void ReactionStoichMgr::writeDestructionRates(ostream& f)
{
    f << "    void getDestructionRates(const doublereal* rf, const doublereal* rb," << endl;
    f << "          doublereal* d) {" << endl;
    map<size_t, string> out;
    m_revproducts.writeIncrementSpecies("rb",out);
    m_reactants.writeIncrementSpecies("rf",out);
    map<size_t, string>::iterator b;
    for (b = out.begin(); b != out.end(); ++b) {
        string rhs = wrapString(b->second);
        rhs[1] = '=';
        f << "     d[" << b->first << "] " << rhs << ";" << endl;
    }
    f << "    }" << endl << endl << endl;
}

void ReactionStoichMgr::writeNetProductionRates(ostream& f)
{
    f << "    void getNetProductionRates(const doublereal* r, doublereal* w) {" << endl;
    map<size_t, string> out;
    m_revproducts.writeIncrementSpecies("r",out);
    m_irrevproducts.writeIncrementSpecies("r",out);
    m_reactants.writeDecrementSpecies("r",out);
    map<size_t, string>::iterator b;
    for (b = out.begin(); b != out.end(); ++b) {
        string rhs = wrapString(b->second);
        rhs[1] = '=';
        f << "     w[" << b->first << "] " << rhs << ";" << endl;
    }
    f << "    }" << endl << endl << endl;
}

void ReactionStoichMgr::writeMultiplyReactants(ostream& f)
{
    f << "    void multiplyReactants(const doublereal* c, doublereal* r) {" << endl;
    map<size_t, string> out;
    m_reactants.writeMultiply("c",out);
    map<size_t, string>::iterator b;
    for (b = out.begin(); b != out.end(); ++b) {
        string rhs = b->second;
        f << "      r[" << b->first << "] *= " << rhs << ";" << endl;
    }
    f << "    }" << endl << endl << endl;
}

void ReactionStoichMgr::writeMultiplyRevProducts(ostream& f)
{
    f << "    void multiplyRevProducts(const doublereal* c, doublereal* r) {" << endl;
    map<size_t, string> out;
    m_revproducts.writeMultiply("c",out);
    map<size_t, string>::iterator b;
    for (b = out.begin(); b != out.end(); ++b) {
        string rhs = b->second;
        f << "      r[" << b->first << "] *= " << rhs << ";" << endl;
    }
    f << "    }" << endl << endl << endl;
}
}
