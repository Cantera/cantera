//------------------------------------------------
///
///  @file ReactionStoichMgr.cpp
///
///
//------------------------------------------------

#include "cantera/kinetics/ReactionStoichMgr.h"

#include "cantera/kinetics/StoichManager.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/kinetics/ReactionData.h"

#include <fstream>

using namespace std;

namespace Cantera
{
ReactionStoichMgr::ReactionStoichMgr()
{
    m_dummy.resize(10,1.0);
}

ReactionStoichMgr::~ReactionStoichMgr()
{
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

void ReactionStoichMgr::
add(size_t rxn, const std::vector<size_t>& reactants,
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

void ReactionStoichMgr::
add(size_t rxn, const ReactionData& r)
{

    std::vector<size_t> rk;
    doublereal frac;
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

    // if the reaction has fractional stoichiometric coefficients
    // or specified reaction orders, then add it in a general reaction
    if (isfrac || r.global || rk.size() > 3) {
        m_reactants.add(rxn, r.reactants, r.rorder, r.rstoich);
    } else {
        m_reactants.add(rxn, rk);
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
        if (isfrac && !r.isReversibleWithFrac) {
            throw CanteraError("ReactionStoichMgr::add",
                               "Fractional product stoichiometric coefficients only allowed "
                               "\nfor irreversible reactions and most reversible reactions");
        }
        if (pk.size() > 3 || r.isReversibleWithFrac) {
            m_revproducts.add(rxn, r.products, r.porder, r.pstoich);
        } else {
            m_revproducts.add(rxn, pk);
        }
    } else if (isfrac || pk.size() > 3) {
        m_irrevproducts.add(rxn, r.products, r.porder,  r.pstoich);
    } else {
        m_irrevproducts.add(rxn, pk);
    }
}

void ReactionStoichMgr::
getCreationRates(size_t nsp, const doublereal* ropf,
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

void ReactionStoichMgr::
getDestructionRates(size_t nsp, const doublereal* ropf,
                    const doublereal* ropr, doublereal* d)
{
    fill(d, d + nsp, 0.0);
    // the reverse direction destroys products in reversible reactions
    m_revproducts.incrementSpecies(ropr, d);
    // the forward direction destroys reactants
    m_reactants.incrementSpecies(ropf, d);
}

void ReactionStoichMgr::
getNetProductionRates(size_t nsp, const doublereal* ropnet, doublereal* w)
{
    fill(w, w + nsp, 0.0);
    // products are created for positive net rate of progress
    m_revproducts.incrementSpecies(ropnet, w);
    m_irrevproducts.incrementSpecies(ropnet, w);
    // reactants are destroyed for positive net rate of progress
    m_reactants.decrementSpecies(ropnet, w);
}

void ReactionStoichMgr::
getReactionDelta(size_t nr, const doublereal* g, doublereal* dg)
{
    fill(dg, dg + nr, 0.0);
    // products add
    m_revproducts.incrementReactions(g, dg);
    m_irrevproducts.incrementReactions(g, dg);
    // reactants subtract
    m_reactants.decrementReactions(g, dg);
}

void ReactionStoichMgr::
getRevReactionDelta(size_t nr, const doublereal* g, doublereal* dg)
{
    fill(dg, dg + nr, 0.0);
    m_revproducts.incrementReactions(g, dg);
    m_reactants.decrementReactions(g, dg);
}

void ReactionStoichMgr::
multiplyReactants(const doublereal* c, doublereal* r)
{
    m_reactants.multiply(c, r);
}

void ReactionStoichMgr::
multiplyRevProducts(const doublereal* c, doublereal* r)
{
    m_revproducts.multiply(c, r);
}

void ReactionStoichMgr::
write(const string& filename)
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

void ReactionStoichMgr::
writeCreationRates(ostream& f)
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

void ReactionStoichMgr::
writeDestructionRates(ostream& f)
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

void ReactionStoichMgr::
writeNetProductionRates(ostream& f)
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

void ReactionStoichMgr::
writeMultiplyReactants(ostream& f)
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

void ReactionStoichMgr::
writeMultiplyRevProducts(ostream& f)
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
