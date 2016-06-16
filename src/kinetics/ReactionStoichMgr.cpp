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
#include "cantera/kinetics/RxnRates.h"
#include "cantera/thermo/SpeciesThermo.h"

#include <fstream>
#include <iomanip>

using namespace std;

namespace Cantera
{
ReactionStoichMgr::ReactionStoichMgr()
{
  m_ii     = 0;
  m_kk     = 0;
  m_nrev   = 0;
  m_nirrev = 0;
  m_dummy.resize(10,1.0);
}

ReactionStoichMgr::~ReactionStoichMgr()
{
}

ReactionStoichMgr::ReactionStoichMgr(const  ReactionStoichMgr& right) :
    m_kk(right.m_kk),
    m_ii(right.m_ii),
    m_nrev(right.m_nrev),
    m_nirrev(right.m_nirrev),
    m_reactants(right.m_reactants),
    m_revproducts(right.m_revproducts),
    m_irrevproducts(right.m_irrevproducts),
    m_dummy(right.m_dummy)
{
}

ReactionStoichMgr& ReactionStoichMgr::operator=(const ReactionStoichMgr& right)
{
    if (this != &right) {
        m_kk = right.m_kk;
	m_ii = right.m_ii;
	m_nrev = right.m_nrev;
	m_nirrev = right.m_nirrev;
        m_reactants = right.m_reactants;
        m_revproducts = right.m_revproducts;
        m_irrevproducts = right.m_irrevproducts;
	m_dn    = right.m_dn;
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

void ReactionStoichMgr::add(size_t rxn, const ReactionData& r)
{

    doublereal nsFlt;
    doublereal reactantGlobalOrder = 0.0;
    doublereal productGlobalOrder  = 0.0;

    std::vector<size_t> rk;
    doublereal frac;
    bool isfrac = false;
    for (size_t n = 0; n < r.reactants.size(); n++) {
        size_t ns = size_t(r.rstoich[n]);
	nsFlt = r.rstoich[n];
	reactantGlobalOrder += nsFlt;
        frac  = r.rstoich[n] - 1.0*int(r.rstoich[n]);
        if (frac != 0.0) {
            isfrac = true;
        }
        for (size_t m = 0; m < ns; m++) {
            rk.push_back(r.reactants[n]);
        }
    }
    
    // install rate coefficient calculators
    if(r.reactionType == TEDEP_RXN) {
      
      m_tedep_rates.install(rxn, r);
      m_deltaE.push_back(r.deltaE);
      m_tedep_index.push_back(rxn);
      
    } else if(r.reactionType == VIBREL_RXN) {
      
      m_vibrel_rates.install(rxn, r);
      m_vibrel_index.push_back(rxn);
      
    } else {
      
      m_rates.install(rxn, r);
      m_crxn_index.push_back(rxn);
      
    }

    if(r.reactionType == THREE_BODY_RXN){
      m_3b_concm.install(rxn, r.thirdBodyEfficiencies, r.default_3b_eff);
      m_3b_index.push_back(rxn);
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
	nsFlt = r.pstoich[n];
	productGlobalOrder += nsFlt;
        frac  = r.pstoich[n] - 1.0*int(r.pstoich[n]);
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

    if (r.reversible) {
      m_revindex.push_back(rxn);
      ++m_nrev;
    } else {
      m_irrevindex.push_back(rxn);
      ++m_nirrev;
    }

    m_dn.push_back(productGlobalOrder - reactantGlobalOrder);
          
    ++m_ii;
    
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

void ReactionStoichMgr::
writeMech(size_t ns, const string& filename, SpeciesThermo& sp)
{

  m_kk = ns;
  SpeciesThermo& m_sp = sp;
  ofstream f(filename.c_str());
  
  f << "namespace mech {" << endl;
  f << endl;
  f << "  template <class Type>" << endl;
  writeSpeciesThermo(m_sp, f);
  f << endl;
  f << "  template <class Type>" << endl;
  writeRateConstants(f);
  f << endl;
  f << "  template <class Type>" << endl;
  writeSpeciesSourceTerms(f);
  f << endl;
  f << "} // namespace mech" << endl;

  f.close();
  
}

void ReactionStoichMgr::
writeRateConstants(ostream& f)
{

  f << "  void updateRateConstants(int& ns, int& nr, double& p, Type& T, vector<Type>& C," << endl;
  f << "			   vector<Type>& kfwd, vector<Type>& krev) {" << endl;
  f << endl;
  f << "    Type         tlog = log(T);"  << endl;
  f << "    Type         rt   = 1.0 / T;" << endl;
  f << "    vector<Type> keqs(nr, 0.0);"     << endl;
  f << endl;
  f << "    getEqConstants(ns, p, T, keqs);" << endl;
  f << endl;
  
  size_t crxnindx   = 0;
  size_t tedepindx  = 0;
  size_t vibrelindx = 0;
  for(int i = 0; i < m_ii; ++i) {
    f << "    kfwd[" << i << "] = ";
    if(m_tedep_index[tedepindx] == i) {
      m_tedep_rates.m_rates[tedepindx].writeUpdateRHS(f);
      ++tedepindx;
    } else if(m_vibrel_index[vibrelindx] == i) {
      m_vibrel_rates.m_rates[vibrelindx].writeUpdateRHS(f);
      ++vibrelindx;
    } else if(m_crxn_index[crxnindx] == i) {
      m_rates.m_rates[crxnindx].writeUpdateRHS(f);
      ++crxnindx;
    }
  }
  f << endl;

  // check for third-body reactions
  if(m_3b_concm.m_n != 0) {
    size_t         j = 0;
    size_t         m = 0;
    vector<size_t> index;
    m_3b_eff.resize(m_kk);
    
    for(size_t r = 0; r < m_3b_concm.m_n; r++) {
      j     = 1;
      m     = m_3b_index[r];
      index = m_3b_index;
      m_3b_concm.m_concm[r].getEfficiencies(m_3b_eff);
      
      f << "    kfwd[" << m << "] = (";
      if(index[0] == 0){
      	f << scientific << setprecision(1) << m_3b_eff[0] << " * N[0]";
      }
      else{
      	f << scientific << setprecision(1) << m_3b_concm.m_concm[r].m_deflt << " * N[0]";
      }
      for(size_t i = 1; i < m_kk; i++) {
      	if( i == index[j]) {
      	  f << " + " << scientific << setprecision(1) << m_3b_eff[i]
      	    << " * N[" << i << "]";
      	  j = j + 1;
      	}
      	else{
      	  f << " + " << scientific << setprecision(1) << m_3b_concm.m_concm[r].m_deflt
      	    << " * N[" << i << "]"; 
      	}
      }
      f << ") * kfwd[" << m << "];" << endl;
      
      }
  }
  f << endl;

  // reverse rate coefficients
  f << "    for(int i = 0; i < nr; ++i) { krev[i] = kfwd[i] * keqs[i]; }" << endl;
  f << endl;
  
  // end subroutine
  f << "  }" << endl;
  
}

void ReactionStoichMgr::
writeSpeciesSourceTerms(ostream& f)
{

  f << "  void updateSpeciesSourceTerm(int& ns, int& nr, double* mw, double& p, double& h,"
    << endl;
  f << "			       double& Told, vector<Type>& Y,"
    << "vector<Type>& omega) {" << endl;
  f << endl;
  f << endl;
  f << "    Type           rho;"           << endl;
  f << "    Type           W;"             << endl;
  f << "    Type           T;"             << endl;
  f << "    vector<Type>   C(ns,    0.0);" << endl;
  f << "    vector<Type>   kfwd(nr, 0.0);" << endl;
  f << "    vector<Type>   krev(nr, 0.0);" << endl;
  f << "    vector<Type>   Rfwd(nr, 0.0);" << endl;
  f << "    vector<Type>   Rrev(nr, 0.0);" << endl;
  f << "    vector<Type>   Rnet(nr, 0.0);" << endl;
  f << endl;
  f << "    mech::getTemperature(ns, mw, h, Told, Y, T);" << endl;
  f << endl;
  f << "    W   = 0.0;" << endl;
  f << "    for(int i = 0; i < ns; ++i) { W += Y[i] / mw[i]; }" << endl;
  f << "    W   = 1.0/W;" << endl;
  f << "    rho = (p * W)/(GasConstant * T);" << endl;
  f << "    for(int i = 0; i < ns; ++i) { C[i] = rho * Y[i] / mw[i]; }" << endl;
  f << "    mech::updateRateConstants(ns, nr, p, T, C, kfwd, krev);" << endl;
  f << endl;

  map<size_t, string> outfwd;
  m_reactants.writeMultiply("N",outfwd);
  map<size_t, string>::iterator bfwd;
  for (bfwd = outfwd.begin(); bfwd != outfwd.end(); ++bfwd) {
    string rhsfwd = bfwd->second;
    f << "    Rfwd[" << bfwd->first << "] = kfwd[" << bfwd->first
	<< "] * " << rhsfwd << ";" << endl;
  }
  f << endl;
  map<size_t, string> outrev;
  m_revproducts.writeMultiply("N",outrev);
  map<size_t, string>::iterator brev;
  for (brev = outrev.begin(); brev != outrev.end(); ++brev) {
    string rhsrev = brev->second;
    f << "    Rrev[" << brev->first << "] = krev[" << brev->first
      << "] * " << rhsrev << ";" << endl;
  }
  f << endl;
  f << "    for(int i = 0; i < nr; ++i) { Rnet[i] = Rfwd[i] - Rrev[i]; }" << endl;
  f << endl;
  map<size_t, string> out3;
  m_revproducts.writeIncrementSpecies("Rnet",out3);
  m_irrevproducts.writeIncrementSpecies("Rnet",out3);
  m_reactants.writeDecrementSpecies("Rnet",out3);
  map<size_t, string>::iterator b3;
  for (b3 = out3.begin(); b3 != out3.end(); ++b3) {
    string rhs = wrapString(b3->second);
    f << "    omega[" << b3->first << "] = " << rhs << ";" << endl;
  }

  f << endl;
  f << "  }" << endl;
   
}

void ReactionStoichMgr::
writeSpeciesThermo(SpeciesThermo& sp, ostream& f)
{

  f.setf(ios::scientific);
  
  int        type;
  doublereal midTemp;
  doublereal minTemp;
  doublereal maxTemp;
  doublereal refPressure;
  doublereal c[15];

  f << "  void getSpecificHeats_R(Type& T, vector<Type>& cp0_R) {" << endl;
  f << endl;
  f << "    Type tt0 = T;" << endl;
  f << "    Type tt1 = T * tt0;" << endl;
  f << "    Type tt2 = T * tt1;" << endl;
  f << "    Type tt3 = T * tt2;" << endl;
  f << endl;
  for (size_t i = 0; i < m_kk; i++) {
    sp.reportParams(i, type, c, minTemp, maxTemp, refPressure);
    midTemp = c[0];
    f << "    if(tt0 < " << setprecision(4) << midTemp << ") {" << endl;
    f << "      cp0_R[" << i << "] = " << setprecision(15)
      << c[1] << " + "
      << c[2] << " * tt0 + "
      << c[3] << " * tt1 + "
      << c[4] << " * tt2 + "
      << c[5] << " * tt3;"
      << endl;
    f << "    } else {" << endl;
    f << "      cp0_R[" << i << "] = "
      << c[8]  << " + "
      << c[9]  << " * tt0 + "
      << c[10] << " * tt1 + "
      << c[11] << " * tt2 + "
      << c[12] << " * tt3;"
      << endl;
    f << "    };" << endl;
    f << endl;
  };
  f << endl;
  f << "  };" << endl;
  f << endl;

  f << "  template <class Type>" << endl;
  f << "  void getEnthalpies_RT(Type& T, vector<Type>& h0_RT) {" << endl;
  f << endl;
  f << "    Type tt0 = T;" << endl;
  f << "    Type tt1 = T * tt0;" << endl;
  f << "    Type tt2 = T * tt1;" << endl;
  f << "    Type tt3 = T * tt2;" << endl;
  f << "    Type tt4 = 1.0 / T;" << endl;
  f << endl;
  for (size_t i = 0; i < m_kk; i++) {
    sp.reportParams(i, type, c, minTemp, maxTemp, refPressure);
    midTemp = c[0];
    f << "    if(tt0 < " << setprecision(4) << midTemp << ") {" << endl;
    f << "      h0_RT[" << i << "] = " << setprecision(15)
      << c[1] << " + "
      << c[2] << " * tt0 * 0.50 + "
      << c[3] << " * tt1 * OneThird + "
      << c[4] << " * tt2 * 0.25 + "
      << c[5] << " * tt3 * 0.20 + "
      << c[6] << " * tt4;"
      << endl;
    f << "    } else {" << endl;
    f << "      h0_RT[" << i << "] = "
      << c[8]  << " + "
      << c[9]  << " * tt0 * 0.50 + "
      << c[10] << " * tt1 * OneThird + "
      << c[11] << " * tt2 * 0.25 + "
      << c[12] << " * tt3 * 0.20 + "
      << c[13] << " * tt4;"
      << endl;
    f << "    };" << endl;
  f << endl;
  };
  f << "  };" << endl;
  f << endl;

  f << "  template <class Type>" << endl;
  f << "  void getEnthalpiesDerivatives(Type& T, vector<Type>& dh0dT) {" << endl;
  f << endl;
  f << "    Type tt0 = T;" << endl;
  f << "    Type tt1 = T * tt0;" << endl;
  f << "    Type tt2 = T * tt1;" << endl;
  f << "    Type tt3 = T * tt2;" << endl;
  f << endl;
  for (size_t i = 0; i < m_kk; i++) {
    sp.reportParams(i, type, c, minTemp, maxTemp, refPressure);
    midTemp = c[0];
    f << "    if(tt0 < " << setprecision(4) << midTemp << ") {" << endl;
    f << "      dh0dT[" << i << "] = " << setprecision(15)
      << c[1] << " + "
      << c[2] << " * tt0 + "
      << c[3] << " * tt1 + "
      << c[4] << " * tt2 + "
      << c[5] << " * tt3;"
      << endl;
    f << "    } else {" << endl;
    f << "      dh0dT[" << i << "] = "
      << c[8]  << " + "
      << c[9]  << " * tt0 + "
      << c[10] << " * tt1 + "
      << c[11] << " * tt2 + "
      << c[12] << " * tt3;"
      << endl;
    f << "    };" << endl;
  f << endl;
  };
  f << "  };" << endl;
  f << endl;

  f << "  template <class Type>" << endl;
  f << "  void getEntropies_R(Type& T, vector<Type>& s0_R) {" << endl;
  f << endl;
  f << "    Type tt0 = T;"        << endl;
  f << "    Type tt1 = T * tt0;" << endl;
  f << "    Type tt2 = T * tt1;"  << endl;
  f << "    Type tt3 = T * tt2;"  << endl;
  f << "    Type tt4 = 1.0 / T;"  << endl;
  f << "    Type tt5 = log(T);"   << endl;
  f << endl;
  for (size_t i = 0; i < m_kk; i++) {
    sp.reportParams(i, type, c, minTemp, maxTemp, refPressure);
    midTemp = c[0];
    f << "    if(tt0 < " << setprecision(4) << midTemp << ") {" << endl;
    f << "      s0_R[" << i << "] = " << setprecision(15)
      << c[1] << " * tt5 + "
      << c[2] << " * tt0 + "
      << c[3] << " * tt1 * 0.50 + "
      << c[4] << " * tt2 * OneThird + "
      << c[5] << " * tt3 * 0.25 + "
      << c[7] << ";"
      << endl;
    f << "    } else {" << endl;
    f << "      s0_R[" << i << "] = "
      << c[8]  << " * tt5 +  "
      << c[9]  << " * tt0 + "
      << c[10] << " * tt1 * 0.50 + "
      << c[11] << " * tt2 * OneThird + "
      << c[12] << " * tt3 * 0.25 + "
      << c[14] << ";"
      << endl;
    f << "    };" << endl;
  f << endl;
  };
  f << "  };" << endl;
  f << endl;

  f << "  template <class Type>" << endl;
  f << "  void getGibbsFunctions_RT(int& ns, double& p, Type& T,"
    <<   "vector<Type>& g0_RT) {"
    << endl;
  f << endl;
  f << "    vector<Type> h0_RT(ns, 0.0);" << endl;
  f << "    vector<Type> s0_R(ns,  0.0);" << endl;
  f << endl;
  f << "    getEnthalpies_RT(T, h0_RT);" << endl;
  f << "    getEntropies_R(T, s0_R);"    << endl;
  f << "    for(int i = 0; i < ns; ++i) { g0_RT[i] = h0_RT[i] - s0_R[i]; }" << endl;
  f << endl;
  f << "  };" << endl;
  f << endl;

  f << "  template <class Type>" << endl;
  f << "  void getTemperature(int& ns, double* mw, double& h, double& Told," << endl;
  f <<   "                    vector<Type>& Y, Type& T) {"
    << endl;
  f << endl;
  f << "    double       tol   = 1.0e-06;"  << endl;
  f << "    int          niter = 500;"      << endl;
  f << "    Type         m_RT;"             << endl;
  f << "    Type         m_To;"             << endl;
  f << "    Type         m_Tp;"             << endl;
  f << "    Type         m_dT = 1.0;"       << endl;
  f << "    Type         m_FY = 1.0;"       << endl;
  f << "    Type         m_JY = 1.0;"       << endl;
  f << "    vector<Type> m_hi(ns,    0.0);" << endl;
  f << "    vector<Type> m_dhidT(ns, 0.0);" << endl;
  f << endl;
  f << "    m_To = Told;" << endl;
  f << "    m_Tp = Told;" << endl;
  f << endl;
  f << "    for(int i = 0; i < ns; ++i) { " << endl;
  f << endl;
  f << "      m_RT = GasConstant * m_To;" << endl;
  f << "      mech::getEnthalpies_RT(m_To, m_hi);" << endl;
  f << "      mech::getEnthalpiesDerivatives(m_To, m_dhidT);" << endl;
  f << "      for(int i = 0; i < ns; ++i) { "
    << "m_hi[i]    = m_RT * m_hi[i] / mw[i]; }" << endl;
  f << "      for(int i = 0; i < ns; ++i) { "
    << "m_dhidT[i] = GasConstant * m_dhidT[i] / mw[i]; }" << endl;
  f << "      for(int i = 0; i < ns; ++i) { m_FY -= m_hi[i] * Y[i]; }" << endl;
  f << "      for(int i = 0; i < ns; ++i) { m_JY -= m_dhidT[i] * Y[i]; }" << endl;
  f << "      m_FY = h + m_FY;" << endl;
  f << "      m_JY = 1.0 / m_JY;" << endl;
  f << "      m_dT = -m_FY * m_JY;" << endl;
  f << "      m_Tp = m_To + m_dT;" << endl;
  f << "      m_To = m_Tp;" << endl;
  f << endl;
  f << "      if( (fabs(m_dT) < tol)) {" << endl;
  f << "	T = m_Tp;" << endl;
  f << "	return;" << endl;
  f << "      }" << endl;
  f << endl;
  f << "      m_FY = 0.0;" << endl;
  f << "      m_JY = 0.0;" << endl;
  f << endl;
  f << "    }" << endl;
  f << endl;
  f << "    T = m_Tp;" << endl;
  f << endl;
  f << "  };" << endl;
  f << endl;

  f << "  template <class Type>" << endl;
  f << "  void getEqConstants(int& ns, double& p, Type& T,"
    <<   "vector<Type>& keqs) {" << endl;
  f << endl;
  f << "    double       p0 = OneAtm;"    << endl;
  f << "    Type         RT = GasConstant * T;" << endl;
  f << "    Type         C0 = p0 / RT;"   << endl;            
  f << "    vector<Type> g0_RT(ns, 0.0);" << endl;
  f << endl;
  f << "    getGibbsFunctions_RT(ns, p, T, g0_RT);" << endl;
  f << "    for(int i = 0; i < ns; ++i) { g0_RT[i] = exp(g0_RT[i]); }"
    << endl;
  f << endl;
  
  map<size_t, string> outfwd;
  map<size_t, string> outrev;
  m_reactants.writeMultiply("g0_RT",outfwd);
  m_revproducts.writeMultiply("g0_RT",outrev);
  map<size_t, string>::iterator bfwd;
  map<size_t, string>::iterator brev;
  for(brev = outrev.begin(); brev != outrev.end(); ++brev) {
    string rhsrev = brev->second;
    if (m_dn[brev->first] < 0) {
      f << "    keqs[" << brev->first << "] = (C0 * "
  	<< rhsrev << ");" << endl;
    } else {
      f << "    keqs[" << brev->first << "] = (" << rhsrev
  	<< ");" << endl;
    }
  }
  f << endl;
  int j = 0;
  for(bfwd = outfwd.begin(); bfwd != outfwd.end(); ++bfwd) {
    size_t irxn = m_revindex[j];
    if(bfwd->first == irxn) {
      string rhsfwd = bfwd->second;
      if (m_dn[bfwd->first] > 0) {
  	f << "    keqs[" << bfwd->first << "] /= (C0 * "
  	  << rhsfwd << ");" << endl;
      } else {
  	f << "    keqs[" << bfwd->first << "] /= (" << rhsfwd
  	  << ");" << endl;
      }
      ++j;
    }
  }
  f << endl;;
  for (int i = 0; i < m_nirrev; ++i) {
    size_t irxn = m_irrevindex[i];
    f << "    keqs[" << irxn << "] = 0.0;" << endl;
  }
  f << endl;
  
  f << "  };" << endl;
  f << endl;
  
}
  
}
