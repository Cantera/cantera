/**
 *  @file StoichManager.h
 *
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_STOICH_MGR_H
#define CT_STOICH_MGR_H

#include <vector>
#include <map>
using namespace std;

#include "stringUtils.h"

namespace Cantera {

  /** 
   * @defgroup Stoichiometry Stoichiometry
   *
   * Note: these classes are designed for internal use in class
   * ReactionStoichManager.
   * 
   * Operations on reactions that require knowing the reaction
   * stoichiometry.  This module consists of class StoichManager, and
   * classes C1, C2, and C3.  Classes C1, C2, and C3 handle operations
   * involving one, two, or three species, respectively, in a
   * reaction. Instances are instantiated with a reaction number, and n
   * species numbers (n = 1 for C1, etc.).  All three classes have the
   * same interface.
   * 
   * These classes are designed for use by StoichManager, and the
   * operations implemented are those needed to efficiently compute
   * quantities such as rates of progress, species production rates,
   * reaction thermochemistry, etc. The compiler will inline these
   * methods into the body of the corresponding StoichManager method,
   * and so there is no performance penalty (unless inlining is turned
   * off).
   * 
   * To describe the methods, consider class C3 and suppose an instance
   * is created with reaction number irxn and species numbers k0, k1,
   * and k2.
   *
   *  - multiply(in, out) : out[irxn] is multiplied by 
   *    in[k0] * in[k1] * in[k2]
   *
   *  - power(in, out) : out[irxn] is multiplied by 
   *     (in[k0]^order0) * (in[k1]^order1) * (in[k2]^order2)
   * 
   *  - incrementReaction(in, out) : out[irxn] is incremented by 
   *    in[k0] + in[k1] + in[k2]
   *
   *  - decrementReaction(in, out) : out[irxn] is decremented by 
   *    in[k0] + in[k1] + in[k2]
   *
   *  - incrementSpecies(in, out)  : out[k0], out[k1], and out[k2]
   *    are all incremented by in[irxn]
   *
   *  - decrementSpecies(in, out)  : out[k0], out[k1], and out[k2]
   *    are all decremented by in[irxn]
   *
   * The function multiply() is usually used when evaluating the
   * forward and reverse rates of progress of reactions.
   * The rate constants are usually loaded into out[]. Then 
   * multply() is called to add in the dependence of the 
   * species concentrations to yield a forward and reverse rop.
   *
   * The function incrementSpecies() and its cousin decrementSpecies()
   * is used to translate from rates of progress to species production 
   * rates. The vector in[] is preloaed with the rates of progess of 
   * all reactions. Then incrementSpecies() is called to 
   * increment the species production vector, out[], with the rates
   * of progress. 
   *
   * The functions incrementReaction() and decrementReaction() are
   * used to find the standard state equilibrium constant for
   * a reaction. Here, output[] is a vector of length
   * number of reactions, usually the standard gibbs free energies
   * of reaction, while input, usually the standard state
   * gibbs free energies of species, is a vector of length number of
   * species.
   * 
   * Note the stoichiometric coefficient for a species in a reaction
   * is handled by always assuming it is equal to one and then 
   * treating reactants and products for a reaction separately.
   * Bimolecular reactions involving the identical species are 
   * treated as involving separate species.
   *
   * @internal This class should be upgraded to include cases where
   * real stoichiometric coefficients are used. Shouldn't be that
   * hard to do, and they occur in engineering simulations with some
   * regularity.
   * 
   */


  /**
   * Handles one species in a reaction.
   * @ingroup Stoichiometry
   * @internal
   */
  class C1 {

  public:

    C1( int rxn = 0, int ic0 = 0, doublereal order = 1.0 ) 
      : m_rxn (rxn),  m_ic0 (ic0), m_order(order) {}

    int data(vector<int>& ic) {
      ic.resize(3);
      ic[0] = m_ic0;
      return m_rxn;
    }

    void multiply(const doublereal* input, doublereal* output) const { 
      *(output + m_rxn) *= *(input + m_ic0); 
    }
    void power(const doublereal* input, doublereal* output) const { 
      output[m_rxn] *= pow(input[m_ic0], m_order); 
    }
    void incrementSpecies(const doublereal* input, 
			  doublereal* output) const { 
      *(output + m_ic0) += *(input + m_rxn); 
    }
    void decrementSpecies(const doublereal* input, 
			  doublereal* output) const { 
      *(output + m_ic0) -= *(input + m_rxn); 
    }
    void incrementReaction(const doublereal* input, 
			   doublereal* output) const { 
      *(output + m_rxn) += *(input + m_ic0); 
    }
    void decrementReaction(const doublereal* input, 
			   doublereal* output) const { 
      *(output + m_rxn) -= *(input + m_ic0); 
    }
  private:
    int m_rxn, m_ic0;
    doublereal m_order;
  };
    


  /**
   * Handles two species in a single reaction.
   * @ingroup Stoichiometry
   */
  class C2 {
  public:
    C2( int rxn = 0, int ic0 = 0, int ic1 = 0,
	doublereal order0 = 1.0, doublereal order1 = 1.0 ) 
      : m_rxn (rxn), m_ic0 (ic0), m_ic1 (ic1), 
      m_order0(order0), m_order1(order1) {}

    int data(vector<int>& ic) {
      ic.resize(2);
      ic[0] = m_ic0;
      ic[1] = m_ic1;
      return m_rxn;
    }

    void multiply(const doublereal* input, doublereal* output) const {
      output[m_rxn] *= input[m_ic0] * input[m_ic1]; 
    }

    void power(const doublereal* input, doublereal* output) const { 
      output[m_rxn] *= pow(input[m_ic0],m_order0) * 
	pow(input[m_ic1],m_order1); 
    }

    void incrementSpecies(const doublereal* input, 
			  doublereal* output) const {
      doublereal x = input[m_rxn]; 
      output[m_ic0] += x;
      output[m_ic1] += x;
    }

    void decrementSpecies(const doublereal* input, 
			  doublereal* output) const {
      doublereal x = input[m_rxn]; 
      output[m_ic0] -= x;
      output[m_ic1] -= x; 
    }

    void incrementReaction(const doublereal* input, 
			   doublereal* output) const { 
      *(output + m_rxn) += *(input + m_ic0) + *(input + m_ic1); 
    }

    void decrementReaction(const doublereal* input, 
			   doublereal* output) const { 
      *(output + m_rxn) -= (*(input + m_ic0) + *(input + m_ic1)); 
    }

  private:

    /**
     * Reaction index -> index into the ROP vector
     */
    int m_rxn;

    /**
     * Species indecise -> index into the species vector for the
     * two species.
     */
    int m_ic0, m_ic1;
    doublereal m_order0, m_order1;
  };
  

  /**
   * Handles three species in a reaction.
   * @ingroup Stoichiometry
   */  
  class C3 {
  public:
    C3( int rxn = 0, int ic0 = 0, int ic1 = 0, int ic2 = 0,
	doublereal order0 = 1.0, doublereal order1 = 1.0, 
	doublereal order2 = 1.0) 
      : m_rxn (rxn), m_ic0 (ic0), m_ic1 (ic1), m_ic2 (ic2),
      m_order0(order0), m_order1(order1), m_order2(order2) {}

    int data(vector<int>& ic) {
      ic.resize(3);
      ic[0] = m_ic0;
      ic[1] = m_ic1;
      ic[2] = m_ic2;
      return m_rxn;
    }

    void multiply(const doublereal* input, doublereal* output) const { 
      *(output + m_rxn) *= (*(input + m_ic0)) * (*(input + m_ic1)) 
	* (*(input + m_ic2)); 
    }
    void power(const doublereal* input, doublereal* output) const { 
      output[m_rxn] *= pow(input[m_ic0],m_order0) * 
	pow(input[m_ic1],m_order1) * 
	pow(input[m_ic2],m_order2);
    }
    void incrementSpecies(const doublereal* input, 
			  doublereal* output) const {
      doublereal x = *(input + m_rxn); 
      *(output + m_ic0) += x;
      *(output + m_ic1) += x;
      *(output + m_ic2) += x;
    }
    void decrementSpecies(const doublereal* input, 
			  doublereal* output) const {
      doublereal x = *(input + m_rxn); 
      *(output + m_ic0) -= x;
      *(output + m_ic1) -= x; 
      *(output + m_ic2) -= x; 
    }
    void incrementReaction(const doublereal* input, 
			   doublereal* output) const { 
      *(output + m_rxn) += *(input + m_ic0) + *(input + m_ic1) 
	+ *(input + m_ic2); 
    }
    void decrementReaction(const doublereal* input, 
			   doublereal* output) const { 
      *(output + m_rxn) -= (*(input + m_ic0) + *(input + m_ic1) 
			    + *(input + m_ic2)); 
    }
  private:
    int m_rxn, m_ic0, m_ic1, m_ic2;    
    doublereal m_order0, m_order1, m_order2;
  };


  /**
   * Handles any number of species in a reaction.
   * @ingroup Stoichiometry
   */
  class C_AnyN {
  public:
    C_AnyN() : m_rxn (-1) {}

    C_AnyN( int rxn, const vector_int& ic, const vector_fp& order) 
      : m_rxn (rxn) {
      m_n = ic.size();
      m_ic.resize(m_n);
      m_order.resize(m_n);
      for (int n = 0; n < m_n; n++) {
	m_ic[n] = ic[n];
	m_order[n] = order[n];
      }
    }

    int data(vector<int>& ic) {
      ic.resize(m_n);
      int n;
      for (n = 0; n < m_n; n++) ic[n] = m_ic[n];
      return m_rxn;
    }

    void power(const doublereal* input, doublereal* output) const {
      for (int n = 0; n < m_n; n++) output[m_rxn] 
	*= pow(input[m_ic[n]],m_order[n]); 
    }

    void multiply(const doublereal* input, doublereal* output) const {
      for (int n = 0; n < m_n; n++) output[m_rxn] *= input[m_ic[n]]; 
    }

    void incrementSpecies(const doublereal* input, 
			  doublereal* output) const {
      doublereal x = input[m_rxn];
      for (int n = 0; n < m_n; n++) output[m_ic[n]] += x;
    }
    void decrementSpecies(const doublereal* input, 
			  doublereal* output) const {
      doublereal x = input[m_rxn];
      for (int n = 0; n < m_n; n++) output[m_ic[n]] -= x;
    }
    void incrementReaction(const doublereal* input, 
			   doublereal* output) const { 
      for (int n = 0; n < m_n; n++) output[m_rxn] += input[m_ic[n]];
    }
    void decrementReaction(const doublereal* input, 
			   doublereal* output) const { 
      for (int n = 0; n < m_n; n++) output[m_rxn] -= input[m_ic[n]];
    }
  private:
    int m_n, m_rxn, m_ic0, m_ic1, m_ic2;
    vector_int m_ic;
    vector_fp m_order;
  };


  template<class _InputIter, class Vec1, class Vec2>
    inline static void _multiply(_InputIter __begin, _InputIter __end, 
				 const Vec1& input, Vec2& output) {
    for (; __begin != __end; ++__begin) 
      __begin->multiply(input, output);
  }
  template<class _InputIter, class Vec1, class Vec2>
    inline static void _power(_InputIter __begin, _InputIter __end, 
			      const Vec1& input, Vec2& output) {
    for (; __begin != __end; ++__begin) 
      __begin->power(input, output);
  }
  template<class _InputIter, class Vec1, class Vec2>
    inline static void _incrementSpecies(_InputIter __begin, 
					 _InputIter __end, const Vec1& input, Vec2& output) {
    for (; __begin != __end; ++__begin) 
      __begin->incrementSpecies(input, output);
  }
  template<class _InputIter, class Vec1, class Vec2>
    inline static void _decrementSpecies(_InputIter __begin, 
					 _InputIter __end, const Vec1& input, Vec2& output) {
    for (; __begin != __end; ++__begin) 
      __begin->decrementSpecies(input, output);
  }
  template<class _InputIter, class Vec1, class Vec2>
    inline static void _incrementReactions(_InputIter __begin, 
					   _InputIter __end, const Vec1& input, Vec2& output) {
    for (; __begin != __end; ++__begin) 
      __begin->incrementReaction(input, output);
  }
  template<class _InputIter, class Vec1, class Vec2>
    inline static void _decrementReactions(_InputIter __begin, 
					   _InputIter __end, const Vec1& input, Vec2& output) {
    for (; __begin != __end; ++__begin) 
      __begin->decrementReaction(input, output);
  }


  /**
   * This class handles operations involving the stoichiometric
   * coefficients on one side of a reaction (reactant or product) for
   * a set of reactions comprising a reaction mechanism. This class is
   * used by class ReactionStoichMgr, which contains three instances
   * of this class (one to handle operations on the reactions, one for
   * the products of reversible reactions, and one for the products of
   * irreversible reactions).
   *
   * This class is designed for use with elementary reactions, or at
   * least ones with integral stoichiometric coefficients. Let \f$ M(i) \f$
   * be the number of molecules on the product or reactant side of
   * reaction number i.
   * \f[
   * r_i = \sum_m^{M_i} s_{k_{m,i}}
   * \f] 
   * To understand the operations performed by this class, let
   * $N_{k,i}$ denote the stoichiometric coefficient of species k on
   * one side (reactant or product) in reaction i. Then \b N is a sparse
   * K by I matrix of stoichiometric coefficients.
   * 
   * The following matrix operations may be carried out with a vector
   * S of length K, and a vector R of length I:
   * 
   * - \f$ S = S + N R\f$   (incrementSpecies)
   * - \f$ S = S - N R\f$   (decrementSpecies)
   * - \f$ R = R + N^T S \f$ (incrementReaction)
   * - \f$ R = R - N^T S \f$ (deccrementReaction)
   *
   * The actual implementation, however, does not compute these
   * quantities by matrix multiplication. A faster algorithm is used
   * that makes use of the fact that the \b integer-valued N matrix is
   * very sparse, and the non-zero terms are small positive integers.
   * \f[
   * S_k = R_{i1} + \dots + R_{iM}
   * \f]
   * where M is the number of molecules, and $\f i(m) \f$ is the 
   * @ingroup Stoichiometry
   */
  class StoichManagerN {
  public:

    /**
     * Constructor for the StoichManagerN class.
     *
     * @internal Consider adding defaulted entries here that supply
     * the total number of reactions in the mechanism and the total
     * number of species in the species list. Then, we could use those
     * numbers to provide error checks during the construction of the
     * object. Those numbers would also provide some clarity to the
     * purpose and utility of this class.
     *
     * DGG - the problem is that the number of reactions and species
     * are not known initially.
     */
    StoichManagerN()  {}

    /**
     * Add a single reaction to the list of reactions that this
     * stoichiometric manager object handles.
     *
     * This function is the same as the add() function below. However,
     * the order of each species in the power list expression is
     * set to one automatically.
     */
    void add(int rxn, const vector_int& k) {
      vector_fp order(k.size(), 1.0);
      add(rxn, k, order);
    }

    /**
     * Add a single reaction to the list of reactions that this
     * stoichiometric manager object handles.
     *
     * @param rxn  Reaction index of the current reaction. This is used
     *             as an index into vectors which have length n_total_rxn.
     * @param k    This is a vector of integer values specifying the
     *             species indecises. The length of this vector species
     *             the number of different species in the description.
     *             The value of the entries are the species indices.
     *             These are used as indexes into vectors which have
     *             length n_total_species.
     *  @param order This is a vector of the same length as vector k.
     *         The order is used for the routine order(), which produces
     *         a power law expression involving the species vector.
     */
    void add(int rxn, const vector_int& k, const vector_fp& order) {
      m_n[rxn] = k.size();
      switch (k.size()) {
      case 1:
	m_loc[rxn] = m_c1_list.size();
	m_c1_list.push_back(C1(rxn, k[0], order[0])); 
	break; 
      case 2:
	m_loc[rxn] = m_c2_list.size();  
	m_c2_list.push_back(C2(rxn, k[0], k[1], order[0], order[1])); 
	break; 
      case 3:
	m_loc[rxn] = m_c3_list.size();  
	m_c3_list.push_back(C3(rxn, k[0], k[1], k[2], 
			       order[0], order[1], order[2])); 
	break; 
      default:
	m_loc[rxn] = m_cn_list.size(); 
	m_cn_list.push_back(C_AnyN(rxn, k, order));
      }
    }
    
    void multiply(const doublereal* input, doublereal* output) const {
      _multiply(m_c1_list.begin(), m_c1_list.end(), input, output);
      _multiply(m_c2_list.begin(), m_c2_list.end(), input, output);
      _multiply(m_c3_list.begin(), m_c3_list.end(), input, output);
      _multiply(m_cn_list.begin(), m_cn_list.end(), input, output);
    }

    void power(const doublereal* input, doublereal* output) const {
      _power(m_c1_list.begin(), m_c1_list.end(), input, output);
      _power(m_c2_list.begin(), m_c2_list.end(), input, output);
      _power(m_c3_list.begin(), m_c3_list.end(), input, output);
      _power(m_cn_list.begin(), m_cn_list.end(), input, output);
    }
    
    void incrementSpecies(const doublereal* input, doublereal* output) const {
      _incrementSpecies(m_c1_list.begin(), m_c1_list.end(), input, output);
      _incrementSpecies(m_c2_list.begin(), m_c2_list.end(), input, output);
      _incrementSpecies(m_c3_list.begin(), m_c3_list.end(), input, output);
      _incrementSpecies(m_cn_list.begin(), m_cn_list.end(), input, output);
    }

    void decrementSpecies(const doublereal* input, doublereal* output) const {
      _decrementSpecies(m_c1_list.begin(), m_c1_list.end(), input, output);
      _decrementSpecies(m_c2_list.begin(), m_c2_list.end(), input, output);
      _decrementSpecies(m_c3_list.begin(), m_c3_list.end(), input, output);
      _decrementSpecies(m_cn_list.begin(), m_cn_list.end(), input, output);
    }

    void incrementReactions(const doublereal* input, doublereal* output) const {
      _incrementReactions(m_c1_list.begin(), m_c1_list.end(), input, output);
      _incrementReactions(m_c2_list.begin(), m_c2_list.end(), input, output);
      _incrementReactions(m_c3_list.begin(), m_c3_list.end(), input, output);
      _incrementReactions(m_cn_list.begin(), m_cn_list.end(), input, output);
    }

    void decrementReactions(const doublereal* input, doublereal* output) const {
      _decrementReactions(m_c1_list.begin(), m_c1_list.end(), input, output);
      _decrementReactions(m_c2_list.begin(), m_c2_list.end(), input, output);
      _decrementReactions(m_c3_list.begin(), m_c3_list.end(), input, output);
      _decrementReactions(m_cn_list.begin(), m_cn_list.end(), input, output);
    }
 
  private:

    vector<C1>     m_c1_list;
    vector<C2>     m_c2_list;
    vector<C3>     m_c3_list;
    vector<C_AnyN> m_cn_list;
    /**
     * Mapping with the Reaction Number as key and the Number of species
     * as the value.
     */
    map<int, int>  m_n;
    /**
     * Mapping with the Reaction Number as key and the placement in the
     * vector of reactions list( i.e., m_c1_list[]) as key
     */
    map<int, int>  m_loc;
  };




  class StoichWriter {
  public:

    StoichWriter()  {}

    void add(int rxn, const vector_int& k) {
      int n, nn = k.size();
      for (n = 0; n < nn; n++) {
	if (m_mult[rxn] != "") m_mult[rxn] += " * ";
	m_mult[rxn] += "c[" + int2str(k[n]) + "]";
	m_is[k[n]] += " + rop[" + int2str(rxn) + "]";
	m_ds[k[n]] += " - rop[" + int2str(rxn) + "]";
	m_ir[rxn] += " + grt[" + int2str(k[n]) + "]";
	m_dr[rxn] += " - grt[" + int2str(k[n]) + "]";
      }
    }

    void writeIncSpec(ostream& s, int nsp) {
      int k; 
      for (k = 0; k < nsp; k++) {
	s << "out[" << k << "] = " << m_is[k] << ";" << endl;
      }
    }

    string mult(int rxn) { return m_mult[rxn]; }
    string incrSpec(int k) { return m_is[k]; }
    string decrSpec(int k) { return m_ds[k]; }
    string incrRxn(int rxn) { return m_ir[rxn]; }    
    string decrRxn(int rxn) { return m_dr[rxn]; }

  private:
    map<int, string> m_mult, m_ir, m_dr, m_is, m_ds;
  };
}

#endif
