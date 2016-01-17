/**
 *  @file StoichManager.h
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_STOICH_MGR_H
#define CT_STOICH_MGR_H

#include "cantera/base/stringUtils.h"

namespace Cantera
{

/**
 * @defgroup Stoichiometry Stoichiometry
 *
 * Note: these classes are designed for internal use in class
 * ReactionStoichManager.
 *
 * The classes defined here implement simple operations that are used by class
 * ReactionStoichManager to compute things like rates of progress, species
 * production rates, etc. In general, a reaction mechanism may involve many
 * species and many reactions, but any given reaction typically only involves
 * a few species as reactants, and a few as products. Therefore, the matrix of
 * stoichiometric coefficients is very sparse. Not only is it sparse, but the
 * non-zero matrix elements often have the value 1, and in many cases no more
 * than three coefficients are non-zero for the reactants and/or the products.
 *
 * For the present purposes, we will consider each direction of a reversible
 * reaction to be a separate reaction. We often need to compute quantities
 * that can formally be written as a matrix product of a stoichiometric
 * coefficient matrix and a vector of reaction rates. For example, the species
 * creation rates are given by
 *
 * \f[
 *  \dot C_k = \sum_k \nu^{(p)}_{k,i} R_i
 * \f]
 *
 * where \f$ \nu^{(p)_{k,i}} \f$ is the product-side stoichiometric
 * coefficient of species \a k in reaction \a i. This could be done be
 * straightforward matrix multiplication, but would be inefficient, since most
 * of the matrix elements of \f$ \nu^{(p)}_{k,i} \f$ are zero. We could do
 * better by using sparse-matrix algorithms to compute this product.
 *
 * If the reactions are general ones, with non-integral stoichiometric
 * coefficients, this is about as good as we can do. But we are particularly
 * concerned here with the performance for very large reaction mechanisms,
 * which are usually composed of elementary reactions, which have integral
 * stoichiometric coefficients. Furthermore, very few elementary reactions
 * involve more than 3 product or reactant molecules.
 *
 * But we can do even better if we take account of the special structure of
 * this matrix for elementary reactions involving three or fewer product
 * molecules (or reactant molecules).
 *
 * To take advantage of this structure, reactions are divided into four
 * groups. These classes are designed to take advantage of this sparse
 * structure when computing quantities that can be written as matrix
 * multiplies.
 *
 * They are designed to explicitly unroll loops over species or reactions for
 * Operations on reactions that require knowing the reaction stoichiometry.
 *
 * This module consists of class StoichManager, and classes C1, C2, and C3.
 * Classes C1, C2, and C3 handle operations involving one, two, or three
 * species, respectively, in a reaction. Instances are instantiated with a
 * reaction number, and n species numbers (n = 1 for C1, etc.).  All three
 * classes have the same interface.
 *
 * These classes are designed for use by StoichManager, and the operations
 * implemented are those needed to efficiently compute quantities such as
 * rates of progress, species production rates, reaction thermochemistry, etc.
 * The compiler will inline these methods into the body of the corresponding
 * StoichManager method, and so there is no performance penalty (unless
 * inlining is turned off).
 *
 * To describe the methods, consider class C3 and suppose an instance is
 * created with reaction number irxn and species numbers k0, k1, and k2.
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
 * The function multiply() is usually used when evaluating the forward and
 * reverse rates of progress of reactions. The rate constants are usually
 * loaded into out[]. Then multiply() is called to add in the dependence of
 * the species concentrations to yield a forward and reverse rop.
 *
 * The function incrementSpecies() and its cousin decrementSpecies() is used
 * to translate from rates of progress to species production rates. The vector
 * in[] is preloaded with the rates of progress of all reactions. Then
 * incrementSpecies() is called to increment the species production vector,
 * out[], with the rates of progress.
 *
 * The functions incrementReaction() and decrementReaction() are used to find
 * the standard state equilibrium constant for a reaction. Here, output[] is a
 * vector of length number of reactions, usually the standard gibbs free
 * energies of reaction, while input, usually the standard state gibbs free
 * energies of species, is a vector of length number of species.
 *
 * Note the stoichiometric coefficient for a species in a reaction is handled
 * by always assuming it is equal to one and then treating reactants and
 * products for a reaction separately. Bimolecular reactions involving the
 * identical species are treated as involving separate species.
 *
 * @internal This class should be upgraded to include cases where
 * real stoichiometric coefficients are used. Shouldn't be that
 * hard to do, and they occur in engineering simulations with some
 * regularity.
 *
 */

static doublereal ppow(doublereal x, doublereal order)
{
    if (x > 0.0) {
        return std::pow(x, order);
    } else {
        return 0.0;
    }
}

inline static std::string fmt(const std::string& r, size_t n)
{
    return r + "[" + int2str(n) + "]";
}

/**
 * Handles one species in a reaction.
 * See @ref Stoichiometry
 * @ingroup Stoichiometry
 * @internal
 */
class C1
{
public:
    C1(size_t rxn = 0, size_t ic0 = 0) :
        m_rxn(rxn),
        m_ic0(ic0) {
    }

    C1(const C1& right) :
        m_rxn(right.m_rxn),
        m_ic0(right.m_ic0) {
    }

    C1& operator=(const C1& right) {
        if (this != &right) {
            m_rxn = right.m_rxn;
            m_ic0 = right.m_ic0;
        }
        return *this;
    }

    size_t data(std::vector<size_t>& ic) {
        ic.resize(1);
        ic[0] = m_ic0;
        return m_rxn;
    }

    void incrementSpecies(const doublereal* R, doublereal* S) const {
        S[m_ic0] += R[m_rxn];
    }

    void decrementSpecies(const doublereal* R, doublereal* S) const {
        S[m_ic0] -= R[m_rxn];
    }

    void multiply(const doublereal* S, doublereal* R) const {
        R[m_rxn] *= S[m_ic0];
    }

    void incrementReaction(const doublereal* S, doublereal* R) const {
        R[m_rxn] += S[m_ic0];
    }

    void decrementReaction(const doublereal* S, doublereal* R) const {
        R[m_rxn] -= S[m_ic0];
    }

    size_t rxnNumber() const {
        return m_rxn;
    }
    size_t speciesIndex(size_t n) const {
        return m_ic0;
    }
    size_t nSpecies() {
        return 1;
    }

    void writeMultiply(const std::string& r, std::map<size_t, std::string>& out) {
        out[m_rxn] = fmt(r, m_ic0);
    }

    void writeIncrementReaction(const std::string& r, std::map<size_t, std::string>& out) {
        out[m_rxn] += " + "+fmt(r, m_ic0);
    }
    void writeDecrementReaction(const std::string& r, std::map<size_t, std::string>& out) {
        out[m_rxn] += " - "+fmt(r, m_ic0);
    }

    void writeIncrementSpecies(const std::string& r, std::map<size_t, std::string>& out) {
        out[m_ic0] += " + "+fmt(r, m_rxn);
    }
    void writeDecrementSpecies(const std::string& r, std::map<size_t, std::string>& out) {
        out[m_ic0] += " - "+fmt(r, m_rxn);
    }

private:
    //! Reaction number
    size_t m_rxn;
    //! Species number
    size_t m_ic0;
};

/**
 * Handles two species in a single reaction.
 * See @ref Stoichiometry
 * @ingroup Stoichiometry
 */
class C2
{
public:
    C2(size_t rxn = 0, size_t ic0 = 0, size_t ic1 = 0)
        : m_rxn(rxn), m_ic0(ic0), m_ic1(ic1) {}

    C2(const C2& right) :
        m_rxn(right.m_rxn),
        m_ic0(right.m_ic0),
        m_ic1(right.m_ic1) {
    }

    C2& operator=(const C2& right) {
        if (this != &right) {
            m_rxn = right.m_rxn;
            m_ic0 = right.m_ic0;
            m_ic1 = right.m_ic1;
        }
        return *this;
    }

    size_t data(std::vector<size_t>& ic) {
        ic.resize(2);
        ic[0] = m_ic0;
        ic[1] = m_ic1;
        return m_rxn;
    }

    void incrementSpecies(const doublereal* R, doublereal* S) const {
        S[m_ic0] += R[m_rxn];
        S[m_ic1] += R[m_rxn];
    }

    void decrementSpecies(const doublereal* R, doublereal* S) const {
        S[m_ic0] -= R[m_rxn];
        S[m_ic1] -= R[m_rxn];
    }

    void multiply(const doublereal* S, doublereal* R) const {
        R[m_rxn] *= S[m_ic0] * S[m_ic1];
    }

    void incrementReaction(const doublereal* S, doublereal* R) const {
        R[m_rxn] += S[m_ic0] + S[m_ic1];
    }

    void decrementReaction(const doublereal* S, doublereal* R) const {
        R[m_rxn] -= (S[m_ic0] + S[m_ic1]);
    }

    size_t rxnNumber() const {
        return m_rxn;
    }
    size_t speciesIndex(size_t n) const {
        return (n == 0 ? m_ic0 : m_ic1);
    }
    size_t nSpecies() {
        return 2;
    }

    void writeMultiply(const std::string& r, std::map<size_t, std::string>& out) {
        out[m_rxn] = fmt(r, m_ic0) + " * " + fmt(r, m_ic1);
    }
    void writeIncrementReaction(const std::string& r, std::map<size_t, std::string>& out) {
        out[m_rxn] += " + "+fmt(r, m_ic0)+" + "+fmt(r, m_ic1);
    }
    void writeDecrementReaction(const std::string& r, std::map<size_t, std::string>& out) {
        out[m_rxn] += " - "+fmt(r, m_ic0)+" - "+fmt(r, m_ic1);
    }

    void writeIncrementSpecies(const std::string& r, std::map<size_t, std::string>& out) {
        std::string s = " + "+fmt(r, m_rxn);
        out[m_ic0] += s;
        out[m_ic1] += s;
    }
    void writeDecrementSpecies(const std::string& r, std::map<size_t, std::string>& out) {
        std::string s = " - "+fmt(r, m_rxn);
        out[m_ic0] += s;
        out[m_ic1] += s;
    }

private:
    //! Reaction index -> index into the ROP vector
    size_t m_rxn;

    //! Species index -> index into the species vector for the two species.
    size_t m_ic0, m_ic1;
};

/**
 * Handles three species in a reaction.
 * See @ref Stoichiometry
 * @ingroup Stoichiometry
 */
class C3
{
public:
    C3(size_t rxn = 0, size_t ic0 = 0, size_t ic1 = 0, size_t ic2 = 0)
        : m_rxn(rxn), m_ic0(ic0), m_ic1(ic1), m_ic2(ic2) {}

    C3(const C3& right) :
        m_rxn(right.m_rxn),
        m_ic0(right.m_ic0),
        m_ic1(right.m_ic1),
        m_ic2(right.m_ic2) {
    }

    C3& operator=(const C3& right) {
        if (this != &right) {
            m_rxn = right.m_rxn;
            m_ic0 = right.m_ic0;
            m_ic1 = right.m_ic1;
            m_ic2 = right.m_ic2;
        }
        return *this;
    }

    size_t data(std::vector<size_t>& ic) {
        ic.resize(3);
        ic[0] = m_ic0;
        ic[1] = m_ic1;
        ic[2] = m_ic2;
        return m_rxn;
    }

    void incrementSpecies(const doublereal* R, doublereal* S) const {
        S[m_ic0] += R[m_rxn];
        S[m_ic1] += R[m_rxn];
        S[m_ic2] += R[m_rxn];
    }

    void decrementSpecies(const doublereal* R, doublereal* S) const {
        S[m_ic0] -= R[m_rxn];
        S[m_ic1] -= R[m_rxn];
        S[m_ic2] -= R[m_rxn];
    }

    void multiply(const doublereal* S, doublereal* R) const {
        R[m_rxn] *= S[m_ic0] * S[m_ic1] * S[m_ic2];
    }

    void incrementReaction(const doublereal* S, doublereal* R) const {
        R[m_rxn] += S[m_ic0] + S[m_ic1] + S[m_ic2];
    }

    void decrementReaction(const doublereal* S, doublereal* R) const {
        R[m_rxn] -= (S[m_ic0] + S[m_ic1] + S[m_ic2]);
    }

    size_t rxnNumber() const {
        return m_rxn;
    }
    size_t speciesIndex(size_t n) const {
        return (n == 0 ? m_ic0 : (n == 1 ? m_ic1 : m_ic2));
    }
    size_t nSpecies() {
        return 3;
    }

    void writeMultiply(const std::string& r, std::map<size_t, std::string>& out) {
        out[m_rxn] = fmt(r, m_ic0) + " * " + fmt(r, m_ic1) + " * " + fmt(r, m_ic2);
    }
    void writeIncrementReaction(const std::string& r, std::map<size_t, std::string>& out) {
        out[m_rxn] += " + "+fmt(r, m_ic0)+" + "+fmt(r, m_ic1)+" + "+fmt(r, m_ic2);
    }
    void writeDecrementReaction(const std::string& r, std::map<size_t, std::string>& out) {
        out[m_rxn] += " - "+fmt(r, m_ic0)+" - "+fmt(r, m_ic1)+" - "+fmt(r, m_ic2);
    }
    void writeIncrementSpecies(const std::string& r, std::map<size_t, std::string>& out) {
        std::string s = " + "+fmt(r, m_rxn);
        out[m_ic0] += s;
        out[m_ic1] += s;
        out[m_ic2] += s;
    }
    void writeDecrementSpecies(const std::string& r, std::map<size_t, std::string>& out) {
        std::string s = " - "+fmt(r, m_rxn);
        out[m_ic0] += s;
        out[m_ic1] += s;
        out[m_ic2] += s;
    }
private:
    size_t m_rxn;
    size_t m_ic0;
    size_t m_ic1;
    size_t m_ic2;
};

/**
 * Handles any number of species in a reaction, including fractional
 * stoichiometric coefficients, and arbitrary reaction orders.
 * See @ref Stoichiometry
 * @ingroup Stoichiometry
 */
class C_AnyN
{
public:
    C_AnyN() :
        m_n(0),
        m_rxn(npos) {
    }

    C_AnyN(size_t rxn, const std::vector<size_t>& ic, const vector_fp& order_, const vector_fp& stoich_) :
        m_n(0),
        m_rxn(rxn) {
        m_n = ic.size();
        m_ic.resize(m_n);
        m_order.resize(m_n);
        m_stoich.resize(m_n);
        for (size_t n = 0; n < m_n; n++) {
            m_ic[n] = ic[n];
            m_order[n] = order_[n];
            m_stoich[n] = stoich_[n];
        }
    }

    C_AnyN(const C_AnyN& right) :
        m_n(right.m_n),
        m_rxn(right.m_rxn),
        m_ic(right.m_ic),
        m_order(right.m_order),
        m_stoich(right.m_stoich) {
    }

    C_AnyN& operator=(const C_AnyN& right) {
        if (this != &right) {
            m_n = right.m_n;
            m_rxn = right.m_rxn;
            m_ic = right.m_ic;
            m_order = right.m_order;
            m_stoich = right.m_stoich;
        }
        return *this;
    }

    size_t data(std::vector<size_t>& ic) {
        ic.resize(m_n);
        for (size_t n = 0; n < m_n; n++) {
            ic[n] = m_ic[n];
        }
        return m_rxn;
    }

    doublereal order(size_t n) const {
        return m_order[n];
    }
    doublereal stoich(size_t n) const {
        return m_stoich[n];
    }
    size_t speciesIndex(size_t n) const {
        return m_ic[n];
    }

    void multiply(const doublereal* input, doublereal* output) const {
        doublereal oo;
        for (size_t n = 0; n < m_n; n++) {
            oo = m_order[n];
            if (oo != 0.0) {
                output[m_rxn] *= ppow(input[m_ic[n]], oo);
            }
        }
    }

    void incrementSpecies(const doublereal* input,
                          doublereal* output) const {
        doublereal x = input[m_rxn];
        for (size_t n = 0; n < m_n; n++) {
            output[m_ic[n]] += m_stoich[n]*x;
        }
    }

    void decrementSpecies(const doublereal* input,
                          doublereal* output) const {
        doublereal x = input[m_rxn];
        for (size_t n = 0; n < m_n; n++) {
            output[m_ic[n]] -= m_stoich[n]*x;
        }
    }

    void incrementReaction(const doublereal* input,
                           doublereal* output) const {
        for (size_t n = 0; n < m_n; n++) output[m_rxn]
            += m_stoich[n]*input[m_ic[n]];
    }

    void decrementReaction(const doublereal* input,
                           doublereal* output) const {
        for (size_t n = 0; n < m_n; n++) output[m_rxn]
            -= m_stoich[n]*input[m_ic[n]];
    }

    void writeMultiply(const std::string& r, std::map<size_t, std::string>& out) {
        out[m_rxn] = "";
        for (size_t n = 0; n < m_n; n++) {
            if (m_order[n] == 1.0) {
                out[m_rxn] += fmt(r, m_ic[n]);
            } else {
                out[m_rxn] += "pow("+fmt(r, m_ic[n])+","+fp2str(m_order[n])+")";
            }
            if (n < m_n-1) {
                out[m_rxn] += " * ";
            }
        }
    }
    void writeIncrementReaction(const std::string& r, std::map<size_t, std::string>& out) {
        for (size_t n = 0; n < m_n; n++) {
            out[m_rxn] += " + "+fp2str(m_stoich[n]) + "*" + fmt(r, m_ic[n]);
        }
    }
    void writeDecrementReaction(const std::string& r, std::map<size_t, std::string>& out) {
        for (size_t n = 0; n < m_n; n++) {
            out[m_rxn] += " - "+fp2str(m_stoich[n]) + "*" + fmt(r, m_ic[n]);
        }
    }
    void writeIncrementSpecies(const std::string& r, std::map<size_t, std::string>& out) {
        std::string s = fmt(r, m_rxn);
        for (size_t n = 0; n < m_n; n++) {
            out[m_ic[n]] += " + "+fp2str(m_stoich[n]) + "*" + s;
        }
    }

    void writeDecrementSpecies(const std::string& r, std::map<size_t, std::string>& out) {
        std::string s = fmt(r, m_rxn);
        for (size_t n = 0; n < m_n; n++) {
            out[m_ic[n]] += " - "+fp2str(m_stoich[n]) + "*" + s;
        }
    }

private:
    //! Length of the m_ic vector
    /*!
     *   This is the number of species which have non-zero entries in either the
     *   reaction order matrix or the stoichiometric order matrix for this reaction.
     */
    size_t m_n;

    //!  ID of the reaction corresponding to this stoichiometric manager
    /*!
     *  This is used within the interface to select the
     */
    size_t m_rxn;

    //! Vector of species which are involved with this stoichiometric manager calculations
    /*!
     *  This is an integer list of species which have non-zero entries in either the
     *  reaction order matrix or the stoichiometric order matrix for this reaction, m_rxn.
     *  It's used as the index into the arrays m_order[] and m_stoich[].
     */
    std::vector<size_t> m_ic;
    vector_fp m_order;
    vector_fp m_stoich;
};

template<class InputIter, class Vec1, class Vec2>
inline static void _multiply(InputIter begin, InputIter end,
                             const Vec1& input, Vec2& output)
{
    for (; begin != end; ++begin) {
        begin->multiply(input, output);
    }
}

template<class InputIter, class Vec1, class Vec2>
inline static void _incrementSpecies(InputIter begin,
                                     InputIter end, const Vec1& input, Vec2& output)
{
    for (; begin != end; ++begin) {
        begin->incrementSpecies(input, output);
    }
}

template<class InputIter, class Vec1, class Vec2>
inline static void _decrementSpecies(InputIter begin,
                                     InputIter end, const Vec1& input, Vec2& output)
{
    for (; begin != end; ++begin) {
        begin->decrementSpecies(input, output);
    }
}

template<class InputIter, class Vec1, class Vec2>
inline static void _incrementReactions(InputIter begin,
                                       InputIter end, const Vec1& input, Vec2& output)
{
    for (; begin != end; ++begin) {
        begin->incrementReaction(input, output);
    }
}

template<class InputIter, class Vec1, class Vec2>
inline static void _decrementReactions(InputIter begin,
                                       InputIter end, const Vec1& input, Vec2& output)
{
    for (; begin != end; ++begin) {
        begin->decrementReaction(input, output);
    }
}


template<class InputIter>
inline static void _writeIncrementSpecies(InputIter begin, InputIter end,
        const std::string& r, std::map<size_t, std::string>& out)
{
    for (; begin != end; ++begin) {
        begin->writeIncrementSpecies(r, out);
    }
}

template<class InputIter>
inline static void _writeDecrementSpecies(InputIter begin, InputIter end,
        const std::string& r, std::map<size_t, std::string>& out)
{
    for (; begin != end; ++begin) {
        begin->writeDecrementSpecies(r, out);
    }
}

template<class InputIter>
inline static void _writeIncrementReaction(InputIter begin, InputIter end,
        const std::string& r, std::map<size_t, std::string>& out)
{
    for (; begin != end; ++begin) {
        begin->writeIncrementReaction(r, out);
    }
}

template<class InputIter>
inline static void _writeDecrementReaction(InputIter begin, InputIter end,
        const std::string& r, std::map<size_t, std::string>& out)
{
    for (; begin != end; ++begin) {
        begin->writeDecrementReaction(r, out);
    }
}

template<class InputIter>
inline static void _writeMultiply(InputIter begin, InputIter end,
                                  const std::string& r, std::map<size_t, std::string>& out)
{
    for (; begin != end; ++begin) {
        begin->writeMultiply(r, out);
    }
}

/*
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
 * \f$ N_{k,i}\f$ denote the stoichiometric coefficient of species k on
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
 * See @ref Stoichiometry
 * @ingroup Stoichiometry
 */
class StoichManagerN
{
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
    StoichManagerN() {
    }

    StoichManagerN(const StoichManagerN& right) :
        m_c1_list(right.m_c1_list),
        m_c2_list(right.m_c2_list),
        m_c3_list(right.m_c3_list),
        m_cn_list(right.m_cn_list),
        m_n(right.m_n),
        m_loc(right.m_loc) {
    }

    StoichManagerN& operator=(const StoichManagerN& right) {
        if (this != &right) {
            m_c1_list = right.m_c1_list;
            m_c2_list = right.m_c2_list;
            m_c3_list = right.m_c3_list;
            m_cn_list = right.m_cn_list;
            m_n = right.m_n;
            m_loc = right.m_loc;
        }
        return *this;
    }

    /**
     * Add a single reaction to the list of reactions that this
     * stoichiometric manager object handles.
     *
     * This function is the same as the add() function below. However,
     * the order of each species in the power list expression is
     * set to one automatically.
     */
    void add(size_t rxn, const std::vector<size_t>& k) {
        vector_fp order(k.size(), 1.0);
        vector_fp stoich(k.size(), 1.0);
        add(rxn, k, order, stoich);
    }

    void add(size_t rxn, const std::vector<size_t>& k, const vector_fp& order) {
        vector_fp stoich(k.size(), 1.0);
        add(rxn, k, order, stoich);
    }

    //! Add a single reaction to the list of reactions that this
    //! stoichiometric manager object handles.
    /*!
     * @param rxn  Reaction index of the current reaction. This is used
     *             as an index into vectors which have length n_total_rxn.
     * @param k    This is a vector of integer values specifying the
     *             species indices. The length of this vector species
     *             the number of different species in the description.
     *             The value of the entries are the species indices.
     *             These are used as indexes into vectors which have
     *             length n_total_species.
     *  @param order This is a vector of the same length as vector k.
     *         The order is used for the routine power(), which produces
     *         a power law expression involving the species vector.
     *  @param stoich  This is used to handle fractional stoichiometric coefficients
     *                 on the product side of irreversible reactions.
     */
    void add(size_t rxn, const std::vector<size_t>& k, const vector_fp& order,
             const vector_fp& stoich) {
        m_n[rxn] = k.size();
        bool frac = false;
        for (size_t n = 0; n < stoich.size(); n++) {
            if (stoich[n] != 1.0 || order[n] != 1.0) {
                frac = true;
                break;
            }
        }
        if (frac) {
            m_loc[rxn] = m_cn_list.size();
            m_cn_list.push_back(C_AnyN(rxn, k, order, stoich));
        } else {
            switch (k.size()) {
            case 1:
                m_loc[rxn] = m_c1_list.size();
                m_c1_list.push_back(C1(rxn, k[0]));
                break;
            case 2:
                m_loc[rxn] = m_c2_list.size();
                m_c2_list.push_back(C2(rxn, k[0], k[1]));
                break;
            case 3:
                m_loc[rxn] = m_c3_list.size();
                m_c3_list.push_back(C3(rxn, k[0], k[1], k[2]));
                break;
            default:
                m_loc[rxn] = m_cn_list.size();
                m_cn_list.push_back(C_AnyN(rxn, k, order, stoich));
            }
        }
    }

    void multiply(const doublereal* input, doublereal* output) const {
        _multiply(m_c1_list.begin(), m_c1_list.end(), input, output);
        _multiply(m_c2_list.begin(), m_c2_list.end(), input, output);
        _multiply(m_c3_list.begin(), m_c3_list.end(), input, output);
        _multiply(m_cn_list.begin(), m_cn_list.end(), input, output);
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

    void writeIncrementSpecies(const std::string& r, std::map<size_t, std::string>& out) {
        _writeIncrementSpecies(m_c1_list.begin(), m_c1_list.end(), r, out);
        _writeIncrementSpecies(m_c2_list.begin(), m_c2_list.end(), r, out);
        _writeIncrementSpecies(m_c3_list.begin(), m_c3_list.end(), r, out);
        _writeIncrementSpecies(m_cn_list.begin(), m_cn_list.end(), r, out);
    }

    void writeDecrementSpecies(const std::string& r, std::map<size_t, std::string>& out) {
        _writeDecrementSpecies(m_c1_list.begin(), m_c1_list.end(), r, out);
        _writeDecrementSpecies(m_c2_list.begin(), m_c2_list.end(), r, out);
        _writeDecrementSpecies(m_c3_list.begin(), m_c3_list.end(), r, out);
        _writeDecrementSpecies(m_cn_list.begin(), m_cn_list.end(), r, out);
    }

    void writeIncrementReaction(const std::string& r, std::map<size_t, std::string>& out) {
        _writeIncrementReaction(m_c1_list.begin(), m_c1_list.end(), r, out);
        _writeIncrementReaction(m_c2_list.begin(), m_c2_list.end(), r, out);
        _writeIncrementReaction(m_c3_list.begin(), m_c3_list.end(), r, out);
        _writeIncrementReaction(m_cn_list.begin(), m_cn_list.end(), r, out);
    }

    void writeDecrementReaction(const std::string& r, std::map<size_t, std::string>& out) {
        _writeDecrementReaction(m_c1_list.begin(), m_c1_list.end(), r, out);
        _writeDecrementReaction(m_c2_list.begin(), m_c2_list.end(), r, out);
        _writeDecrementReaction(m_c3_list.begin(), m_c3_list.end(), r, out);
        _writeDecrementReaction(m_cn_list.begin(), m_cn_list.end(), r, out);
    }

    void writeMultiply(const std::string& r, std::map<size_t, std::string>& out) {
        _writeMultiply(m_c1_list.begin(), m_c1_list.end(), r, out);
        _writeMultiply(m_c2_list.begin(), m_c2_list.end(), r, out);
        _writeMultiply(m_c3_list.begin(), m_c3_list.end(), r, out);
        _writeMultiply(m_cn_list.begin(), m_cn_list.end(), r, out);
    }

private:
    std::vector<C1>     m_c1_list;
    std::vector<C2>     m_c2_list;
    std::vector<C3>     m_c3_list;
    std::vector<C_AnyN> m_cn_list;
    /**
     * Map with the Reaction Number as key and the Number of species
     * as the value.
     */
    std::map<size_t, size_t>  m_n;
    /**
     * Map with the Reaction Number as key and the placement in the
     * vector of reactions list( i.e., m_c1_list[]) as key
     */
    std::map<size_t, size_t>  m_loc;
};

}

#endif
