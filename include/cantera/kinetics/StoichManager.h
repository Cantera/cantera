/**
 *  @file StoichManager.h
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_STOICH_MGR_H
#define CT_STOICH_MGR_H

#include "cantera/numerics/eigen_sparse.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

/**
 * @defgroup Stoichiometry Stoichiometry
 *
 * The classes defined here implement simple operations that are used by class
 * Kinetics to compute things like rates of progress, species
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
 * coefficient of species \a k in reaction \a i. This could be done by
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
 * To take advantage of this structure, reactions are divided into four groups.
 * These classes are designed to take advantage of this sparse structure when
 * computing quantities that can be written as matrix multiplies.
 *
 * They are designed to explicitly unroll loops over species or reactions for
 * operations on reactions that require knowing the reaction stoichiometry.
 *
 * This module consists of class StoichManagerN, and classes C1, C2, and C3.
 * Classes C1, C2, and C3 handle operations involving one, two, or three
 * species, respectively, in a reaction. Instances are instantiated with a
 * reaction number, and n species numbers (n = 1 for C1, etc.). All three
 * classes have the same interface.
 *
 * These classes are designed for use by StoichManagerN, and the operations
 * implemented are those needed to efficiently compute quantities such as
 * rates of progress, species production rates, reaction thermochemistry, etc.
 * The compiler will inline these methods into the body of the corresponding
 * StoichManagerN method, and so there is no performance penalty (unless
 * inlining is turned off).
 *
 * To describe the methods, consider class C3 and suppose an instance is
 * created with reaction number irxn and species numbers k0, k1, and k2.
 *
 *  - multiply(in, out) : out[irxn] is multiplied by
 *    in[k0] * in[k1] * in[k2]
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
 * vector of length number of reactions, usually the standard Gibbs free
 * energies of reaction, while input, usually the standard state Gibbs free
 * energies of species, is a vector of length number of species.
 *
 * Note the stoichiometric coefficient for a species in a reaction is handled
 * by always assuming it is equal to one and then treating reactants and
 * products for a reaction separately. Bimolecular reactions involving the
 * identical species are treated as involving separate species.
 *
 * The methods resizeCoeffs(), derivatives() and scale() are used for the calculation
 * of derivatives with respect to species mole fractions. In this context,
 * resizeCoeffs() is used to establish a mapping between a reaction and corresponding
 * non-zero entries of the sparse derivative matrix of the reaction rate-of-progress
 * vector, which itself is evaluated by the derivatives() method. The scale() method is
 * used to multiply rop entries by reaction order and a user-supplied factor.
 */

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

    void resizeCoeffs(const std::map<std::pair<size_t, size_t>, size_t>& indices)
    {
        m_jc0 = indices.at({m_rxn, m_ic0});
    }

    void derivatives(const double* S, const double* R, vector_fp& jac) const
    {
        jac[m_jc0] += R[m_rxn]; // index (m_ic0, m_rxn)
    }


    void scale(const double* R, double* out, double factor) const
    {
        out[m_rxn] = R[m_rxn] * factor;
    }

private:
    //! Reaction number
    size_t m_rxn;
    //! Species number
    size_t m_ic0;

    size_t m_jc0; //!< Index in derivative triplet vector
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

    void incrementSpecies(const doublereal* R, doublereal* S) const {
        S[m_ic0] += R[m_rxn];
        S[m_ic1] += R[m_rxn];
    }

    void decrementSpecies(const doublereal* R, doublereal* S) const {
        S[m_ic0] -= R[m_rxn];
        S[m_ic1] -= R[m_rxn];
    }

    void multiply(const doublereal* S, doublereal* R) const {
        if (S[m_ic0] < 0 && S[m_ic1] < 0) {
            R[m_rxn] = 0;
        } else {
            R[m_rxn] *= S[m_ic0] * S[m_ic1];
        }
    }

    void incrementReaction(const doublereal* S, doublereal* R) const {
        R[m_rxn] += S[m_ic0] + S[m_ic1];
    }

    void decrementReaction(const doublereal* S, doublereal* R) const {
        R[m_rxn] -= (S[m_ic0] + S[m_ic1]);
    }

    void resizeCoeffs(const std::map<std::pair<size_t, size_t>, size_t>& indices)
    {
        m_jc0 = indices.at({m_rxn, m_ic0});
        m_jc1 = indices.at({m_rxn, m_ic1});
    }

    void derivatives(const double* S, const double* R, vector_fp& jac) const
    {
        if (S[m_ic1] > 0) {
            jac[m_jc0] += R[m_rxn] * S[m_ic1]; // index (m_ic0, m_rxn)
        }
        if (S[m_ic0] > 0) {
            jac[m_jc1] += R[m_rxn] * S[m_ic0]; // index (m_ic1, m_rxn)
        }
    }

    void scale(const double* R, double* out, double factor) const
    {
        out[m_rxn] = 2 * R[m_rxn] * factor;
    }

private:
    //! Reaction index -> index into the ROP vector
    size_t m_rxn;

    //! Species index -> index into the species vector for the two species.
    size_t m_ic0, m_ic1;

    size_t m_jc0, m_jc1; //!< Indices in derivative triplet vector
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
        if ((S[m_ic0] < 0 && (S[m_ic1] < 0 || S[m_ic2] < 0)) ||
            (S[m_ic1] < 0 && S[m_ic2] < 0)) {
            R[m_rxn] = 0;
        } else {
            R[m_rxn] *= S[m_ic0] * S[m_ic1] * S[m_ic2];
        }
    }

    void incrementReaction(const doublereal* S, doublereal* R) const {
        R[m_rxn] += S[m_ic0] + S[m_ic1] + S[m_ic2];
    }

    void decrementReaction(const doublereal* S, doublereal* R) const {
        R[m_rxn] -= (S[m_ic0] + S[m_ic1] + S[m_ic2]);
    }

    void resizeCoeffs(const std::map<std::pair<size_t, size_t>, size_t>& indices)
    {
        m_jc0 = indices.at({m_rxn, m_ic0});
        m_jc1 = indices.at({m_rxn, m_ic1});
        m_jc2 = indices.at({m_rxn, m_ic2});
    }

    void derivatives(const double* S, const double* R, vector_fp& jac) const
    {
        if (S[m_ic1] > 0 && S[m_ic2] > 0) {
            jac[m_jc0] += R[m_rxn] * S[m_ic1] * S[m_ic2];; // index (m_ic0, m_rxn)
        }
        if (S[m_ic0] > 0 && S[m_ic2] > 0) {
            jac[m_jc1] += R[m_rxn] * S[m_ic0] * S[m_ic2]; // index (m_ic1, m_ic1)
        }
        if (S[m_ic0] > 0 && S[m_ic1] > 0) {
            jac[m_jc2] += R[m_rxn] * S[m_ic0] * S[m_ic1]; // index (m_ic2, m_ic2)
        }
    }

    void scale(const double* R, double* out, double factor) const
    {
        out[m_rxn] = 3 * R[m_rxn] * factor;
    }

private:
    size_t m_rxn;
    size_t m_ic0;
    size_t m_ic1;
    size_t m_ic2;

    size_t m_jc0, m_jc1, m_jc2; //!< Indices in derivative triplet vector
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
        m_jc.resize(m_n, 0);
        m_order.resize(m_n);
        m_stoich.resize(m_n);
        for (size_t n = 0; n < m_n; n++) {
            m_ic[n] = ic[n];
            m_order[n] = order_[n];
            m_stoich[n] = stoich_[n];
        }
    }

    void multiply(const doublereal* input, doublereal* output) const {
        for (size_t n = 0; n < m_n; n++) {
            double order = m_order[n];
            if (order != 0.0) {
                double c = input[m_ic[n]];
                if (c > 0.0) {
                    output[m_rxn] *= std::pow(c, order);
                } else {
                    output[m_rxn] = 0.0;
                }
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
        for (size_t n = 0; n < m_n; n++) {
            output[m_rxn] += m_stoich[n]*input[m_ic[n]];
        }
    }

    void decrementReaction(const doublereal* input,
                           doublereal* output) const {
        for (size_t n = 0; n < m_n; n++) {
            output[m_rxn] -= m_stoich[n]*input[m_ic[n]];
        }
    }

    void resizeCoeffs(const std::map<std::pair<size_t, size_t>, size_t>& indices)
    {
        for (size_t i = 0; i < m_n; i++) {
            m_jc[i] = indices.at({m_rxn, m_ic[i]});
        }

        m_sum_order = 0.;
        for (size_t n = 0; n < m_n; ++n) {
            m_sum_order += m_order[n];
        }
    }

    void derivatives(const double* S, const double* R, vector_fp& jac) const
    {
        for (size_t i = 0; i < m_n; i++) {
            // calculate derivative
            double prod = R[m_rxn];
            double order_i = m_order[i];
            if (S[m_ic[i]] > 0. && order_i != 0.) {
                prod *= order_i * std::pow(S[m_ic[i]], order_i - 1);
                for (size_t j = 0; j < m_n; j++) {
                    if (i != j) {
                        if (S[m_ic[j]] > 0) {
                            prod *= std::pow(S[m_ic[j]], m_order[j]);
                        } else {
                            prod = 0.;
                            break;
                        }
                    }
                }
            } else {
                prod = 0.;
            }
            jac[m_jc[i]] += prod;
        }
    }

    void scale(const double* R, double* out, double factor) const
    {
        out[m_rxn] = m_sum_order * R[m_rxn] * factor;
    }

private:
    //! Length of the m_ic vector
    /*!
     *  This is the number of species which participate in the reaction order
     *  and stoichiometric coefficient vectors for the reactant or product
     *  description of the reaction.
     */
    size_t m_n;

    //! ID of the reaction corresponding to this stoichiometric manager
    /*!
     *  This is used within the interface to select the array position to read
     *  and write to Normally this is associated with the reaction number in an
     *  array of quantities indexed by the reaction number, for example, ROP[irxn].
     */
    size_t m_rxn;

    //! Vector of species which are involved with this stoichiometric manager
    //! calculations
    /*!
     *  This is an integer list of species which participate in either the
     *  reaction order matrix or the stoichiometric order matrix for this
     *  reaction, m_rxn.
     */
    std::vector<size_t> m_ic;

    //! Reaction orders for the reaction
    /*!
     * This is either for the reactants or products. Length = m_n. Species
     * number, m_ic[n], has a reaction order of m_order[n].
     */
    vector_fp m_order;

    //! Stoichiometric coefficients for the reaction, reactant or product side.
    /*!
     *  This is either for the reactants or products. Length = m_n. Species
     *  number m_ic[m], has a stoichiometric coefficient of m_stoich[n].
     */
    vector_fp m_stoich;

    std::vector<size_t> m_jc; //!< Indices in derivative triplet vector

    double m_sum_order; //!< Sum of reaction order vector
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

template<class InputIter, class Indices>
inline static void _resizeCoeffs(InputIter begin, InputIter end, Indices& ix)
{
    for (; begin != end; ++begin) {
        begin->resizeCoeffs(ix);
    }
}

template<class InputIter, class Vec1, class Vec2, class Vec3>
inline static void _derivatives(InputIter begin, InputIter end,
                             const Vec1& S, const Vec2& R, Vec3& jac)
{
    for (; begin != end; ++begin) {
        begin->derivatives(S, R, jac);
    }
}

template<class InputIter, class Vec1, class Vec2>
inline static void _scale(InputIter begin, InputIter end,
                          const Vec1& in, Vec2& out, double factor)
{
    for (; begin != end; ++begin) {
        begin->scale(in, out, factor);
    }
}

/*!
 * This class handles operations involving the stoichiometric coefficients on
 * one side of a reaction (reactant or product) for a set of reactions
 * comprising a reaction mechanism. This class is used by class Kinetics, which
 * contains three instances of this class (one to handle operations on the
 * reactions, one for the products of reversible reactions, and one for the
 * products of irreversible reactions).
 *
 * This class is designed for use with elementary reactions, or at least ones
 * with integral stoichiometric coefficients. Let \f$ M(i) \f$ be the number of
 * molecules on the product or reactant side of reaction number i.
 * \f[
 *     r_i = \sum_m^{M_i} s_{k_{m,i}}
 * \f]
 * To understand the operations performed by this class, let \f$ N_{k,i}\f$
 * denote the stoichiometric coefficient of species k on one side (reactant or
 * product) in reaction i. Then \b N is a sparse K by I matrix of stoichiometric
 * coefficients.
 *
 * The following matrix operations may be carried out with a vector S of length
 * K, and a vector R of length I:
 *
 * - \f$ S = S + N R\f$   (incrementSpecies)
 * - \f$ S = S - N R\f$   (decrementSpecies)
 * - \f$ R = R + N^T S \f$ (incrementReaction)
 * - \f$ R = R - N^T S \f$ (decrementReaction)
 *
 * The actual implementation, however, does not compute these quantities by
 * matrix multiplication. A faster algorithm is used that makes use of the fact
 * that the \b integer-valued N matrix is very sparse, and the non-zero terms
 * are small positive integers.
 * \f[
 * S_k = R_{i1} + \dots + R_{iM}
 * \f]
 * where M is the number of molecules, and \f$ i(m) \f$ is the
 * See @ref Stoichiometry
 * @ingroup Stoichiometry
 */
class StoichManagerN
{
public:
    /**
     * Constructor for the StoichManagerN class.
     *
     * @internal Consider adding defaulted entries here that supply the total
     *     number of reactions in the mechanism and the total number of species
     *     in the species list. Then, we could use those numbers to provide
     *     error checks during the construction of the object. Those numbers
     *     would also provide some clarity to the purpose and utility of this
     *     class.
     *
     * DGG - the problem is that the number of reactions and species are not
     * known initially.
     */
    StoichManagerN() : m_ready(true) {
        m_stoichCoeffs.setZero();
        m_stoichCoeffs.resize(0, 0);
    }

    //! Resize the sparse coefficient matrix
    void resizeCoeffs(size_t nSpc, size_t nRxn)
    {
        size_t nCoeffs = m_coeffList.size();

        // Stoichiometric coefficient matrix
        m_stoichCoeffs.resize(nSpc, nRxn);
        m_stoichCoeffs.reserve(nCoeffs);
        m_stoichCoeffs.setFromTriplets(m_coeffList.begin(), m_coeffList.end());

        // Set up outer/inner indices for mapped derivative output
        Eigen::SparseMatrix<double> tmp = m_stoichCoeffs.transpose();
        m_outerIndices.resize(nSpc + 1); // number of columns + 1
        for (int i = 0; i < tmp.outerSize() + 1; i++) {
            m_outerIndices[i] = tmp.outerIndexPtr()[i];
        }
        m_innerIndices.resize(nCoeffs);
        for (size_t n = 0; n < nCoeffs; n++) {
            m_innerIndices[n] = tmp.innerIndexPtr()[n];
        }
        m_values.resize(nCoeffs, 0.);

        // Set up index pairs for derivatives
        std::map<std::pair<size_t, size_t>, size_t> indices;
        size_t n = 0;
        for (int i = 0; i < tmp.outerSize(); i++) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(tmp, i); it; ++it) {
                indices[{static_cast<size_t>(it.row()),
                    static_cast<size_t>(it.col())}] = n++;
            }
        }
        // update reaction setup
        _resizeCoeffs(m_c1_list.begin(), m_c1_list.end(), indices);
        _resizeCoeffs(m_c2_list.begin(), m_c2_list.end(), indices);
        _resizeCoeffs(m_c3_list.begin(), m_c3_list.end(), indices);
        _resizeCoeffs(m_cn_list.begin(), m_cn_list.end(), indices);

        m_ready = true;
    }

    /**
     * Add a single reaction to the list of reactions that this stoichiometric
     * manager object handles.
     *
     * This function is the same as the add() function below. However, the order
     * of each species in the power list expression is set to one automatically.
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
     * @param rxn  Reaction index of the current reaction. This is used as an
     *     index into vectors which have length n_total_rxn.
     * @param k    This is a vector of integer values specifying the species
     *     indices. The length of this vector species the number of different
     *     species in the description. The value of the entries are the species
     *     indices. These are used as indexes into vectors which have length
     *     n_total_species.
     *  @param order This is a vector of the same length as vector k. The order
     *     is used for the routine power(), which produces a power law
     *     expression involving the species vector.
     *  @param stoich  This is used to handle fractional stoichiometric
     *     coefficients on the product side of irreversible reactions.
     */
    void add(size_t rxn, const std::vector<size_t>& k, const vector_fp& order,
             const vector_fp& stoich) {
        if (order.size() != k.size()) {
           throw CanteraError("StoichManagerN::add()", "size of order and species arrays differ");
        }
        if (stoich.size() != k.size()) {
           throw CanteraError("StoichManagerN::add()", "size of stoich and species arrays differ");
        }
        bool frac = false;
        for (size_t n = 0; n < stoich.size(); n++) {
            m_coeffList.emplace_back(
                static_cast<int>(k[n]), static_cast<int>(rxn), stoich[n]);
            if (fmod(stoich[n], 1.0) || stoich[n] != order[n]) {
                frac = true;
            }
        }
        if (frac || k.size() > 3) {
            m_cn_list.emplace_back(rxn, k, order, stoich);
        } else {
            // Try to express the reaction with unity stoichiometric
            // coefficients (by repeating species when necessary) so that the
            // simpler 'multiply' function can be used to compute the rate
            // instead of 'power'.
            std::vector<size_t> kRep;
            for (size_t n = 0; n < k.size(); n++) {
                for (size_t i = 0; i < stoich[n]; i++) {
                    kRep.push_back(k[n]);
                }
            }

            switch (kRep.size()) {
            case 1:
                m_c1_list.emplace_back(rxn, kRep[0]);
                break;
            case 2:
                m_c2_list.emplace_back(rxn, kRep[0], kRep[1]);
                break;
            case 3:
                m_c3_list.emplace_back(rxn, kRep[0], kRep[1], kRep[2]);
                break;
            default:
                m_cn_list.emplace_back(rxn, k, order, stoich);
            }
        }
        m_ready = false;
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

    //! Return matrix containing stoichiometric coefficients
    const Eigen::SparseMatrix<double>& stoichCoeffs() const
    {
        if (!m_ready) {
            // This can happen if a user overrides default behavior:
            // Kinetics::resizeReactions is not called after adding reactions via
            // Kinetics::addReaction with the 'resize' flag set to 'false'
            throw CanteraError("StoichManagerN::stoichCoeffs", "The object "
                "is not fully configured; make sure to call resizeCoeffs().");
        }
        return m_stoichCoeffs;
    }

    //! Calculate derivatives with respect to species concentrations.
    /*!
     * The species derivative is the term of the Jacobian that accounts for
     * species products in the law of mass action. As StoichManagerN does not account
     * for third bodies or rate constants that depend on species concentrations,
     * corresponding derivatives are not included.
     *
     *  @param conc    Species concentration.
     *  @param rates   Rates-of-progress.
     */
    Eigen::SparseMatrix<double> derivatives(const double* conc, const double* rates)
    {
        // calculate derivative entries using known sparse storage order
        std::fill(m_values.begin(), m_values.end(), 0.);
        _derivatives(m_c1_list.begin(), m_c1_list.end(), conc, rates, m_values);
        _derivatives(m_c2_list.begin(), m_c2_list.end(), conc, rates, m_values);
        _derivatives(m_c3_list.begin(), m_c3_list.end(), conc, rates, m_values);
        _derivatives(m_cn_list.begin(), m_cn_list.end(), conc, rates, m_values);

        return Eigen::Map<Eigen::SparseMatrix<double>>(
            m_stoichCoeffs.cols(), m_stoichCoeffs.rows(), m_values.size(),
            m_outerIndices.data(), m_innerIndices.data(), m_values.data());
    }

    //! Scale input by reaction order and factor
    void scale(const double* in, double* out, double factor) const
    {
        _scale(m_c1_list.begin(), m_c1_list.end(), in, out, factor);
        _scale(m_c2_list.begin(), m_c2_list.end(), in, out, factor);
        _scale(m_c3_list.begin(), m_c3_list.end(), in, out, factor);
        _scale(m_cn_list.begin(), m_cn_list.end(), in, out, factor);
    }

private:
    bool m_ready; //!< Boolean flag indicating whether object is fully configured

    std::vector<C1> m_c1_list;
    std::vector<C2> m_c2_list;
    std::vector<C3> m_c3_list;
    std::vector<C_AnyN> m_cn_list;

    //! Sparse matrices for stoichiometric coefficients
    SparseTriplets m_coeffList;
    Eigen::SparseMatrix<double> m_stoichCoeffs;

    //! Storage indicies used to build derivatives
    std::vector<int> m_outerIndices;
    std::vector<int> m_innerIndices;
    vector_fp m_values;
};

}

#endif
