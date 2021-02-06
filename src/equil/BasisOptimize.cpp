/**
 * @file BasisOptimize.cpp Functions which calculation optimized basis of the
 *     stoichiometric coefficient matrix (see /ref equil functions)
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/equil/MultiPhase.h"

using namespace std;

namespace Cantera
{
int BasisOptimize_print_lvl = 0;
static const double USEDBEFORE = -1;

size_t BasisOptimize(int* usedZeroedSpecies, bool doFormRxn, MultiPhase* mphase,
                     std::vector<size_t>& orderVectorSpecies,
                     std::vector<size_t>& orderVectorElements,
                     vector_fp& formRxnMatrix)
{
    // Get the total number of elements defined in the multiphase object
    size_t ne = mphase->nElements();

    // Get the total number of species in the multiphase object
    size_t nspecies = mphase->nSpecies();

    // Perhaps, initialize the element ordering
    if (orderVectorElements.size() < ne) {
        orderVectorElements.resize(ne);
        iota(orderVectorElements.begin(), orderVectorElements.end(), 0);
    }

    // Perhaps, initialize the species ordering
    if (orderVectorSpecies.size() != nspecies) {
        orderVectorSpecies.resize(nspecies);
        iota(orderVectorSpecies.begin(), orderVectorSpecies.end(), 0);
    }

    if (BasisOptimize_print_lvl >= 1) {
        writelog("   ");
        writeline('-', 77);
        writelog("   --- Subroutine BASOPT called to ");
        writelog("calculate the number of components and ");
        writelog("evaluate the formation matrix\n");
        if (BasisOptimize_print_lvl > 0) {
            writelog("   ---\n");

            writelog("   ---      Formula Matrix used in BASOPT calculation\n");
            writelog("   ---      Species | Order | ");
            for (size_t j = 0; j < ne; j++) {
                size_t jj = orderVectorElements[j];
                writelog(" {:>4.4s}({:1d})", mphase->elementName(jj), j);
            }
            writelog("\n");
            for (size_t k = 0; k < nspecies; k++) {
                size_t kk = orderVectorSpecies[k];
                writelog("   --- {:>11.11s} |   {:4d} |",
                         mphase->speciesName(kk), k);
                for (size_t j = 0; j < ne; j++) {
                    size_t jj = orderVectorElements[j];
                    double num = mphase->nAtoms(kk,jj);
                    writelogf("%6.1g  ", num);
                }
                writelog("\n");
            }
            writelog("   --- \n");
        }
    }

    // Calculate the maximum value of the number of components possible. It's
    // equal to the minimum of the number of elements and the number of total
    // species.
    size_t nComponents = std::min(ne, nspecies);
    size_t nNonComponents = nspecies - nComponents;

    // Set this return variable to false
    *usedZeroedSpecies = false;

    // Create an array of mole numbers
    vector_fp molNum(nspecies,0.0);
    mphase->getMoles(molNum.data());

    // Other workspace
    DenseMatrix sm(ne, ne);
    vector_fp ss(ne, 0.0);
    vector_fp sa(ne, 0.0);
    if (formRxnMatrix.size() < nspecies*ne) {
        formRxnMatrix.resize(nspecies*ne, 0.0);
    }

    // For debugging purposes keep an unmodified copy of the array.
    vector_fp molNumBase = molNum;
    double molSave = 0.0;
    size_t jr = 0;

    // Top of a loop of some sort based on the index JR. JR is the current
    // number of component species found.
    while (jr < nComponents) {
        // Top of another loop point based on finding a linearly independent
        // species
        size_t k = npos;
        while (true) {
            // Search the remaining part of the mole number vector, molNum for
            // the largest remaining species. Return its identity. kk is the raw
            // number. k is the orderVectorSpecies index.
            size_t kk = max_element(molNum.begin(), molNum.end()) - molNum.begin();
            size_t j;
            for (j = 0; j < nspecies; j++) {
                if (orderVectorSpecies[j] == kk) {
                    k = j;
                    break;
                }
            }
            if (j == nspecies) {
                throw CanteraError("BasisOptimize", "orderVectorSpecies contains an error");
            }

            if (molNum[kk] == 0.0) {
                *usedZeroedSpecies = true;
            }
            // If the largest molNum is negative, then we are done.
            if (molNum[kk] == USEDBEFORE) {
                nComponents = jr;
                nNonComponents = nspecies - nComponents;
                break;
            }

            // Assign a small negative number to the component that we have
            // just found, in order to take it out of further consideration.
            molSave = molNum[kk];
            molNum[kk] = USEDBEFORE;

            // CHECK LINEAR INDEPENDENCE WITH PREVIOUS SPECIES

            // Modified Gram-Schmidt Method, p. 202 Dalquist
            // QR factorization of a matrix without row pivoting.
            size_t jl = jr;
            for (j = 0; j < ne; ++j) {
                size_t jj = orderVectorElements[j];
                sm(j, jr) = mphase->nAtoms(kk,jj);
            }
            if (jl > 0) {
                // Compute the coefficients of JA column of the the upper
                // triangular R matrix, SS(J) = R_J_JR (this is slightly
                // different than Dalquist) R_JA_JA = 1
                for (j = 0; j < jl; ++j) {
                    ss[j] = 0.0;
                    for (size_t i = 0; i < ne; ++i) {
                        ss[j] += sm(i, jr) * sm(i, j);
                    }
                    ss[j] /= sa[j];
                }

                // Now make the new column, (*,JR), orthogonal to the previous
                // columns
                for (j = 0; j < jl; ++j) {
                    for (size_t i = 0; i < ne; ++i) {
                        sm(i, jr) -= ss[j] * sm(i, j);
                    }
                }
            }

            // Find the new length of the new column in Q.
            // It will be used in the denominator in future row calcs.
            sa[jr] = 0.0;
            for (size_t ml = 0; ml < ne; ++ml) {
                sa[jr] += pow(sm(ml, jr), 2);
            }

            // IF NORM OF NEW ROW  .LT. 1E-3 REJECT
            if (sa[jr] > 1.0e-6) {
                break;
            }
        }

        // REARRANGE THE DATA
        if (jr != k) {
            if (BasisOptimize_print_lvl >= 1) {
                size_t kk = orderVectorSpecies[k];
                writelogf("   ---   %-12.12s", mphase->speciesName(kk));
                size_t jj = orderVectorSpecies[jr];
                writelogf("(%9.2g) replaces %-12.12s",
                          molSave, mphase->speciesName(jj));
                writelogf("(%9.2g) as component %3d\n", molNum[jj], jr);
            }
            std::swap(orderVectorSpecies[jr], orderVectorSpecies[k]);
        }

        // If we haven't found enough components, go back and find some more
        jr++;
    }

    if (! doFormRxn) {
        return nComponents;
    }

    // EVALUATE THE STOICHIOMETRY
    //
    // Formulate the matrix problem for the stoichiometric
    // coefficients. CX + B = 0
    //
    // C will be an nc x nc matrix made up of the formula vectors for the
    // components. Each component's formula vector is a column. The rows are the
    // elements.
    //
    // n RHS's will be solved for. Thus, B is an nc x n matrix.
    //
    // BIG PROBLEM 1/21/99:
    //
    // This algorithm makes the assumption that the first nc rows of the formula
    // matrix aren't rank deficient. However, this might not be the case. For
    // example, assume that the first element in FormulaMatrix[] is argon.
    // Assume that no species in the matrix problem actually includes argon.
    // Then, the first row in sm[], below will be identically zero. bleh.
    //
    // What needs to be done is to perform a rearrangement of the ELEMENTS ->
    // i.e. rearrange, FormulaMatrix, sp, and gai, such that the first nc
    // elements form in combination with the nc components create an invertible
    // sm[]. not a small project, but very doable.
    //
    // An alternative would be to turn the matrix problem below into an ne x nc
    // problem, and do QR elimination instead of Gauss-Jordan elimination.
    //
    // Note the rearrangement of elements need only be done once in the problem.
    // It's actually very similar to the top of this program with ne being the
    // species and nc being the elements!!

    sm.resize(nComponents, nComponents);
    for (size_t k = 0; k < nComponents; ++k) {
        size_t kk = orderVectorSpecies[k];
        for (size_t j = 0; j < nComponents; ++j) {
            size_t jj = orderVectorElements[j];
            sm(j, k) = mphase->nAtoms(kk, jj);
        }
    }

    for (size_t i = 0; i < nNonComponents; ++i) {
        size_t k = nComponents + i;
        size_t kk = orderVectorSpecies[k];
        for (size_t j = 0; j < nComponents; ++j) {
            size_t jj = orderVectorElements[j];
            formRxnMatrix[j + i * ne] = - mphase->nAtoms(kk, jj);
        }
    }

    // // Use LU factorization to calculate the reaction matrix
    solve(sm, formRxnMatrix.data(), nNonComponents, ne);

    if (BasisOptimize_print_lvl >= 1) {
        writelog("   ---\n");
        writelogf("   ---  Number of Components = %d\n", nComponents);
        writelog("   ---  Formula Matrix:\n");
        writelog("   ---                      Components:    ");
        for (size_t k = 0; k < nComponents; k++) {
            size_t kk = orderVectorSpecies[k];
            writelogf(" %3d (%3d) ", k, kk);
        }
        writelog("\n   ---                Components Moles:       ");
        for (size_t k = 0; k < nComponents; k++) {
            size_t kk = orderVectorSpecies[k];
            writelogf("%-11.3g", molNumBase[kk]);
        }
        writelog("\n   ---        NonComponent |   Moles  |       ");
        for (size_t i = 0; i < nComponents; i++) {
            size_t kk = orderVectorSpecies[i];
            writelogf("%-11.10s", mphase->speciesName(kk));
        }
        writelog("\n");

        for (size_t i = 0; i < nNonComponents; i++) {
            size_t k = i + nComponents;
            size_t kk = orderVectorSpecies[k];
            writelogf("   --- %3d (%3d) ", k, kk);
            writelogf("%-10.10s", mphase->speciesName(kk));
            writelogf("|%10.3g|", molNumBase[kk]);
            // Print the negative of formRxnMatrix[]; it's easier to interpret.
            for (size_t j = 0; j < nComponents; j++) {
                writelogf("     %6.2f", - formRxnMatrix[j + i * ne]);
            }
            writelog("\n");
        }
        writelog("   ");
        writeline('-', 77);
    }

    return nComponents;
} // basopt()


void ElemRearrange(size_t nComponents, const vector_fp& elementAbundances,
                   MultiPhase* mphase,
                   std::vector<size_t>& orderVectorSpecies,
                   std::vector<size_t>& orderVectorElements)
{
    size_t nelements = mphase->nElements();
    // Get the total number of species in the multiphase object
    size_t nspecies = mphase->nSpecies();

    if (BasisOptimize_print_lvl > 0) {
        writelog("   ");
        writeline('-', 77);
        writelog("   --- Subroutine ElemRearrange() called to ");
        writelog("check stoich. coefficient matrix\n");
        writelog("   ---    and to rearrange the element ordering once\n");
    }

    // Perhaps, initialize the element ordering
    if (orderVectorElements.size() < nelements) {
        orderVectorElements.resize(nelements);
        for (size_t j = 0; j < nelements; j++) {
            orderVectorElements[j] = j;
        }
    }

    // Perhaps, initialize the species ordering. However, this is dangerous, as
    // this ordering is assumed to yield the component species for the problem
    if (orderVectorSpecies.size() != nspecies) {
        orderVectorSpecies.resize(nspecies);
        for (size_t k = 0; k < nspecies; k++) {
            orderVectorSpecies[k] = k;
        }
    }

    // If the elementAbundances aren't input, just create a fake one based on
    // summing the column of the stoich matrix. This will force elements with
    // zero species to the end of the element ordering.
    vector_fp eAbund(nelements,0.0);
    if (elementAbundances.size() != nelements) {
        for (size_t j = 0; j < nelements; j++) {
            eAbund[j] = 0.0;
            for (size_t k = 0; k < nspecies; k++) {
                eAbund[j] += fabs(mphase->nAtoms(k, j));
            }
        }
    } else {
        copy(elementAbundances.begin(), elementAbundances.end(),
             eAbund.begin());
    }

    vector_fp sa(nelements,0.0);
    vector_fp ss(nelements,0.0);
    vector_fp sm(nelements*nelements,0.0);

    // Top of a loop of some sort based on the index JR. JR is the current
    // number independent elements found.
    size_t jr = 0;
    while (jr < nComponents) {
        // Top of another loop point based on finding a linearly independent
        // element
        size_t k = nelements;
        while (true) {
            // Search the element vector. We first locate elements that are
            // present in any amount. Then, we locate elements that are not
            // present in any amount. Return its identity in K.
            k = nelements;
            size_t kk;
            for (size_t ielem = jr; ielem < nelements; ielem++) {
                kk = orderVectorElements[ielem];
                if (eAbund[kk] != USEDBEFORE && eAbund[kk] > 0.0) {
                    k = ielem;
                    break;
                }
            }
            for (size_t ielem = jr; ielem < nelements; ielem++) {
                kk = orderVectorElements[ielem];
                if (eAbund[kk] != USEDBEFORE) {
                    k = ielem;
                    break;
                }
            }

            if (k == nelements) {
                // When we are here, there is an error usually.
                // We haven't found the number of elements necessary.
                if (BasisOptimize_print_lvl > 0) {
                    writelogf("Error exit: returning with nComponents = %d\n", jr);
                }
                throw CanteraError("ElemRearrange", "Required number of elements not found.");
            }

            // Assign a large negative number to the element that we have
            // just found, in order to take it out of further consideration.
            eAbund[kk] = USEDBEFORE;

            // CHECK LINEAR INDEPENDENCE OF CURRENT FORMULA MATRIX
            // LINE WITH PREVIOUS LINES OF THE FORMULA MATRIX

            // Modified Gram-Schmidt Method, p. 202 Dalquist
            // QR factorization of a matrix without row pivoting.
            size_t jl = jr;

            // Fill in the row for the current element, k, under consideration
            // The row will contain the Formula matrix value for that element
            // with respect to the vector of component species. (note j and k
            // indices are flipped compared to the previous routine)
            for (size_t j = 0; j < nComponents; ++j) {
                size_t jj = orderVectorSpecies[j];
                kk = orderVectorElements[k];
                sm[j + jr*nComponents] = mphase->nAtoms(jj,kk);
            }
            if (jl > 0) {
                // Compute the coefficients of JA column of the the upper
                // triangular R matrix, SS(J) = R_J_JR (this is slightly
                // different than Dalquist) R_JA_JA = 1
                for (size_t j = 0; j < jl; ++j) {
                    ss[j] = 0.0;
                    for (size_t i = 0; i < nComponents; ++i) {
                        ss[j] += sm[i + jr*nComponents] * sm[i + j*nComponents];
                    }
                    ss[j] /= sa[j];
                }

                // Now make the new column, (*,JR), orthogonal to the
                // previous columns
                for (size_t j = 0; j < jl; ++j) {
                    for (size_t i = 0; i < nComponents; ++i) {
                        sm[i + jr*nComponents] -= ss[j] * sm[i + j*nComponents];
                    }
                }
            }

            // Find the new length of the new column in Q.
            // It will be used in the denominator in future row calcs.
            sa[jr] = 0.0;
            for (size_t ml = 0; ml < nComponents; ++ml) {
                double tmp = sm[ml + jr*nComponents];
                sa[jr] += tmp * tmp;
            }
            // IF NORM OF NEW ROW  .LT. 1E-6 REJECT
            if (sa[jr] > 1.0e-6) {
                break;
            }
        }
        // REARRANGE THE DATA
        if (jr != k) {
            if (BasisOptimize_print_lvl > 0) {
                size_t kk = orderVectorElements[k];
                writelog("   ---   ");
                writelogf("%-2.2s", mphase->elementName(kk));
                writelog("replaces ");
                kk = orderVectorElements[jr];
                writelogf("%-2.2s", mphase->elementName(kk));
                writelogf(" as element %3d\n", jr);
            }
            std::swap(orderVectorElements[jr], orderVectorElements[k]);
        }

        // If we haven't found enough components, go back and find some more
        jr++;
    };
}

}
