/**
 * @file BasisOptimize.cpp Functions which calculation optimized basis of the
 *     stoichiometric coefficient matrix (see /ref equil functions)
 */
#include "cantera/equil/MultiPhase.h"
#include "cantera/numerics/ctlapack.h"

using namespace Cantera;
using namespace std;

namespace Cantera
{
int BasisOptimize_print_lvl = 0;
}

//! Print a string within a given space limit.
/*!
 *  This routine limits the amount of the string that will be printed to a
 *  maximum of "space" characters.
 *    @param       str        String -> must be null terminated.
 *    @param       space      space limit for the printing.
 *    @param       alignment  0 centered
 *                            1 right aligned
 *                            2 left aligned
 */
static void print_stringTrunc(const char* str, int space, int alignment);

size_t Cantera::BasisOptimize(int* usedZeroedSpecies, bool doFormRxn,
                              MultiPhase* mphase, std::vector<size_t>& orderVectorSpecies,
                              std::vector<size_t>& orderVectorElements,
                              vector_fp& formRxnMatrix)
{
    size_t  j, jj, k=0, kk, l, i, jl, ml;
    bool lindep;
    std::string ename;
    std::string sname;
    /*
     * Get the total number of elements defined in the multiphase object
     */
    size_t ne = mphase->nElements();
    /*
     * Get the total number of species in the multiphase object
     */
    size_t nspecies = mphase->nSpecies();
    doublereal tmp;
    doublereal const USEDBEFORE = -1;

    /*
     * Perhaps, initialize the element ordering
     */
    if (orderVectorElements.size() < ne) {
        orderVectorElements.resize(ne);
        for (j = 0; j < ne; j++) {
            orderVectorElements[j] = j;
        }
    }

    /*
     * Perhaps, initialize the species ordering
     */
    if (orderVectorSpecies.size() != nspecies) {
        orderVectorSpecies.resize(nspecies);
        for (k = 0; k < nspecies; k++) {
            orderVectorSpecies[k] = k;
        }
    }

    if (DEBUG_MODE_ENABLED && BasisOptimize_print_lvl >= 1) {
        writelog("   ");
        for (i=0; i<77; i++) {
            writelog("-");
        }
        writelog("\n");
        writelog("   --- Subroutine BASOPT called to ");
        writelog("calculate the number of components and ");
        writelog("evaluate the formation matrix\n");
        if (BasisOptimize_print_lvl > 0) {
            writelog("   ---\n");

            writelog("   ---      Formula Matrix used in BASOPT calculation\n");
            writelog("   ---      Species | Order | ");
            for (j = 0; j < ne; j++) {
                jj = orderVectorElements[j];
                writelog(" ");
                ename = mphase->elementName(jj);
                print_stringTrunc(ename.c_str(), 4, 1);
                writelogf("(%1d)", j);
            }
            writelog("\n");
            for (k = 0; k < nspecies; k++) {
                kk = orderVectorSpecies[k];
                writelog("   --- ");
                sname = mphase->speciesName(kk);
                print_stringTrunc(sname.c_str(), 11, 1);
                writelogf(" |   %4d |", k);
                for (j = 0; j < ne; j++) {
                    jj = orderVectorElements[j];
                    double num = mphase->nAtoms(kk,jj);
                    writelogf("%6.1g  ", num);
                }
                writelog("\n");
            }
            writelog("   --- \n");
        }
    }

    /*
     *  Calculate the maximum value of the number of components possible
     *     It's equal to the minimum of the number of elements and the
     *     number of total species.
     */
    size_t nComponents = std::min(ne, nspecies);
    size_t nNonComponents = nspecies - nComponents;
    /*
     * Set this return variable to false
     */
    *usedZeroedSpecies = false;

    /*
     * Create an array of mole numbers
     */
    vector_fp molNum(nspecies,0.0);
    mphase->getMoles(DATA_PTR(molNum));

    /*
     * Other workspace
     */
    vector_fp sm(ne*ne, 0.0);
    vector_fp ss(ne, 0.0);
    vector_fp sa(ne, 0.0);
    if (formRxnMatrix.size() < nspecies*ne) {
        formRxnMatrix.resize(nspecies*ne, 0.0);
    }

    /*
     * For debugging purposes keep an unmodified copy of the array.
     */
    vector_fp molNumBase;
    if (DEBUG_MODE_ENABLED) {
        molNumBase = molNum;
    }
    double molSave = 0.0;

    size_t jr = npos;
    /*
     *   Top of a loop of some sort based on the index JR. JR is the
     *   current number of component species found.
     */
    do {
        ++jr;
        /* - Top of another loop point based on finding a linearly */
        /* - independent species */
        do {
            /*
             *    Search the remaining part of the mole number vector, molNum
             *    for the largest remaining species. Return its identity.
             *    kk is the raw number. k is the orderVectorSpecies index.
             */
            kk = max_element(molNum.begin(), molNum.end()) - molNum.begin();
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
            /*
             * If the largest molNum is negative, then we are done.
             */
            if (molNum[kk] == USEDBEFORE) {
                nComponents = jr;
                nNonComponents = nspecies - nComponents;
                break;
            }
            /*
             *  Assign a small negative number to the component that we have
             *  just found, in order to take it out of further consideration.
             */
#ifdef DEBUG_MODE
            molSave = molNum[kk];
#endif
            molNum[kk] = USEDBEFORE;

            /* *********************************************************** */
            /* **** CHECK LINEAR INDEPENDENCE WITH PREVIOUS SPECIES ****** */
            /* *********************************************************** */
            /*
             *          Modified Gram-Schmidt Method, p. 202 Dalquist
             *          QR factorization of a matrix without row pivoting.
             */
            jl = jr;
            for (j = 0; j < ne; ++j) {
                jj = orderVectorElements[j];
                sm[j + jr*ne] = mphase->nAtoms(kk,jj);
            }
            if (jl > 0) {
                /*
                 *         Compute the coefficients of JA column of the
                 *         the upper triangular R matrix, SS(J) = R_J_JR
                 *         (this is slightly different than Dalquist)
                 *         R_JA_JA = 1
                 */
                for (j = 0; j < jl; ++j) {
                    ss[j] = 0.0;
                    for (i = 0; i < ne; ++i) {
                        ss[j] += sm[i + jr*ne] * sm[i + j*ne];
                    }
                    ss[j] /= sa[j];
                }
                /*
                 *     Now make the new column, (*,JR), orthogonal to the
                 *     previous columns
                 */
                for (j = 0; j < jl; ++j) {
                    for (l = 0; l < ne; ++l) {
                        sm[l + jr*ne] -= ss[j] * sm[l + j*ne];
                    }
                }
            }
            /*
             *        Find the new length of the new column in Q.
             *        It will be used in the denominator in future row calcs.
             */
            sa[jr] = 0.0;
            for (ml = 0; ml < ne; ++ml) {
                tmp = sm[ml + jr*ne];
                sa[jr] += tmp * tmp;
            }
            /* **************************************************** */
            /* **** IF NORM OF NEW ROW  .LT. 1E-3 REJECT ********** */
            /* **************************************************** */
            if (sa[jr] < 1.0e-6) {
                lindep = true;
            } else {
                lindep = false;
            }
        } while (lindep);
        /* ****************************************** */
        /* **** REARRANGE THE DATA ****************** */
        /* ****************************************** */
        if (jr != k) {
            if (DEBUG_MODE_ENABLED && BasisOptimize_print_lvl >= 1) {
                kk = orderVectorSpecies[k];
                sname = mphase->speciesName(kk);
                writelogf("   ---   %-12.12s", sname.c_str());
                jj = orderVectorSpecies[jr];
                ename = mphase->speciesName(jj);
                writelogf("(%9.2g) replaces %-12.12s", molSave, ename.c_str());
                writelogf("(%9.2g) as component %3d\n", molNum[jj], jr);
            }
            std::swap(orderVectorSpecies[jr], orderVectorSpecies[k]);
        }

        /*
         *      If we haven't found enough components, go back
         *      and find some more. (nc -1 is used below, because
         *      jr is counted from 0, via the C convention.
         */
    } while (jr < (nComponents-1));


    if (! doFormRxn) {
        return nComponents;
    }

    /* ****************************************************** */
    /* **** EVALUATE THE STOICHIOMETRY ********************** */
    /* ****************************************************** */
    /*
     *  Formulate the matrix problem for the stoichiometric
     *  coefficients. CX + B = 0
     *      C will be an nc x nc matrix made up of the formula
     * vectors for the components. Each component's formula
     * vector is a column. The rows are the elements.
     *      n RHS's will be solved for. Thus, B is an nc x n
     * matrix.
     *
     * BIG PROBLEM 1/21/99:
     *
     *    This algorithm makes the assumption that the
     * first nc rows of the formula matrix aren't rank deficient.
     * However, this might not be the case. For example, assume
     * that the first element in FormulaMatrix[] is argon. Assume that
     * no species in the matrix problem actually includes argon.
     * Then, the first row in sm[], below will be identically
     * zero. bleh.
     *    What needs to be done is to perform a rearrangement
     * of the ELEMENTS -> i.e. rearrange, FormulaMatrix, sp, and gai, such
     * that the first nc elements form in combination with the
     * nc components create an invertible sm[]. not a small
     * project, but very doable.
     *    An alternative would be to turn the matrix problem
     * below into an ne x nc problem, and do QR elimination instead
     * of Gauss-Jordan elimination.
     *    Note the rearrangement of elements need only be done once
     * in the problem. It's actually very similar to the top of
     * this program with ne being the species and nc being the
     * elements!!
     */
    for (k = 0; k < nComponents; ++k) {
        kk = orderVectorSpecies[k];
        for (j = 0; j < nComponents; ++j) {
            jj = orderVectorElements[j];
            sm[j + k*ne] = mphase->nAtoms(kk, jj);
        }
    }

    for (i = 0; i < nNonComponents; ++i) {
        k = nComponents + i;
        kk = orderVectorSpecies[k];
        for (j = 0; j < nComponents; ++j) {
            jj = orderVectorElements[j];
            formRxnMatrix[j + i * ne] = - mphase->nAtoms(kk, jj);
        }
    }
    // Use LU factorization to calculate the reaction matrix
    int info;
    vector_int ipiv(nComponents);
    ct_dgetrf(nComponents, nComponents, &sm[0], ne, &ipiv[0], info);
    if (info) {
        throw CanteraError("basopt", "factorization returned an error condition");
    }
    ct_dgetrs(ctlapack::NoTranspose, nComponents, nNonComponents, &sm[0], ne,
              &ipiv[0], &formRxnMatrix[0], ne, info);

    if (DEBUG_MODE_ENABLED && Cantera::BasisOptimize_print_lvl >= 1) {
        writelog("   ---\n");
        writelogf("   ---  Number of Components = %d\n", nComponents);
        writelog("   ---  Formula Matrix:\n");
        writelog("   ---                      Components:    ");
        for (k = 0; k < nComponents; k++) {
            kk = orderVectorSpecies[k];
            writelogf(" %3d (%3d) ", k, kk);
        }
        writelog("\n   ---                Components Moles:       ");
        for (k = 0; k < nComponents; k++) {
            kk = orderVectorSpecies[k];
            writelogf("%-11.3g", molNumBase[kk]);
        }
        writelog("\n   ---        NonComponent |   Moles  |       ");
        for (i = 0; i < nComponents; i++) {
            kk = orderVectorSpecies[i];
            sname = mphase->speciesName(kk);
            writelogf("%-11.10s", sname.c_str());
        }
        writelog("\n");

        for (i = 0; i < nNonComponents; i++) {
            k = i + nComponents;
            kk = orderVectorSpecies[k];
            writelogf("   --- %3d (%3d) ", k, kk);
            sname = mphase->speciesName(kk);
            writelogf("%-10.10s", sname.c_str());
            writelogf("|%10.3g|", molNumBase[kk]);
            /*
             * Print the negative of formRxnMatrix[]; it's easier to interpret.
             */
            for (j = 0; j < nComponents; j++) {
                writelogf("     %6.2f", - formRxnMatrix[j + i * ne]);
            }
            writelog("\n");
        }
        writelog("   ");
        for (i=0; i<77; i++) {
            writelog("-");
        }
        writelog("\n");
    }

    return nComponents;
} /* basopt() ************************************************************/

static void print_stringTrunc(const char* str, int space, int alignment)

/***********************************************************************
 *  vcs_print_stringTrunc():
 *
 *     Print a string within a given space limit. This routine
 *     limits the amount of the string that will be printed to a
 *     maximum of "space" characters.
 *
 *     str = String -> must be null terminated.
 *     space = space limit for the printing.
 *     alignment = 0 centered
 *           1 right aligned
 *           2 left aligned
 ***********************************************************************/
{
    int i, ls=0, rs=0;
    int len = static_cast<int>(strlen(str));
    if ((len) >= space) {
        for (i = 0; i < space; i++) {
            writelogf("%c", str[i]);
        }
    } else {
        if (alignment == 1) {
            ls = space - len;
        } else if (alignment == 2) {
            rs = space - len;
        } else {
            ls = (space - len) / 2;
            rs = space - len - ls;
        }
        if (ls != 0) {
            for (i = 0; i < ls; i++) {
                writelog(" ");
            }
        }
        writelogf("%s", str);
        if (rs != 0) {
            for (i = 0; i < rs; i++) {
                writelog(" ");
            }
        }
    }
}

size_t Cantera::ElemRearrange(size_t nComponents, const vector_fp& elementAbundances,
                              MultiPhase* mphase,
                              std::vector<size_t>& orderVectorSpecies,
                              std::vector<size_t>& orderVectorElements)
{
    size_t j, k, l, i, jl, ml, jr, ielem, jj, kk=0;

    bool lindep = false;
    size_t nelements = mphase->nElements();
    std::string ename;
    /*
     * Get the total number of species in the multiphase object
     */
    size_t nspecies = mphase->nSpecies();

    double test = -1.0E10;
    if (DEBUG_MODE_ENABLED && BasisOptimize_print_lvl > 0) {
        writelog("   ");
        for (i=0; i<77; i++) {
            writelog("-");
        }
        writelog("\n");
        writelog("   --- Subroutine ElemRearrange() called to ");
        writelog("check stoich. coefficient matrix\n");
        writelog("   ---    and to rearrange the element ordering once\n");
    }

    /*
     * Perhaps, initialize the element ordering
     */
    if (orderVectorElements.size() < nelements) {
        orderVectorElements.resize(nelements);
        for (j = 0; j < nelements; j++) {
            orderVectorElements[j] = j;
        }
    }

    /*
     * Perhaps, initialize the species ordering. However, this is
     * dangerous, as this ordering is assumed to yield the
     * component species for the problem
     */
    if (orderVectorSpecies.size() != nspecies) {
        orderVectorSpecies.resize(nspecies);
        for (k = 0; k < nspecies; k++) {
            orderVectorSpecies[k] = k;
        }
    }

    /*
     * If the elementAbundances aren't input, just create a fake one
     * based on summing the column of the stoich matrix.
     * This will force elements with zero species to the
     * end of the element ordering.
     */
    vector_fp eAbund(nelements,0.0);
    if (elementAbundances.size() != nelements) {
        for (j = 0; j < nelements; j++) {
            eAbund[j] = 0.0;
            for (k = 0; k < nspecies; k++) {
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

    /*
     *        Top of a loop of some sort based on the index JR. JR is the
     *       current number independent elements found.
     */
    jr = npos;
    do {
        ++jr;
        /*
         *     Top of another loop point based on finding a linearly
         *     independent element
         */
        do {
            /*
             *    Search the element vector. We first locate elements that
             *    are present in any amount. Then, we locate elements that
             *    are not present in any amount.
             *    Return its identity in K.
             */
            k = nelements;
            for (ielem = jr; ielem < nelements; ielem++) {
                kk = orderVectorElements[ielem];
                if (eAbund[kk] != test && eAbund[kk] > 0.0) {
                    k = ielem;
                    break;
                }
            }
            for (ielem = jr; ielem < nelements; ielem++) {
                kk = orderVectorElements[ielem];
                if (eAbund[kk] != test) {
                    k = ielem;
                    break;
                }
            }

            if (k == nelements) {
                // When we are here, there is an error usually.
                // We haven't found the number of elements necessary.
                if (DEBUG_MODE_ENABLED && BasisOptimize_print_lvl > 0) {
                    writelogf("Error exit: returning with nComponents = %d\n", jr);
                }
                throw CanteraError("ElemRearrange", "Required number of elements not found.");
            }

            /*
             *  Assign a large negative number to the element that we have
             *  just found, in order to take it out of further consideration.
             */
            eAbund[kk] = test;

            /* *********************************************************** */
            /* **** CHECK LINEAR INDEPENDENCE OF CURRENT FORMULA MATRIX    */
            /* **** LINE WITH PREVIOUS LINES OF THE FORMULA MATRIX  ****** */
            /* *********************************************************** */
            /*
             *          Modified Gram-Schmidt Method, p. 202 Dalquist
             *          QR factorization of a matrix without row pivoting.
             */
            jl = jr;
            /*
             *   Fill in the row for the current element, k, under consideration
             *   The row will contain the Formula matrix value for that element
             *   with respect to the vector of component species.
             *   (note j and k indices are flipped compared to the previous routine)
             */
            for (j = 0; j < nComponents; ++j) {
                jj = orderVectorSpecies[j];
                kk = orderVectorElements[k];
                sm[j + jr*nComponents] = mphase->nAtoms(jj,kk);
            }
            if (jl > 0) {
                /*
                 *         Compute the coefficients of JA column of the
                 *         the upper triangular R matrix, SS(J) = R_J_JR
                 *         (this is slightly different than Dalquist)
                 *         R_JA_JA = 1
                 */
                for (j = 0; j < jl; ++j) {
                    ss[j] = 0.0;
                    for (i = 0; i < nComponents; ++i) {
                        ss[j] += sm[i + jr*nComponents] * sm[i + j*nComponents];
                    }
                    ss[j] /= sa[j];
                }
                /*
                 *     Now make the new column, (*,JR), orthogonal to the
                 *     previous columns
                 */
                for (j = 0; j < jl; ++j) {
                    for (l = 0; l < nComponents; ++l) {
                        sm[l + jr*nComponents] -= ss[j] * sm[l + j*nComponents];
                    }
                }
            }

            /*
             *        Find the new length of the new column in Q.
             *        It will be used in the denominator in future row calcs.
             */
            sa[jr] = 0.0;
            for (ml = 0; ml < nComponents; ++ml) {
                double tmp = sm[ml + jr*nComponents];
                sa[jr] += tmp * tmp;
            }
            /* **************************************************** */
            /* **** IF NORM OF NEW ROW  .LT. 1E-6 REJECT ********** */
            /* **************************************************** */
            if (sa[jr] < 1.0e-6) {
                lindep = true;
            } else {
                lindep = false;
            }
        } while (lindep);
        /* ****************************************** */
        /* **** REARRANGE THE DATA ****************** */
        /* ****************************************** */
        if (jr != k) {
            if (DEBUG_MODE_ENABLED && BasisOptimize_print_lvl > 0) {
                kk = orderVectorElements[k];
                ename = mphase->elementName(kk);
                writelog("   ---   ");
                writelogf("%-2.2s", ename.c_str());
                writelog("replaces ");
                kk = orderVectorElements[jr];
                ename = mphase->elementName(kk);
                writelogf("%-2.2s", ename.c_str());
                writelogf(" as element %3d\n", jr);
            }
            std::swap(orderVectorElements[jr], orderVectorElements[k]);
        }

        /*
         *      If we haven't found enough components, go back
         *      and find some more. (nc -1 is used below, because
         *      jr is counted from 0, via the C convention.
         */
    } while (jr < (nComponents-1));
    return nComponents;
}
