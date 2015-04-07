/**
 *  @file vcs_util.cpp
 *  Internal definitions for utility functions for the VCSnonideal package
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/equil/vcs_internal.h"
#include <cassert>

using namespace std;

namespace VCSnonideal
{
#ifndef USE_MEMSET
void vcs_dzero(double* vector, int length)
{
    int i;
    for (i = 0; i < length; i++) {
        vector[i] = 0.0;
    }
}
#endif

#ifndef USE_MEMSET
void vcs_izero(int* vector, int length)
{
    int i;
    for (i = 0; i < length; i++) {
        vector[i] = 0;
    }
}
#endif

#ifndef USE_MEMSET
void vcs_dcopy(double* const vec_to, const double* const vec_from, int length)
{
    int i;
    for (i = 0; i < length; i++) {
        vec_to[i] = vec_from[i];
    }
}
#endif

#ifndef USE_MEMSET
void vcs_icopy(int* vec_to, int* vec_from, int length)
{
    int i;
    for (i = 0; i < length; i++) {
        vec_to[i] = vec_from[i];
    }
}
#endif

#ifndef USE_MEMSET
void vcs_vdzero(std::vector<double> &vvv, int len)
{
    if (len < 0) {
        std::fill(vvv.begin(), vvv.end(), 0.0);
    } else {
        std::fill_n(vvv.begin(), len, 0.0);
    }
}
#endif

double vcs_l2norm(const std::vector<double> vec)
{
    size_t len = vec.size();
    if (len == 0) {
        return 0.0;
    }
    double sum = 0.0;
    std::vector<double>::const_iterator pos;
    for (pos = vec.begin(); pos != vec.end(); ++pos) {
        sum += (*pos) * (*pos);
    }
    return std::sqrt(sum / len);
}

#ifndef USE_MEMSET
void vcs_vizero(std::vector<int> &vvv, int len)
{
    if (len < 0) {
        std::fill(vvv.begin(), vvv.end(), 0.0);
    } else {
        std::fill_n(vvv.begin(), len, 0.0);
    }
}
#endif

#ifndef USE_MEMSET
void vcs_vdcopy(std::vector<double> &vec_to,
                const std::vector<double> & vec_from, int length)
{
    std::copy(vec_from.begin(), vec_from.begin() + length, vec_to.begin());
}
#endif

#ifndef USE_MEMSET
void vcs_vicopy(std::vector<int> &vec_to,
                const std::vector<int> & vec_from, int length)
{
    std::copy(vec_from.begin(), vec_from.begin() + length, vec_to.begin());
}
#endif

size_t vcs_optMax(const double* x, const double* xSize, size_t j, size_t n)
{
    size_t i;
    size_t largest = j;
    double big = x[j];
    if (xSize) {
        assert(xSize[j] > 0.0);
        big *= xSize[j];
        for (i = j + 1; i < n; ++i) {
            assert(xSize[i] > 0.0);
            if ((x[i] * xSize[i]) > big) {
                largest = i;
                big = x[i] * xSize[i];
            }
        }
    } else {
        for (i = j + 1; i < n; ++i) {
            if (x[i] > big) {
                largest = i;
                big = x[i];
            }
        }
    }
    return largest;
}

int vcs_max_int(const int* vector, int length)
{
    int i, retn;
    if (vector == NULL || length <= 0) {
        return 0;
    }
    retn = vector[0];
    for (i = 1; i < length; i++) {
        retn = std::max(retn, vector[i]);
    }
    return retn;
}

#ifdef DEBUG_HKM
static void mlequ_matrixDump(double* c, int idem, int n)
{
    int i, j;
    printf("vcsUtil_mlequ()     MATRIX DUMP --------------------------------------------------\n");
    printf("      ");
    for (j = 0; j < n; ++j) {
        printf("     % 3d   ", j);
    }
    printf("\n");
    for (j = 0; j < n; ++j) {
        printf("-----------");
    }
    printf("\n");
    for (i = 0; i < n; ++i) {
        printf(" %3d | ", i);
        for (j = 0; j < n; ++j) {
            printf("% 10.3e ", c[i + j * idem]);
        }
        printf("\n");
    }
    for (j = 0; j < n; ++j) {
        printf("-----------");
    }
    printf("\n");
    printf("vcsUtil_mlequ() END MATRIX DUMP --------------------------------------------------\n");

}
#endif

//!  Swap rows in the c matrix and the b rhs matrix
/*!
 *  @param c          Matrix of size nxn, row first
 *  @param idem       C storage dimension for the number of rows
 *  @param n          Size of the matrix
 *  @param b          RHS of the Ax=b problem to solve
 *  @param m          Number of rhs to solve
 *  @param irowa      first row to swap
 *  @param irowb      second row to swap
 */
static void vcsUtil_swapRows(double* c, size_t idem, size_t n, double* b,
                             size_t m, size_t irowa, size_t irowb)
{
    if (irowa == irowb) {
        return;
    }
    for (size_t j = 0; j < n; j++) {
        std::swap(c[irowa + j * idem], c[irowb + j * idem]);
    }
    for (size_t j = 0; j < m; j++) {
        std::swap(b[irowa + j * idem], b[irowb + j * idem]);
    }
}

//!  Swap rows in the c matrix and the b rhs matrix to lower the condition number of the matrix
/*!
 *  @param c          Matrix of size nxn, row first
 *  @param idem       C storage dimension for the number of rows
 *  @param n          Size of the matrix
 *  @param b          RHS of the Ax=b problem to solve
 *  @param m          Number of rhs to solve
 */
static void vcsUtil_mlequ_preprocess(double* c, size_t idem, size_t n,
                                     double* b, size_t m)
{
    size_t j = 0;
    std::vector<int> irowUsed(n, 0);

    for (j = 0; j < n; j++) {
        int numNonzero = 0;
        size_t inonzero = npos;
        for (size_t i = 0; i < n; i++) {
            if (c[i + j * idem] != 0.0) {
                numNonzero++;
                inonzero = i;
            }
        }
        if (numNonzero == 1) {
            if (inonzero != j) {
                if (irowUsed[inonzero] == 0) {
                    vcsUtil_swapRows(c, idem, n, b, m, j, inonzero);
#ifdef DEBUG_HKM
                    // mlequ_matrixDump(c, idem, n);
#endif
                }
            }
            irowUsed[j] = 1;
        }
    }

    for (j = 0; j < n; j++) {
        if (c[j + j * idem] == 0.0) {
            int numNonzero = 0;
            size_t inonzero = npos;
            for (size_t i = 0; i < n; i++) {
                if (!irowUsed[i]) {
                    if (c[i + j * idem] != 0.0) {
                        if ((c[i + i * idem] == 0.0)
                                || (c[j + i * idem] != 0.0)) {
                            numNonzero++;
                            inonzero = i;
                        }
                    }
                }
            }
            if (numNonzero == 1) {
                if (inonzero != j) {
                    if (irowUsed[inonzero] == 0) {
                        vcsUtil_swapRows(c, idem, n, b, m, j, inonzero);
#ifdef DEBUG_HKM
                        // mlequ_matrixDump(c, idem, n);
#endif
                    }
                }
                irowUsed[j] = 1;
            }
        }
    }

    for (j = 0; j < n; j++) {
        if (c[j + j * idem] == 0.0) {
            int numNonzero = 0;
            size_t inonzero = npos;
            for (size_t i = 0; i < n; i++) {
                if (!irowUsed[i]) {
                    if (c[i + j * idem] != 0.0) {
                        if ((c[i + i * idem] == 0.0)
                                || (c[j + i * idem] != 0.0)) {
                            numNonzero++;
                            inonzero = i;
                        }
                    }
                }
            }
            if (inonzero != npos) {
                if (inonzero != j) {
                    if (irowUsed[inonzero] == 0) {
                        vcsUtil_swapRows(c, idem, n, b, m, j, inonzero);
#ifdef DEBUG_HKM
                        // mlequ_matrixDump(c, idem, n);
#endif
                    }
                }
            }
        }
    }
}

int vcsUtil_mlequ(double* c, size_t idem, size_t n, double* b, size_t m)
{
    size_t k;
#ifdef DEBUG_HKM
    // mlequ_matrixDump(c, idem, n);
#endif
    vcsUtil_mlequ_preprocess(c, idem, n, b, m);
#ifdef DEBUG_HKM
    // mlequ_matrixDump(c, idem, n);
    static int s_numCalls = 0;
    s_numCalls++;
#endif

    double R;
    if (n > idem || n <= 0) {
        plogf("vcsUtil_mlequ ERROR: badly dimensioned matrix: %d %d\n", n, idem);
        return 1;
    }

#ifdef DEBUG_HKM
    int dmatrix = 0;
    for (size_t i = 0; i < n; ++i) {
        bool notFound = true;
        for (size_t j = 0; j < n; ++j) {
            if (c[i + j * idem] != 0.0) {
                notFound = false;
            }
        }
        if (notFound) {
            printf(" vcsUtil_mlequ ERROR(): row %d is identically zero\n", i);
        }
    }
    for (size_t j = 0; j < n; ++j) {
        bool notFound = true;
        for (size_t i = 0; i < n; ++i) {
            if (c[i + j * idem] != 0.0) {
                notFound = false;
            }
        }
        if (notFound) {
            printf(" vcsUtil_mlequ ERROR(): column %d is identically zero\n", j);
        }
    }
    //  if (s_numCalls >= 32) {
    // printf("vcsUtil_mlequ: we are here\n");
    //  dmatrix = 1;
    // }

    if (dmatrix) {
        mlequ_matrixDump(c, idem, n);
    }
#endif
    /*
     * Loop over the rows
     *    -> At the end of each loop, the only nonzero entry in the column
     *       will be on the diagonal. We can therfore just invert the
     *       diagonal at the end of the program to solve the equation system.
     */
    for (size_t i = 0; i < n; ++i) {
        if (c[i + i * idem] == 0.0) {
            /*
             *   Do a simple form of row pivoting to find a non-zero pivot
             */
            for (k = i + 1; k < n; ++k) {
                if (c[k + i * idem] != 0.0) {
                    goto FOUND_PIVOT;
                }
            }
            plogf("vcsUtil_mlequ ERROR: Encountered a zero column: %d\n", i);
#ifdef DEBUG_HKM
            plogf("                     call # %d\n", s_numCalls);
#endif
#ifdef DEBUG_HKM
            mlequ_matrixDump(c, idem, n);
#endif
            return 1;
FOUND_PIVOT:
            ;
            for (size_t j = 0; j < n; ++j) {
                c[i + j * idem] += c[k + j * idem];
            }
            for (size_t j = 0; j < m; ++j) {
                b[i + j * idem] += b[k + j * idem];
            }
        }

        for (size_t l = 0; l < n; ++l) {
            if (l != i && c[l + i * idem] != 0.0) {
                R = c[l + i * idem] / c[i + i * idem];
                c[l + i * idem] = 0.0;
                for (size_t j = i + 1; j < n; ++j) {
                    c[l + j * idem] -= c[i + j * idem] * R;
                }
                for (size_t j = 0; j < m; ++j) {
                    b[l + j * idem] -= b[i + j * idem] * R;
                }
            }
        }
    }
    /*
     *  The negative in the last expression is due to the form of B upon
     *  input
     */
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            b[i + j * idem] = -b[i + j * idem] / c[i + i * idem];
        }
    }
    return 0;
}

int vcsUtil_gaussj(double* c, size_t idem, size_t n, double* b, size_t m)
{
    size_t i, j, k, l, ll;
    size_t irow = npos;
    size_t icol = npos;
    bool needInverse = false;
    double pivinv;
#ifdef DEBUG_HKM
    static int s_numCalls = 0;
    s_numCalls++;
#endif
#ifdef DEBUG_HKM
    // mlequ_matrixDump(c, idem, n);
#endif
    /*
     *  Preprocess the problem
     */
    vcsUtil_mlequ_preprocess(c, idem, n, b, m);

#ifdef DEBUG_HKM
    // mlequ_matrixDump(c, idem, n);
#endif

    std::vector<size_t> indxc(n);
    std::vector<size_t> indxr(n);
    std::vector<int> ipiv(n, 0);
    doublereal big = 0.0;
    /*
     *  This is the main loop over the columns to be reduced.
     */
    for (i = 0; i < n; i++) {
        big = 0.0;
        for (j = 0; j < n; j++) {
            if (ipiv[j] != 1) {
                for (k = 0; k < n; k++) {
                    if (ipiv[k] == 0) {
                        if (fabs(c[j + idem * k]) >= big) {
                            big = fabs(c[j + idem * k]);
                            irow = j;
                            icol = k;
                        }
                    }
                }
            }
        }
        ++(ipiv[icol]);
        if (irow != icol) {
            vcsUtil_swapRows(c, idem, n, b, m, irow, icol);
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (c[icol + idem * icol] == 0.0) {
            plogf("vcsUtil_gaussj ERROR: Encountered a zero column: %d\n", i);
            return 1;
        }
        pivinv = 1.0 / c[icol + idem * icol];
        c[icol + idem * icol] = 1.0;
        for (l = 0; l < n; l++) {
            c[icol + idem * l] *= pivinv;
        }
        for (l = 0; l < m; l++) {
            b[icol + idem * l] *= pivinv;
        }
        for (ll = 0; ll < n; ll++) {
            if (ll != icol) {
                double dum = c[ll + idem * icol];
                c[ll + idem * icol] = 0;
                for (l = 0; l < n; l++) {
                    c[ll + idem * l] -= c[icol + idem * l] * dum;
                }
                for (l = 0; l < m; l++) {
                    b[ll + idem * l] -= b[icol + idem * l] * dum;
                }
            }
        }
    }
    if (needInverse) {
        for (l = n - 1; l != npos; l--) {
            if (indxr[l] != indxc[l]) {
                for (k = 0; k < n; k++) {
                    std::swap(c[k + idem * indxr[l]], c[k + idem * indxr[l]]);
                }
            }
        }
    }

    /*
     *  The negative in the last expression is due to the form of B upon
     *  input
     */
    for (i = 0; i < n; ++i) {
        for (j = 0; j < m; ++j) {
            b[i + j * idem] = -b[i + j * idem];
        }
    }
    return 0;
}

double vcsUtil_gasConstant(int mu_units)
{
    double r;
    switch (mu_units) {
    case VCS_UNITS_KCALMOL:
        r = Cantera::GasConst_cal_mol_K * 1e-3;
        break;
    case VCS_UNITS_UNITLESS:
        r = 1.0;
        break;
    case VCS_UNITS_KJMOL:
        r = Cantera::GasConstant * 1e-6;
        break;
    case VCS_UNITS_KELVIN:
        r = 1.0;
        break;
    case VCS_UNITS_MKS:
        /* joules / kg-mol K = kg m2 / s2 kg-mol K */
        r = Cantera::GasConstant;
        break;
    default:
        plogf("vcs_gasConstant error: uknown units: %d\n",
              mu_units);
        exit(EXIT_FAILURE);
    }
    return r;
}

void vcs_print_line(const char* string, int num)
{
    if (string) {
        for (int j = 0; j < num; j++) {
            plogf("%s", string);
        }
    }
    plogendl();
}

const char* vcs_speciesType_string(int speciesStatus, int length)
{
    const char* sss;
    switch (speciesStatus) {
    case VCS_SPECIES_COMPONENT:
        sss = "Component Species";
        break;
    case VCS_SPECIES_MAJOR:
        sss = "Major Species";
        break;
    case VCS_SPECIES_MINOR:
        sss = "Minor Species";
        break;
    case VCS_SPECIES_ZEROEDPHASE:
        if (length < 48) {
            sss = "Set Zeroed-Phase";
        } else {
            sss = "Purposely Zeroed-Phase Species (not in problem)";
        }
        break;
    case VCS_SPECIES_ZEROEDMS:
        if (length < 23) {
            sss = "Zeroed-MS Phase";
        } else {
            sss = "Zeroed-MS Phase Species";
        }
        break;
    case VCS_SPECIES_ZEROEDSS:
        if (length < 23) {
            sss = "Zeroed-SS Phase";
        } else {
            sss = "Zeroed-SS Phase Species";
        }
        break;
    case VCS_SPECIES_DELETED:
        if (length < 22) {
            sss = "Deleted Species";
        } else if (length < 40) {
            sss = "Deleted-Small Species";
        } else {
            sss = "Deleted-Small Species in a MS phase";
        }
        break;
    case VCS_SPECIES_ACTIVEBUTZERO:
        if (length < 47) {
            sss = "Tmp Zeroed in MS";
        } else {
            sss = "Zeroed Species in an active MS phase (tmp)";
        }
        break;
    case VCS_SPECIES_STOICHZERO:
        if (length < 56) {
            sss = "Stoich Zeroed in MS";
        } else {
            sss = "Zeroed Species in an active MS phase (Stoich Constraint)";
        }
        break;
    case VCS_SPECIES_INTERFACIALVOLTAGE:
        if (length < 29) {
            sss = "InterfaceVoltage";
        } else {
            sss = "InterfaceVoltage Species";
        }
        break;
    default:
        sss = "unknown species type";
    }
    return sss;
}

void vcs_print_stringTrunc(const char* str, size_t space, int alignment)
{
    size_t i, ls = 0, rs = 0;
    size_t len = strlen(str);
    if ((len) >= space) {
        for (i = 0; i < space; i++) {
            plogf("%c", str[i]);
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
                plogf(" ");
            }
        }
        plogf("%s", str);
        if (rs != 0) {
            for (i = 0; i < rs; i++) {
                plogf(" ");
            }
        }
    }
}

bool vcs_doubleEqual(double d1, double d2)
{
    double denom = fabs(d1) + fabs(d2) + 1.0;
    double fac = fabs(d1 - d2) / denom;
    if (fac > 1.0E-10) {
        return false;
    }
    return true;
}

void vcs_heapsort(std::vector<int> & x)
{
    int n = x.size();
    if (n < 2) {
        return;
    }
    doublereal rra;
    int ll = n / 2;
    int iret = n - 1;

    while (1 > 0) {
        if (ll > 0) {
            ll--;
            rra = x[ll];
        } else {
            rra = x[iret];
            x[iret] = x[0];
            iret--;
            if (iret == 0) {
                x[0] = rra;
                return;
            }
        }

        int i = ll;
        int j = ll + ll + 1;

        while (j <= iret) {
            if (j < iret) {
                if (x[j] < x[j + 1]) {
                    j++;
                }
            }
            if (rra < x[j]) {
                x[i] = x[j];
                i = j;
                j = j + j + 1;
            } else {
                j = iret + 1;
            }
        }
        x[i] = rra;
    }
}

void vcs_orderedUnique(std::vector<int> & xOrderedUnique, const std::vector<int> & x)
{
    std::vector<int> xordered(x);
    vcs_heapsort(xordered);
    int lastV = x[0] - 1;
    xOrderedUnique.clear();
    for (int i = 0; i < (int) xordered.size(); i++) {
        if (lastV != xordered[i]) {
            xOrderedUnique.push_back(xordered[i]);
            lastV = xordered[i];
        }
    }
}

}
