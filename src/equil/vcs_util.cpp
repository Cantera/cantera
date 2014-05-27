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
    for (int i = 0; i < length; i++) {
        vector[i] = 0.0;
    }
}
#endif

#ifndef USE_MEMSET
void vcs_izero(int* vector, int length)
{
    for (int i = 0; i < length; i++) {
        vector[i] = 0;
    }
}
#endif

#ifndef USE_MEMSET
void vcs_dcopy(double* const vec_to, const double* const vec_from, int length)
{
    for (int i = 0; i < length; i++) {
        vec_to[i] = vec_from[i];
    }
}
#endif

#ifndef USE_MEMSET
void vcs_icopy(int* vec_to, int* vec_from, int length)
{
    for (int i = 0; i < length; i++) {
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
    size_t largest = j;
    double big = x[j];
    if (xSize) {
        assert(xSize[j] > 0.0);
        big *= xSize[j];
        for (size_t i = j + 1; i < n; ++i) {
            assert(xSize[i] > 0.0);
            if ((x[i] * xSize[i]) > big) {
                largest = i;
                big = x[i] * xSize[i];
            }
        }
    } else {
        for (size_t i = j + 1; i < n; ++i) {
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
    int retn;
    if (vector == NULL || length <= 0) {
        return 0;
    }
    retn = vector[0];
    for (int i = 1; i < length; i++) {
        retn = std::max(retn, vector[i]);
    }
    return retn;
}

double vcsUtil_gasConstant(int mu_units)
{
    switch (mu_units) {
    case VCS_UNITS_KCALMOL:
        return Cantera::GasConst_cal_mol_K * 1e-3;
    case VCS_UNITS_UNITLESS:
        return 1.0;
    case VCS_UNITS_KJMOL:
        return Cantera::GasConstant * 1e-6;
    case VCS_UNITS_KELVIN:
        return 1.0;
    case VCS_UNITS_MKS:
        /* joules / kg-mol K = kg m2 / s2 kg-mol K */
        return Cantera::GasConstant;
    default:
        plogf("vcs_gasConstant error: uknown units: %d\n",
              mu_units);
        exit(EXIT_FAILURE);
    }
}

const char* vcs_speciesType_string(int speciesStatus, int length)
{
    switch (speciesStatus) {
    case VCS_SPECIES_COMPONENT:
        return "Component Species";
    case VCS_SPECIES_MAJOR:
        return "Major Species";
    case VCS_SPECIES_MINOR:
        return "Minor Species";
    case VCS_SPECIES_ZEROEDPHASE:
        if (length < 48) {
            return "Set Zeroed-Phase";
        } else {
            return "Purposely Zeroed-Phase Species (not in problem)";
        }
    case VCS_SPECIES_ZEROEDMS:
        if (length < 23) {
            return "Zeroed-MS Phase";
        } else {
            return "Zeroed-MS Phase Species";
        }
    case VCS_SPECIES_ZEROEDSS:
        if (length < 23) {
            return "Zeroed-SS Phase";
        } else {
            return "Zeroed-SS Phase Species";
        }
    case VCS_SPECIES_DELETED:
        if (length < 22) {
            return "Deleted Species";
        } else if (length < 40) {
            return "Deleted-Small Species";
        } else {
            return "Deleted-Small Species in a MS phase";
        }
    case VCS_SPECIES_ACTIVEBUTZERO:
        if (length < 47) {
            return "Tmp Zeroed in MS";
        } else {
            return "Zeroed Species in an active MS phase (tmp)";
        }
    case VCS_SPECIES_STOICHZERO:
        if (length < 56) {
            return "Stoich Zeroed in MS";
        } else {
            return "Zeroed Species in an active MS phase (Stoich Constraint)";
        }
    case VCS_SPECIES_INTERFACIALVOLTAGE:
        if (length < 29) {
            return "InterfaceVoltage";
        } else {
            return "InterfaceVoltage Species";
        }
    default:
        return "unknown species type";
    }
}

void vcs_print_stringTrunc(const char* str, size_t space, int alignment)
{
    size_t ls = 0, rs = 0;
    size_t len = strlen(str);
    if ((len) >= space) {
        for (size_t i = 0; i < space; i++) {
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
            for (size_t i = 0; i < ls; i++) {
                plogf(" ");
            }
        }
        plogf("%s", str);
        if (rs != 0) {
            for (size_t i = 0; i < rs; i++) {
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
