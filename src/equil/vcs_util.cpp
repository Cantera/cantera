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
#include "cantera/equil/vcs_defs.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"
#include <cassert>
#include <cstring>

using namespace std;

namespace Cantera
{

double vcs_l2norm(const vector_fp& vec)
{
    if (vec.empty()) {
        return 0.0;
    }
    double sum = 0.0;
    for (const auto& val : vec) {
        sum += val * val;
    }
    return std::sqrt(sum / vec.size());
}

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

}
