/**
 *  @file vcs_util.cpp
 *  Internal definitions for utility functions for the VCSnonideal package
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

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
