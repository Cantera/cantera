// This file is part of Cantera. See License.txt in the top-level directory or
// at https://www.cantera.org/license.txt for license and copyright information.

#include "cantera/electron/ElectronCrossSection.h"

namespace Cantera {

ElectronCrossSection::ElectronCrossSection()
    : mass_ratio(Undef)
    , threshold(0.0)
{
}

ElectronCrossSection::~ElectronCrossSection()
{
}

void ElectronCrossSection::validate()
{
    if (kind == "EFFECTIVE") {
        if (data.size() > 0 && data[0].size() > 0) {
            if (data[0][0] != 0.0) {
                throw CanteraError("ElectronCrossSection::validate",
                    "Invalid energy value of type '{}' for '{}'. "
                    "Energy must starts at zero.", kind, target);
            }
        }
        if (mass_ratio >= 1.0 || mass_ratio < 0.0) {
            throw CanteraError("ElectronCrossSection::validate",
                "Invalid mass ratio of type '{}' for '{}'. "
                "Mass ratio of electron to target must be in the range of 0 to 1.", kind, target);
        }
    } else if (kind == "ELASTIC") {
        if (data[0][0] != 0.0) {
            throw CanteraError("ElectronCrossSection::validate",
                "Invalid energy value of type '{}' for '{}'. "
                "Energy must starts at zero.", kind, target);
        }
        if (mass_ratio >= 1.0 || mass_ratio < 0.0) {
            throw CanteraError("ElectronCrossSection::validate",
                "Invalid mass ratio of type '{}' for '{}'. "
                "Mass ratio of electron to target must be in the range of 0 to 1.", kind, target);
        }
    } else if (kind != "IONIZATION" && kind != "ATTACHMENT" && kind != "EXCITATION"){
        throw CanteraError("ElectronCrossSection::validate",
            "'{}' is an unknown type of cross section data.", kind);
    }
}

}
