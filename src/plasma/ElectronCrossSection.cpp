// This file is part of Cantera. See License.txt in the top-level directory or
// at https://www.cantera.org/license.txt for license and copyright information.

#include "cantera/plasma/ElectronCrossSection.h"

namespace Cantera {

ElectronCrossSection::ElectronCrossSection()
    : threshold(0.0)
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
    } else if (kind == "ELASTIC") {
        if (data[0][0] != 0.0) {
            throw CanteraError("ElectronCrossSection::validate",
                "Invalid energy value of type '{}' for '{}'. "
                "Energy must starts at zero.", kind, target);
        }
    } else if (kind != "IONIZATION" && kind != "ATTACHMENT" && kind != "EXCITATION"){
        throw CanteraError("ElectronCrossSection::validate",
            "'{}' is an unknown type of cross section data.", kind);
    }
}

}
