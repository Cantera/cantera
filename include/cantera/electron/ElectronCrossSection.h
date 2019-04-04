//! @file Species.h Declaration for class Cantera::Species.

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_ELECTRONCROSSSECTION_H
#define CT_ELECTRONCROSSSECTION_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

//! Contains data about the cross sections of electron collision
/*!
 *  This class stores the cross section data of electron collision
 */
class ElectronCrossSection
{
public:
    ElectronCrossSection();

    //! Constructor
    ElectronCrossSection(const std::string& kind, const std::string& target,
                         const std::vector<std::vector<double>> data, double mass_ratio=Undef);

    //! ElectronCrossSection objects are not copyable or assignable
    ElectronCrossSection(const ElectronCrossSection&) = delete;
    ElectronCrossSection& operator=(const ElectronCrossSection& other) = delete;
    ~ElectronCrossSection();

    void validate();

    //! The name of the kind of electron collision
    std::string kind;

    //! The name of the target of electron collision
    std::string target;

    //! Data
    std::vector<vector_fp> data;

    //! The mass ratio of molecule to electron
    double mass_ratio;

    //! Extra data used for specific models
    AnyMap extra;
};

//! Create a new ElectronCrossSection object from an AnyMap specification
unique_ptr<ElectronCrossSection> newElectronCrossSection(const AnyMap& node);

//! Generate ElectronCrossSection objects for each item (an AnyMap) in `items`.
std::vector<shared_ptr<ElectronCrossSection>> getElectronCrossSection(const AnyValue& items);

}

#endif
