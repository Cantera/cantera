//! @file ElectronCrossSection.h Declaration for class Cantera::ElectronCrossSection.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_ELECTRONCROSSSECTION_H
#define CT_ELECTRONCROSSSECTION_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

//! Contains data about the cross sections of electron collision
class ElectronCrossSection
{
public:
    ElectronCrossSection();

    //! ElectronCrossSection objects are not copyable or assignable
    ElectronCrossSection(const ElectronCrossSection&) = delete;
    ElectronCrossSection& operator=(const ElectronCrossSection& other) = delete;
    ~ElectronCrossSection();

    //! Validate the cross-section data.
    // void validate();

    //! The name of the kind of electron collision
    string kind;

    //! The name of the target of electron collision
    string target;

    //! The product of electron collision
    string product;

    vector<string> products;

    //! Data of cross section. [m^2]
    vector<double> crossSection;

    //! The energy level corresponding to the cross section. [eV]
    vector<double> energyLevel;

    //! The threshold of a process in [eV]
    double threshold;

    //! Extra data used for specific models
    // AnyMap extra;
};

//! create an ElectronCrossSection object to store data.
unique_ptr<ElectronCrossSection> newElectronCrossSection(const AnyMap& node);

//! Get a vector of ElectronCrossSection objects to access the data.
// std::vector<shared_ptr<ElectronCrossSection>> getElectronCrossSections(const AnyValue& items);

}

#endif
