/**
 * @file Electron.h
 * Header file for class Electron.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_ELECTRON_H
#define CT_ELECTRON_H

#include "cantera/electron/ElectronCrossSection.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/ValueCache.h"
#include <boost/math/interpolators/barycentric_rational.hpp>

namespace Cantera
{

class Electron
{
public:
    Electron(); //!< Default constructor.

    virtual ~Electron();

    // Electron objects are not copyable or assignable
    Electron(const Electron&) = delete;
    Electron& operator=(const Electron&) = delete;

    //! Add a electron corss section to this Electron. Returns `true` if the electron cross section was
    //! successfully added, or `false` if the electron cross section was ignored.
    virtual bool addElectronCrossSection(shared_ptr<ElectronCrossSection> ecs);

    size_t nElectronCrossSections() const;

    std::vector<std::string> m_electronCrossSectionTargets;
    std::vector<std::string> m_electronCrossSectionKinds;
    vector_fp m_massRatios;
    std::vector<boost::math::barycentric_rational<double>> m_electronCrossSectionFunctions;

protected:
    //! Cached for saved calculations within each Electron.
    /*!
     *   For more information on how to use this, see examples within the source
     *   code and documentation for this within ValueCache class itself.
     */
    mutable ValueCache m_cache;

    size_t m_ncs;
};

}

#endif
