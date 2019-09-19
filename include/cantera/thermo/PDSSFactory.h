//! @file PDSSFactory.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef PDSS_FACTORY_H
#define PDSS_FACTORY_H

#include "PDSS.h"
#include "cantera/base/FactoryBase.h"
#include "cantera/thermo/VPStandardStateTP.h"

namespace Cantera
{

class PDSSFactory : public Factory<PDSS>
{
public:
    //! Static function that creates a static instance of the factory.
    static PDSSFactory* factory() {
        std::unique_lock<std::mutex> lock(thermo_mutex);
        if (!s_factory) {
            s_factory = new PDSSFactory;
        }
        return s_factory;
    }

    //! delete the static instance of this factory
    virtual void deleteFactory() {
        std::unique_lock<std::mutex> lock(thermo_mutex);
        delete s_factory;
        s_factory = 0;
    }

    //! Create a new thermodynamic property manager.
    /*!
     * @param model  String to look up the model against
     * @returns a pointer to a new PDSS instance matching the model string.
     *   Returns NULL if something went wrong. Throws an exception if the string
     *   wasn't matched.
     */
    virtual PDSS* newPDSS(const std::string& model);

private:
    //! static member of a single instance
    static PDSSFactory* s_factory;

    //! Private constructors prevents usage
    PDSSFactory();

    //! Decl for locking mutex for thermo factory singleton
    static std::mutex thermo_mutex;
};

PDSS* newPDSS(const std::string& model);

}

#endif
