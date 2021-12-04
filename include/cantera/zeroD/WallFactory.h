//! @file WallFactory.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef WALL_FACTORY_H
#define WALL_FACTORY_H

#include "cantera/base/FactoryBase.h"
#include "cantera/zeroD/Wall.h"

namespace Cantera
{

//! Factory class to create WallBase objects
//!
//! This class is mainly used via the newWall() function, for example:
//!
//! ```cpp
//!     unique_ptr<WallBase> piston(newWall("Wall"));
//! ```
//!
//! @ingroup ZeroD
class WallFactory : public Factory<WallBase>
{
public:
    static WallFactory* factory() {
        std::unique_lock<std::mutex> lock(wall_mutex);
        if (!s_factory) {
            s_factory = new WallFactory;
        }
        return s_factory;
    }

    virtual void deleteFactory() {
        std::unique_lock<std::mutex> lock(wall_mutex);
        delete s_factory;
        s_factory = 0;
    }

    //! Create a new wall by type name.
    /*!
     * @param wallType the type to be created.
     */
    virtual WallBase* newWall(const std::string& wallType);

private:
    static WallFactory* s_factory;
    static std::mutex wall_mutex;
    WallFactory();
};

//! Create a WallBase object of the specified type
//! @ingroup ZeroD
inline WallBase* newWall(const std::string& model)
{
    return WallFactory::factory()->newWall(model);
}

}

#endif
