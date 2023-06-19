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
//! This class is mainly used via the newWall3() function, for example:
//!
//! ```cpp
//!     shared_ptr<WallBase> piston = newWall3("Wall");
//! ```
//!
//! @ingroup ZeroD
class WallFactory : public Factory<WallBase>
{
public:
    static WallFactory* factory();

    virtual void deleteFactory();

    //! Create a new wall by type name.
    /*!
     * @param wallType the type to be created.
     * @deprecated  To be removed after Cantera 3.0; replaceable by newWall3.
     */
    virtual WallBase* newWall(const std::string& wallType);

private:
    static WallFactory* s_factory;
    static std::mutex wall_mutex;
    WallFactory();
};

//! Create a WallBase object of the specified type
//! @deprecated  To be changed after Cantera 3.0; for new behavior, see newWall3.
//! @ingroup ZeroD
WallBase* newWall(const string& model);

//! Create a WallBase object of the specified type
//! @ingroup ZeroD
//! @since New in Cantera 3.0.
//! @todo Transition back to newWall after Cantera 3.0
shared_ptr<WallBase> newWall3(const string& model);

}

#endif
