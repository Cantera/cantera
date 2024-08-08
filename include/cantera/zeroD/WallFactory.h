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
//!     shared_ptr<WallBase> piston = newWall("Wall");
//! ```
class WallFactory : public Factory<WallBase, const string&>
{
public:
    static WallFactory* factory();

    void deleteFactory() override;

private:
    static WallFactory* s_factory;
    static std::mutex wall_mutex;
    WallFactory();
};

//! @defgroup wallGroup Walls
//! Zero-dimensional objects adjacent to reactors.
//! Wall objects should be instantiated via the newWall function, for
//! example:
//!
//! ```cpp
//!     shared_ptr<WallBase> piston = newWall("Wall", "my_piston");
//! ```
//! @ingroup zerodGroup
//! @{

//! Create a WallBase object of the specified type
//! @since Starting in %Cantera 3.1, this method returns a `shared_ptr<WallBase>`
shared_ptr<WallBase> newWall(const string& model, const string& name="(none)");

//! Create a WallBase object of the specified type
//! @since New in %Cantera 3.0.
//! @deprecated Replaced by newWall. To be removed after %Cantera 3.1.
shared_ptr<WallBase> newWall3(const string& model);

//! @}

}

#endif
