//! @file WallFactory.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/WallFactory.h"
#include "cantera/zeroD/Wall.h"

namespace Cantera
{

WallFactory* WallFactory::s_factory = 0;
std::mutex WallFactory::wall_mutex;

WallFactory::WallFactory()
{
    reg("Wall", []() { return new Wall(); });
}

WallFactory* WallFactory::factory() {
    std::unique_lock<std::mutex> lock(wall_mutex);
    if (!s_factory) {
        s_factory = new WallFactory;
    }
    return s_factory;
}

void WallFactory::deleteFactory() {
    std::unique_lock<std::mutex> lock(wall_mutex);
    delete s_factory;
    s_factory = 0;
}

WallBase* WallFactory::newWall(const string& wallType)
{
    warn_deprecated("WallFactory::newWall",
        "To be removed after Cantera 3.0; for new behavior, see 'newWall3'.");
    return create(wallType);
}

WallBase* newWall(const string& model)
{
    warn_deprecated("newWall",
        "To be changed after Cantera 3.0; for new behavior, see 'newWall3'.");
    return WallFactory::factory()->newWall(model);
}

shared_ptr<WallBase> newWall3(const string& model)
{
    shared_ptr<WallBase> wptr(WallFactory::factory()->create(model));
    return wptr;
}

}
