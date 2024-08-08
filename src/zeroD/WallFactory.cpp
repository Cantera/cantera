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
    reg("Wall", [](const string& name) { return new Wall(name); });
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

shared_ptr<WallBase> newWall(const string& model, const string& name)
{
    return shared_ptr<WallBase>(WallFactory::factory()->create(model, name));
}

shared_ptr<WallBase> newWall3(const string& model)
{
    warn_deprecated("newWall3",
        "Use newWall instead; to be removed after Cantera 3.1.");
    return newWall(model);
}

}
