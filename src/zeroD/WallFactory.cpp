//! @file WallFactory.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/WallFactory.h"
#include "cantera/zeroD/Wall.h"

using namespace std;
namespace Cantera
{

WallFactory* WallFactory::s_factory = 0;
std::mutex WallFactory::wall_mutex;

WallFactory::WallFactory()
{
    reg("Wall", []() { return new Wall(); });

    // only used by clib
    reg_type("Wall", WallType);
}

WallBase* WallFactory::newWall(const std::string& wallType)
{
    return create(wallType);
}

WallBase* WallFactory::newWall(int ir)
{
    try {
        return create(m_types.at(ir));
    } catch (out_of_range&) {
        throw CanteraError("WallFactory::newWall",
                           "unknown wall type!");
    }
}

}
