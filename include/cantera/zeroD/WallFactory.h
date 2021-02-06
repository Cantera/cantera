//! @file WallFactory.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef WALL_FACTORY_H
#define WALL_FACTORY_H

#include "cantera/base/FactoryBase.h"
#include "cantera/zeroD/Wall.h"

namespace Cantera
{

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

    //! Create a new wall by type identifier.
    /*!
     * @param n the type to be created.
     */
    virtual WallBase* newWall(int n);

    //! Create a new wall by type name.
    /*!
     * @param wallType the type to be created.
     */
    virtual WallBase* newWall(const std::string& wallType);

    //! Register a new wall type identifier.
    /*!
     * @param name the name of the wall type.
     * @param type the type identifier of the wall.
     * Integer type identifiers are used by clib and matlab interfaces.
     *
     * @deprecated To be removed after Cantera 2.5.
     */
    void reg_type(const std::string& name, const int type) {
        m_types[type] = name;
    }

protected:
    //! Map containing wall type identifier / wall type name pairs.
    //! @deprecated To be removed after Cantera 2.5.
    std::unordered_map<int, std::string> m_types;

private:
    static WallFactory* s_factory;
    static std::mutex wall_mutex;
    WallFactory();
};

//! Create a Wall object of the specified type
inline WallBase* newWall(const std::string& model)
{
    return WallFactory::factory()->newWall(model);
}

}

#endif
