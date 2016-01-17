/**
 *  @file FactoryBase.h
 *  File contains the FactoryBase class declarations.
 */
#ifndef CT_FACTORY_BASE
#define CT_FACTORY_BASE

#include <vector>

namespace Cantera
{

//! Base class for factories.
/*!   This class maintains a registry of
 *    all factories that derive from it, and deletes them all when
 *    its static method deleteFactories is invoked.
 */
class FactoryBase
{
public:

    //! destructor
    virtual ~FactoryBase() {
    }

    //! static function that deletes all factories
    //! in the internal registry maintained in a static variable
    static void deleteFactories() {
        std::vector< FactoryBase* >::iterator iter ;
        for (iter = s_vFactoryRegistry.begin();
                iter != s_vFactoryRegistry.end();
                ++iter) {
            (*iter)->deleteFactory() ;
        }
        s_vFactoryRegistry.clear() ;
    }

protected:

    //! Constructor.
    /*!
     * Adds the current object to the current static list
     */
    FactoryBase() {
        s_vFactoryRegistry.push_back(this) ;
    }

    //! Virtual abstract function that deletes the factory
    /*!
     *  This must be properly defined in child objects.
     */
    virtual void deleteFactory() = 0 ;

private:

    //! statically held list of Factories.
    static std::vector<FactoryBase*> s_vFactoryRegistry ;
};
}

#endif

