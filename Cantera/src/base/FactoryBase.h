#ifndef CT_FACTORY_BASE
#define CT_FACTORY_BASE

#include <vector>

#if defined(THREAD_SAFE_CANTERA)
#include <boost/thread/mutex.hpp>
#endif

namespace Cantera {

    /**
     * Base class for factories. This class maintains a registry of
     * all factories that derive from it, and deletes them all when
     * its static method deleteFactories is invoked.
     */
    class FactoryBase
    {
    public:
        virtual ~FactoryBase()
            {
            }
        
        static void deleteFactories()
            {
                std::vector< FactoryBase* >::iterator iter ;
                for ( iter = s_vFactoryRegistry.begin(); 
                      iter != s_vFactoryRegistry.end(); 
                      ++iter)
                {
                    (*iter)->deleteFactory() ;
                }
                s_vFactoryRegistry.clear() ;
            }
    
    protected:
        FactoryBase()
            {
                s_vFactoryRegistry.push_back(this) ;
            }

        virtual void deleteFactory() = 0 ;
        
    private:
        static std::vector<FactoryBase*> s_vFactoryRegistry ;
    };
}

#endif

