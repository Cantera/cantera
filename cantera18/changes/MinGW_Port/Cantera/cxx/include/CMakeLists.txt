SET(CXX_H  Cantera.h  equilibrium.h  IncompressibleSolid.h 
        kinetics.h onedim.h  surface.h GRI30.h  integrators.h 
        Metal.h   PureFluid.h   transport.h Edge.h  
        IdealGasMix.h  Interface.h    numerics.h  
        reactionpaths.h  zerodim.h importPhase.h thermo.h 
        radiation.h spectra.h)

INSTALL_FILES(/include/cantera FILES ${CXX_H})
