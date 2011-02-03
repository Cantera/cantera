cd ..\..\..\Cantera\cxx\include

copy   Cantera.h        ..\..\..\build\include\cantera
copy   Edge.h       ..\..\..\build\include\cantera
copy   equilibrium.h      ..\..\..\build\include\cantera
copy   GRI30.h      ..\..\..\build\include\cantera
copy   IdealGasMix.h      ..\..\..\build\include\cantera
copy   IncompressibleSolid.h      ..\..\..\build\include\cantera
copy   integrators.h      ..\..\..\build\include\cantera
copy   Interface.h      ..\..\..\build\include\cantera
copy   kinetics.h      ..\..\..\build\include\cantera
copy   Metal.h      ..\..\..\build\include\cantera
copy   numerics.h      ..\..\..\build\include\cantera
copy   onedim.h      ..\..\..\build\include\cantera
copy   reactionpaths.h      ..\..\..\build\include\cantera
copy   surface.h      ..\..\..\build\include\cantera
copy   transport.h      ..\..\..\build\include\cantera
copy   zerodim.h      ..\..\..\build\include\cantera

cd ..\..\..\win32\vc9\cxxutils
echo 'ok' > status
