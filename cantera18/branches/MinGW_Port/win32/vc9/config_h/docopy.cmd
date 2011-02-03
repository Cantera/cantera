echo on
cd ..\..\..
copy winconfig.h config.h

cd ext\f2c_libs
copy arith.hwin32 arith.h
copy sysdep1.h0 sysdep1.h
copy signal1.h0 signal1.h
cd ..\..

if not exist build\include\cantera mkdir build\include\cantera
if not exist build\include\cantera\kernel mkdir build\include\cantera\kernel
if not exist build\bin\i686-pc-win32 mkdir build\bin\i686-pc-win32
if not exist build\lib\i686-pc-win32 mkdir build\lib\i686-pc-win32

copy config.h build\include\cantera
copy config.h build\include\cantera\winconfig.h

cd Cantera\cxx\include

copy Cantera.h              ..\..\..\build\include\cantera
copy Edge.h                 ..\..\..\build\include\cantera
copy electrolyteThermo.h    ..\..\..\build\include\cantera
copy equilibrium.h          ..\..\..\build\include\cantera
copy GRI30.h                ..\..\..\build\include\cantera
copy IdealGasMix.h          ..\..\..\build\include\cantera
copy importPhase.h          ..\..\..\build\include\cantera
copy IncompressibleSolid.h  ..\..\..\build\include\cantera
copy integrators.h          ..\..\..\build\include\cantera
copy Interface.h            ..\..\..\build\include\cantera
copy kinetics.h             ..\..\..\build\include\cantera
copy Metal.h                ..\..\..\build\include\cantera
copy numerics.h             ..\..\..\build\include\cantera
copy onedim.h               ..\..\..\build\include\cantera
copy PureFluid.h            ..\..\..\build\include\cantera
copy radiation.h            ..\..\..\build\include\cantera
copy reactionpaths.h        ..\..\..\build\include\cantera
copy spectra.h              ..\..\..\build\include\cantera
copy surface.h              ..\..\..\build\include\cantera
copy thermo.h               ..\..\..\build\include\cantera
copy transport.h            ..\..\..\build\include\cantera
copy zerodim.h              ..\..\..\build\include\cantera

cd ..\..\..

cd win32\vc9\config_h
echo ok
echo off
