cd ..\..\..\

if not exist build\include\cantera\kernel mkdir build\include\cantera\kernel

cd Cantera\src\oneD

copy Domain1D.h  ..\..\..\build\include\cantera\kernel
copy Inlet1D.h  ..\..\..\build\include\cantera\kernel
copy MultiJac.h  ..\..\..\build\include\cantera\kernel
copy MultiNewton.h  ..\..\..\build\include\cantera\kernel
copy OneDim.h  ..\..\..\build\include\cantera\kernel
copy refine.h  ..\..\..\build\include\cantera\kernel
copy Resid1D.h  ..\..\..\build\include\cantera\kernel
copy Sim1D.h  ..\..\..\build\include\cantera\kernel
copy StFlow.h  ..\..\..\build\include\cantera\kernel
copy Surf1D.h  ..\..\..\build\include\cantera\kernel

cd ..\..\..\win32\vc9\oneD
echo 'ok' > status
