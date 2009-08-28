cd ..\..\..\

if not exist build\include\cantera\kernel mkdir build\include\cantera\kernel

cd Cantera\src\equil

copy ChemEquil.h           ..\..\..\build\include\cantera\kernel
copy equil.h           ..\..\..\build\include\cantera\kernel
copy MultiPhase.h           ..\..\..\build\include\cantera\kernel
copy MultiPhaseEquil.h           ..\..\..\build\include\cantera\kernel
copy PropertyCalculator.h           ..\..\..\build\include\cantera\kernel
copy vcs_defs.h           ..\..\..\build\include\cantera\kernel
copy vcs_DoubleStarStar.h           ..\..\..\build\include\cantera\kernel
copy vcs_Exception.h           ..\..\..\build\include\cantera\kernel
copy vcs_internal.h           ..\..\..\build\include\cantera\kernel
copy vcs_IntStarStar.h           ..\..\..\build\include\cantera\kernel
copy vcs_MultiPhaseEquil.h           ..\..\..\build\include\cantera\kernel
copy vcs_nasa_poly.h           ..\..\..\build\include\cantera\kernel
copy vcs_prob.h           ..\..\..\build\include\cantera\kernel
copy vcs_solve.h           ..\..\..\build\include\cantera\kernel
copy vcs_species_thermo.h           ..\..\..\build\include\cantera\kernel
copy vcs_SpeciesProperties.h           ..\..\..\build\include\cantera\kernel
copy vcs_VolPhase.h           ..\..\..\build\include\cantera\kernel

cd ..\..\..\win32\vc9\equil
echo 'ok' 
