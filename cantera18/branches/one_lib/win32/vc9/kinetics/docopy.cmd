echo on
cd ..\..\..\

if not exist build\include\cantera\kernel mkdir build\include\cantera\kernel

cd Cantera\src\kinetics

copy AqueousKinetics.h         ..\..\..\build\include\cantera\kernel
copy EdgeKinetics.h         ..\..\..\build\include\cantera\kernel
copy Enhanced3BConc.h         ..\..\..\build\include\cantera\kernel
copy FalloffFactory.h         ..\..\..\build\include\cantera\kernel
copy FalloffMgr.h         ..\..\..\build\include\cantera\kernel
copy GasKinetics.h         ..\..\..\build\include\cantera\kernel
copy GRI_30_Kinetics.h         ..\..\..\build\include\cantera\kernel
copy Group.h         ..\..\..\build\include\cantera\kernel
copy ImplicitChem.h         ..\..\..\build\include\cantera\kernel
copy ImplicitSurfChem.h         ..\..\..\build\include\cantera\kernel
copy importKinetics.h         ..\..\..\build\include\cantera\kernel
copy InterfaceKinetics.h         ..\..\..\build\include\cantera\kernel
copy Kinetics.h         ..\..\..\build\include\cantera\kernel
copy KineticsFactory.h         ..\..\..\build\include\cantera\kernel
copy RateCoeffMgr.h         ..\..\..\build\include\cantera\kernel
copy reaction_defs.h         ..\..\..\build\include\cantera\kernel
copy ReactionData.h         ..\..\..\build\include\cantera\kernel
copy ReactionPath.h         ..\..\..\build\include\cantera\kernel
copy ReactionStoichMgr.h         ..\..\..\build\include\cantera\kernel
copy RxnRates.h         ..\..\..\build\include\cantera\kernel
copy solveSP.h         ..\..\..\build\include\cantera\kernel
copy StoichManager.h         ..\..\..\build\include\cantera\kernel
copy ThirdBodyMgr.h         ..\..\..\build\include\cantera\kernel

cd ..\..\..\win32\vc9\kinetics

echo off
echo 'ok' 
