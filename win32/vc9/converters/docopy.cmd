cd ..\..\..\

if not exist build\include\cantera\kernel mkdir build\include\cantera\kernel

cd Cantera\src\conver~1

copy ck2ct.h     ..\..\..\build\include\cantera\kernel
copy CKParser.h     ..\..\..\build\include\cantera\kernel
copy ckr_defs.h     ..\..\..\build\include\cantera\kernel
copy ckr_utils.h     ..\..\..\build\include\cantera\kernel
copy CKReader.h     ..\..\..\build\include\cantera\kernel
copy Constituent.h     ..\..\..\build\include\cantera\kernel
copy Element.h     ..\..\..\build\include\cantera\kernel
copy Reaction.h     ..\..\..\build\include\cantera\kernel
copy RxnSpecies.h     ..\..\..\build\include\cantera\kernel
copy Species.h     ..\..\..\build\include\cantera\kernel
copy thermoFunctions.h     ..\..\..\build\include\cantera\kernel
copy writelog.h     ..\..\..\build\include\cantera\kernel

cd ..\..\..\win32\vc9\conver~1
echo 'ok' 
