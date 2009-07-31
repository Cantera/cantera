# echo on
cd ..\..\..\

if not exist build\include\cantera\kernel mkdir build\include\cantera\kernel

cd Cantera\src\zeroD

copy ConstPressureReactor.h ..\..\..\build\include\cantera\kernel
copy flowControllers.h ..\..\..\build\include\cantera\kernel
copy FlowDevice.h ..\..\..\build\include\cantera\kernel
copy FlowReactor.h ..\..\..\build\include\cantera\kernel
copy PID_Controller.h ..\..\..\build\include\cantera\kernel
copy Reactor.h ..\..\..\build\include\cantera\kernel
copy ReactorBase.h ..\..\..\build\include\cantera\kernel
copy ReactorFactory.h ..\..\..\build\include\cantera\kernel
copy ReactorNet.h ..\..\..\build\include\cantera\kernel
copy Reservoir.h ..\..\..\build\include\cantera\kernel
copy Wall.h ..\..\..\build\include\cantera\kernel

cd ..\..\..\win32\vc9\zeroD
# echo off
echo 'ok' 
