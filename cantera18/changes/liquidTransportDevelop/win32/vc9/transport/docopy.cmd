cd ..\..\..\

if not exist build\include\cantera\kernel mkdir build\include\cantera\kernel

cd Cantera\src\transp~1

copy AqueousTransport.h         ..\..\..\build\include\cantera\kernel
copy DustyGasTransport.h    ..\..\..\build\include\cantera\kernel
copy FtnTransport.h    ..\..\..\build\include\cantera\kernel
copy L_matrix.h    ..\..\..\build\include\cantera\kernel
copy LiquidTransport.h    ..\..\..\build\include\cantera\kernel
copy LiquidTransportParams.h    ..\..\..\build\include\cantera\kernel
copy MixTransport.h    ..\..\..\build\include\cantera\kernel
copy MMCollisionInt.h    ..\..\..\build\include\cantera\kernel
copy MultiTransport.h    ..\..\..\build\include\cantera\kernel
copy SolidTransport.h    ..\..\..\build\include\cantera\kernel
copy TransportBase.h    ..\..\..\build\include\cantera\kernel
copy TransportFactory.h    ..\..\..\build\include\cantera\kernel
copy TransportParams.h    ..\..\..\build\include\cantera\kernel

cd ..\..\..\win32\vc9\transp~1
echo 'ok' > status
