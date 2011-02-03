cd ..\..\..\

if not exist build\include\cantera\kernel mkdir build\include\cantera\kernel

cd Cantera\src\numerics

copy ArrayViewer.h        ..\..\..\build\include\cantera\kernel
copy BandMatrix.h        ..\..\..\build\include\cantera\kernel
copy ctlapack.h        ..\..\..\build\include\cantera\kernel
copy CVode.h        ..\..\..\build\include\cantera\kernel  
copy CVodesIntegrator.h        ..\..\..\build\include\cantera\kernel
copy DAE_Solver.h        ..\..\..\build\include\cantera\kernel
copy DASPK.h       ..\..\..\build\include\cantera\kernel
copy DenseMatrix.h        ..\..\..\build\include\cantera\kernel
copy Func1.h        ..\..\..\build\include\cantera\kernel
copy FuncEval.h        ..\..\..\build\include\cantera\kernel
copy funcs.h        ..\..\..\build\include\cantera\kernel
copy IDA_Solver.h        ..\..\..\build\include\cantera\kernel
copy Integrator.h        ..\..\..\build\include\cantera\kernel
copy lapack.h        ..\..\..\build\include\cantera\kernel
copy NonlinearSolver.h        ..\..\..\build\include\cantera\kernel
copy polyfit.h        ..\..\..\build\include\cantera\kernel
copy ResidEval.h        ..\..\..\build\include\cantera\kernel
copy ResidJacEval.h        ..\..\..\build\include\cantera\kernel
copy sort.h        ..\..\..\build\include\cantera\kernel
copy SquareMatrix.h        ..\..\..\build\include\cantera\kernel

cd ..\..\..\win32\vc9\numerics
echo 'ok' > status
