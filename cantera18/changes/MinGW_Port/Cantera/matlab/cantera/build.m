disp('building Phase...');
%buildphase
cd @Phase/private
mex newphase.cpp -lct
mex phase_get.cpp -lct
mex phase_set.cpp -lct
cd ../..
    
disp('building Thermo...');
%buildthermo
cd @Thermo/private
mex newthermo.cpp -lct
mex thermo_get.cpp -lct
mex thermo_set.cpp -lct
cd ../..
    
disp('building Kinetics...');
%buildkinetics
cd @Kinetics/private
mex newkinetics.cpp -lct
mex kin_get.cpp -lct
mex kin_set.cpp -lct
mex rstoich.cpp -lct
mex pstoich.cpp -lct
mex rop.cpp -lct
mex production.cpp -lct
mex isrev.cpp -lct
mex rxnstring.cpp -lct
cd ../..

disp('building Transport...');
%buildtransport
cd @Transport/private
mex newTransport.cpp -lct
mex trans_methods.cpp -lct
cd ../..

%
disp('building Functions...');
cd private
mex addCanteraDirectory.cpp -lct
mex clearStorage.cpp -lct
mex import_from_file.cpp -lct
mex getCanteraError.cpp -lct
cd ..
