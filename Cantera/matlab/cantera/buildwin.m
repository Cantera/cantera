disp('thermodynamic properties...');
%buildthermo
cd @ThermoPhase/private
mex newthermo.cpp ../../../../../lib/ct13r5.lib
mex phase_get.cpp ../../../../../lib/ct13r5.lib
mex phase_set.cpp ../../../../../lib/ct13r5.lib
mex thermo_get.cpp ../../../../../lib/ct13r5.lib
mex thermo_set.cpp ../../../../../lib/ct13r5.lib
cd ../..
    
disp('chemical kinetics...');
%buildkinetics
cd @Kinetics/private
mex newkinetics.cpp ../../../../../lib/ct13r5.lib
mex delkinetics.cpp ../../../../../lib/ct13r5.lib
mex kin_get.cpp ../../../../../lib/ct13r5.lib
mex kin_set.cpp ../../../../../lib/ct13r5.lib
mex rstoich.cpp ../../../../../lib/ct13r5.lib
mex pstoich.cpp ../../../../../lib/ct13r5.lib
mex rop.cpp ../../../../../lib/ct13r5.lib
mex production.cpp ../../../../../lib/ct13r5.lib
mex isrev.cpp ../../../../../lib/ct13r5.lib
mex rxnstring.cpp ../../../../../lib/ct13r5.lib
cd ../..

disp('transport properties...');
%buildtransport
cd @Transport/private
mex newTransport.cpp ../../../../../lib/ct13r5.lib
mex trans_methods.cpp ../../../../../lib/ct13r5.lib
cd ../..

%
disp('XML support...');
cd @XML_Node/private
mex newxml.cpp ../../../../../lib/ct13r5.lib
mex xmlmethods.cpp ../../../../../lib/ct13r5.lib
cd ../..

%
disp('utility functions...');
cd private
mex addCanteraDirectory.cpp ../../../../lib/ct13r5.lib
mex clearStorage.cpp ../../../../lib/ct13r5.lib
%mex import_from_file.cpp ../../../../lib/ct13r5.lib
mex getCanteraError.cpp ../../../../lib/ct13r5.lib
mex ck_to_ctml.cpp ../../../../lib/ct13r5.lib
cd ..

disp('Zero-dimensional reactors...');
cd @Reactor/private
mex reactormethods.cpp ../../../../../lib/ct13r5.lib
cd ../..
cd @Wall/private
mex wallmethods.cpp ../../../../../lib/ct13r5.lib
cd ../..
cd @FlowDevice/private
mex flowdevicemethods.cpp ../../../../../lib/ct13r5.lib
cd ../..    
