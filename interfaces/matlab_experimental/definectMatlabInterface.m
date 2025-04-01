%% About definectMatlabInterface.m
% This file defines the MATLAB interface to the library |ctMatlabInterface|.
%
% Commented sections represent C++ functionality that MATLAB cannot automatically define. To include
% functionality, uncomment a section and provide values for <SHAPE>, <DIRECTION>, etc. For more
% information, see helpview(fullfile(docroot,'matlab','helptargets.map'),'cpp_define_interface') to "Define MATLAB Interface for C++ Library".



%% Setup
% Do not edit this setup section.
function libDef = definectMatlabInterface()
libDef = clibgen.LibraryDefinition("ctMatlabInterfaceData.xml");

%% OutputFolder and Libraries 
libDef.OutputFolder = "C:\Users\ssnit\cantera\interfaces\matlab_experimental";
libDef.Libraries = [ "C:\Users\ssnit\cantera\build\lib\cantera_shared.dll" "C:\Users\ssnit\cantera\build\lib\cantera_shared.lib" ];

%% C++ function type |using LogCallback = void(*)(LogLevel, char const *, char const *)| with MATLAB name |clib.ctMatlabInterface.LogCallback| 
addFunctionType(libDef, "using LogCallback = void(*)(LogLevel, char const *, char const *)", "MATLABName", "clib.ctMatlabInterface.LogCallback", ...
    "Description", "clib.ctMatlabInterface.LogCallback    C++ function type LogCallback."); % Modify help description values as needed.

%% C++ enumeration |LogLevel| with MATLAB name |clib.ctMatlabInterface.LogLevel| 
addEnumeration(libDef, "LogLevel", "int32",...
    [...
      "INFO",...  % 0
      "WARN",...  % 1
      "ERROR",...  % 2
    ],...
    "MATLABName", "clib.ctMatlabInterface.LogLevel", ...
    "Description", "clib.ctMatlabInterface.LogLevel    Representation of C++ enumeration LogLevel."); % Modify help description values as needed.

%% C++ function |ct_appdelete| with MATLAB name |clib.ctMatlabInterface.ct_appdelete|
% C++ Signature: int ct_appdelete()

ct_appdeleteDefinition = addFunction(libDef, ...
    "int ct_appdelete()", ...
    "MATLABName", "clib.ctMatlabInterface.ct_appdelete", ...
    "Description", "clib.ctMatlabInterface.ct_appdelete Representation of C++ function ct_appdelete."); % Modify help description values as needed.
defineOutput(ct_appdeleteDefinition, "RetVal", "int32");
validate(ct_appdeleteDefinition);

%% C++ function |soln_newSolution| with MATLAB name |clib.ctMatlabInterface.soln_newSolution|
% C++ Signature: int soln_newSolution(char const * infile,char const * name,char const * transport)

soln_newSolutionDefinition = addFunction(libDef, ...
    "int soln_newSolution(char const * infile,char const * name,char const * transport)", ...
    "MATLABName", "clib.ctMatlabInterface.soln_newSolution", ...
    "Description", "clib.ctMatlabInterface.soln_newSolution Representation of C++ function soln_newSolution."); % Modify help description values as needed.
defineArgument(soln_newSolutionDefinition, "infile", "string", "input", "nullTerminated");
defineArgument(soln_newSolutionDefinition, "name", "string", "input", "nullTerminated");
defineArgument(soln_newSolutionDefinition, "transport", "string", "input", "nullTerminated");
defineOutput(soln_newSolutionDefinition, "RetVal", "int32");
validate(soln_newSolutionDefinition);

%% C++ function |soln_newInterface| with MATLAB name |clib.ctMatlabInterface.soln_newInterface|
% C++ Signature: int soln_newInterface(char const * infile,char const * name,int na,int const * adjacent)

soln_newInterfaceDefinition = addFunction(libDef, ...
    "int soln_newInterface(char const * infile,char const * name,int na,int const * adjacent)", ...
    "MATLABName", "clib.ctMatlabInterface.soln_newInterface", ...
    "Description", "clib.ctMatlabInterface.soln_newInterface Representation of C++ function soln_newInterface."); % Modify help description values as needed.
defineArgument(soln_newInterfaceDefinition, "infile", "string", "input", "nullTerminated");
defineArgument(soln_newInterfaceDefinition, "name", "string", "input", "nullTerminated");
defineArgument(soln_newInterfaceDefinition, "na", "int32");
defineArgument(soln_newInterfaceDefinition, "adjacent", "clib.array.ctMatlabInterface.Int", "input", "na"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Int", or "int32"
defineOutput(soln_newInterfaceDefinition, "RetVal", "int32");
validate(soln_newInterfaceDefinition);

%% C++ function |soln_del| with MATLAB name |clib.ctMatlabInterface.soln_del|
% C++ Signature: int soln_del(int n)

soln_delDefinition = addFunction(libDef, ...
    "int soln_del(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.soln_del", ...
    "Description", "clib.ctMatlabInterface.soln_del Representation of C++ function soln_del." + newline + ...
    "note that linked objects are deleted as well"); % Modify help description values as needed.
defineArgument(soln_delDefinition, "n", "int32");
defineOutput(soln_delDefinition, "RetVal", "int32");
validate(soln_delDefinition);

%% C++ function |soln_name| with MATLAB name |clib.ctMatlabInterface.soln_name|
% C++ Signature: int soln_name(int n,int buflen,char * buf)

soln_nameDefinition = addFunction(libDef, ...
    "int soln_name(int n,int buflen,char * buf)", ...
    "MATLABName", "clib.ctMatlabInterface.soln_name", ...
    "Description", "clib.ctMatlabInterface.soln_name Representation of C++ function soln_name."); % Modify help description values as needed.
defineArgument(soln_nameDefinition, "n", "int32");
defineArgument(soln_nameDefinition, "buflen", "int32");
defineArgument(soln_nameDefinition, "buf", "clib.array.ctMatlabInterface.Char", "input", "buflen"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(soln_nameDefinition, "RetVal", "int32");
validate(soln_nameDefinition);

%% C++ function |soln_thermo| with MATLAB name |clib.ctMatlabInterface.soln_thermo|
% C++ Signature: int soln_thermo(int n)

soln_thermoDefinition = addFunction(libDef, ...
    "int soln_thermo(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.soln_thermo", ...
    "Description", "clib.ctMatlabInterface.soln_thermo Representation of C++ function soln_thermo."); % Modify help description values as needed.
defineArgument(soln_thermoDefinition, "n", "int32");
defineOutput(soln_thermoDefinition, "RetVal", "int32");
validate(soln_thermoDefinition);

%% C++ function |soln_kinetics| with MATLAB name |clib.ctMatlabInterface.soln_kinetics|
% C++ Signature: int soln_kinetics(int n)

soln_kineticsDefinition = addFunction(libDef, ...
    "int soln_kinetics(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.soln_kinetics", ...
    "Description", "clib.ctMatlabInterface.soln_kinetics Representation of C++ function soln_kinetics."); % Modify help description values as needed.
defineArgument(soln_kineticsDefinition, "n", "int32");
defineOutput(soln_kineticsDefinition, "RetVal", "int32");
validate(soln_kineticsDefinition);

%% C++ function |soln_transport| with MATLAB name |clib.ctMatlabInterface.soln_transport|
% C++ Signature: int soln_transport(int n)

soln_transportDefinition = addFunction(libDef, ...
    "int soln_transport(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.soln_transport", ...
    "Description", "clib.ctMatlabInterface.soln_transport Representation of C++ function soln_transport."); % Modify help description values as needed.
defineArgument(soln_transportDefinition, "n", "int32");
defineOutput(soln_transportDefinition, "RetVal", "int32");
validate(soln_transportDefinition);

%% C++ function |soln_setTransportModel| with MATLAB name |clib.ctMatlabInterface.soln_setTransportModel|
% C++ Signature: int soln_setTransportModel(int n,char const * model)

soln_setTransportModelDefinition = addFunction(libDef, ...
    "int soln_setTransportModel(int n,char const * model)", ...
    "MATLABName", "clib.ctMatlabInterface.soln_setTransportModel", ...
    "Description", "clib.ctMatlabInterface.soln_setTransportModel Representation of C++ function soln_setTransportModel."); % Modify help description values as needed.
defineArgument(soln_setTransportModelDefinition, "n", "int32");
defineArgument(soln_setTransportModelDefinition, "model", "string", "input", "nullTerminated");
defineOutput(soln_setTransportModelDefinition, "RetVal", "int32");
validate(soln_setTransportModelDefinition);

%% C++ function |soln_nAdjacent| with MATLAB name |clib.ctMatlabInterface.soln_nAdjacent|
% C++ Signature: size_t soln_nAdjacent(int n)

soln_nAdjacentDefinition = addFunction(libDef, ...
    "size_t soln_nAdjacent(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.soln_nAdjacent", ...
    "Description", "clib.ctMatlabInterface.soln_nAdjacent Representation of C++ function soln_nAdjacent."); % Modify help description values as needed.
defineArgument(soln_nAdjacentDefinition, "n", "int32");
defineOutput(soln_nAdjacentDefinition, "RetVal", "uint64");
validate(soln_nAdjacentDefinition);

%% C++ function |soln_adjacent| with MATLAB name |clib.ctMatlabInterface.soln_adjacent|
% C++ Signature: int soln_adjacent(int n,int a)

soln_adjacentDefinition = addFunction(libDef, ...
    "int soln_adjacent(int n,int a)", ...
    "MATLABName", "clib.ctMatlabInterface.soln_adjacent", ...
    "Description", "clib.ctMatlabInterface.soln_adjacent Representation of C++ function soln_adjacent."); % Modify help description values as needed.
defineArgument(soln_adjacentDefinition, "n", "int32");
defineArgument(soln_adjacentDefinition, "a", "int32");
defineOutput(soln_adjacentDefinition, "RetVal", "int32");
validate(soln_adjacentDefinition);

%% C++ function |soln_adjacentName| with MATLAB name |clib.ctMatlabInterface.soln_adjacentName|
% C++ Signature: int soln_adjacentName(int n,int a,int lennm,char * nm)

soln_adjacentNameDefinition = addFunction(libDef, ...
    "int soln_adjacentName(int n,int a,int lennm,char * nm)", ...
    "MATLABName", "clib.ctMatlabInterface.soln_adjacentName", ...
    "Description", "clib.ctMatlabInterface.soln_adjacentName Representation of C++ function soln_adjacentName."); % Modify help description values as needed.
defineArgument(soln_adjacentNameDefinition, "n", "int32");
defineArgument(soln_adjacentNameDefinition, "a", "int32");
defineArgument(soln_adjacentNameDefinition, "lennm", "int32");
defineArgument(soln_adjacentNameDefinition, "nm", "clib.array.ctMatlabInterface.Char", "input", "lennm"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(soln_adjacentNameDefinition, "RetVal", "int32");
validate(soln_adjacentNameDefinition);

%% C++ function |thermo_parent| with MATLAB name |clib.ctMatlabInterface.thermo_parent|
% C++ Signature: int thermo_parent(int n)

thermo_parentDefinition = addFunction(libDef, ...
    "int thermo_parent(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_parent", ...
    "Description", "clib.ctMatlabInterface.thermo_parent Representation of C++ function thermo_parent."); % Modify help description values as needed.
defineArgument(thermo_parentDefinition, "n", "int32");
defineOutput(thermo_parentDefinition, "RetVal", "int32");
validate(thermo_parentDefinition);

%% C++ function |thermo_size| with MATLAB name |clib.ctMatlabInterface.thermo_size|
% C++ Signature: int thermo_size()

thermo_sizeDefinition = addFunction(libDef, ...
    "int thermo_size()", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_size", ...
    "Description", "clib.ctMatlabInterface.thermo_size Representation of C++ function thermo_size."); % Modify help description values as needed.
defineOutput(thermo_sizeDefinition, "RetVal", "int32");
validate(thermo_sizeDefinition);

%% C++ function |thermo_del| with MATLAB name |clib.ctMatlabInterface.thermo_del|
% C++ Signature: int thermo_del(int n)

thermo_delDefinition = addFunction(libDef, ...
    "int thermo_del(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_del", ...
    "Description", "clib.ctMatlabInterface.thermo_del Representation of C++ function thermo_del." + newline + ...
    "no-op; object is managed by Solution"); % Modify help description values as needed.
defineArgument(thermo_delDefinition, "n", "int32");
defineOutput(thermo_delDefinition, "RetVal", "int32");
validate(thermo_delDefinition);

%% C++ function |thermo_nElements| with MATLAB name |clib.ctMatlabInterface.thermo_nElements|
% C++ Signature: size_t thermo_nElements(int n)

thermo_nElementsDefinition = addFunction(libDef, ...
    "size_t thermo_nElements(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_nElements", ...
    "Description", "clib.ctMatlabInterface.thermo_nElements Representation of C++ function thermo_nElements."); % Modify help description values as needed.
defineArgument(thermo_nElementsDefinition, "n", "int32");
defineOutput(thermo_nElementsDefinition, "RetVal", "uint64");
validate(thermo_nElementsDefinition);

%% C++ function |thermo_nSpecies| with MATLAB name |clib.ctMatlabInterface.thermo_nSpecies|
% C++ Signature: size_t thermo_nSpecies(int n)

thermo_nSpeciesDefinition = addFunction(libDef, ...
    "size_t thermo_nSpecies(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_nSpecies", ...
    "Description", "clib.ctMatlabInterface.thermo_nSpecies Representation of C++ function thermo_nSpecies."); % Modify help description values as needed.
defineArgument(thermo_nSpeciesDefinition, "n", "int32");
defineOutput(thermo_nSpeciesDefinition, "RetVal", "uint64");
validate(thermo_nSpeciesDefinition);

%% C++ function |thermo_temperature| with MATLAB name |clib.ctMatlabInterface.thermo_temperature|
% C++ Signature: double thermo_temperature(int n)

thermo_temperatureDefinition = addFunction(libDef, ...
    "double thermo_temperature(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_temperature", ...
    "Description", "clib.ctMatlabInterface.thermo_temperature Representation of C++ function thermo_temperature."); % Modify help description values as needed.
defineArgument(thermo_temperatureDefinition, "n", "int32");
defineOutput(thermo_temperatureDefinition, "RetVal", "double");
validate(thermo_temperatureDefinition);

%% C++ function |thermo_setTemperature| with MATLAB name |clib.ctMatlabInterface.thermo_setTemperature|
% C++ Signature: int thermo_setTemperature(int n,double t)

thermo_setTemperatureDefinition = addFunction(libDef, ...
    "int thermo_setTemperature(int n,double t)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_setTemperature", ...
    "Description", "clib.ctMatlabInterface.thermo_setTemperature Representation of C++ function thermo_setTemperature."); % Modify help description values as needed.
defineArgument(thermo_setTemperatureDefinition, "n", "int32");
defineArgument(thermo_setTemperatureDefinition, "t", "double");
defineOutput(thermo_setTemperatureDefinition, "RetVal", "int32");
validate(thermo_setTemperatureDefinition);

%% C++ function |thermo_density| with MATLAB name |clib.ctMatlabInterface.thermo_density|
% C++ Signature: double thermo_density(int n)

thermo_densityDefinition = addFunction(libDef, ...
    "double thermo_density(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_density", ...
    "Description", "clib.ctMatlabInterface.thermo_density Representation of C++ function thermo_density."); % Modify help description values as needed.
defineArgument(thermo_densityDefinition, "n", "int32");
defineOutput(thermo_densityDefinition, "RetVal", "double");
validate(thermo_densityDefinition);

%% C++ function |thermo_setDensity| with MATLAB name |clib.ctMatlabInterface.thermo_setDensity|
% C++ Signature: int thermo_setDensity(int n,double rho)

thermo_setDensityDefinition = addFunction(libDef, ...
    "int thermo_setDensity(int n,double rho)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_setDensity", ...
    "Description", "clib.ctMatlabInterface.thermo_setDensity Representation of C++ function thermo_setDensity."); % Modify help description values as needed.
defineArgument(thermo_setDensityDefinition, "n", "int32");
defineArgument(thermo_setDensityDefinition, "rho", "double");
defineOutput(thermo_setDensityDefinition, "RetVal", "int32");
validate(thermo_setDensityDefinition);

%% C++ function |thermo_molarDensity| with MATLAB name |clib.ctMatlabInterface.thermo_molarDensity|
% C++ Signature: double thermo_molarDensity(int n)

thermo_molarDensityDefinition = addFunction(libDef, ...
    "double thermo_molarDensity(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_molarDensity", ...
    "Description", "clib.ctMatlabInterface.thermo_molarDensity Representation of C++ function thermo_molarDensity."); % Modify help description values as needed.
defineArgument(thermo_molarDensityDefinition, "n", "int32");
defineOutput(thermo_molarDensityDefinition, "RetVal", "double");
validate(thermo_molarDensityDefinition);

%% C++ function |thermo_meanMolecularWeight| with MATLAB name |clib.ctMatlabInterface.thermo_meanMolecularWeight|
% C++ Signature: double thermo_meanMolecularWeight(int n)

thermo_meanMolecularWeightDefinition = addFunction(libDef, ...
    "double thermo_meanMolecularWeight(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_meanMolecularWeight", ...
    "Description", "clib.ctMatlabInterface.thermo_meanMolecularWeight Representation of C++ function thermo_meanMolecularWeight."); % Modify help description values as needed.
defineArgument(thermo_meanMolecularWeightDefinition, "n", "int32");
defineOutput(thermo_meanMolecularWeightDefinition, "RetVal", "double");
validate(thermo_meanMolecularWeightDefinition);

%% C++ function |thermo_moleFraction| with MATLAB name |clib.ctMatlabInterface.thermo_moleFraction|
% C++ Signature: double thermo_moleFraction(int n,size_t k)

thermo_moleFractionDefinition = addFunction(libDef, ...
    "double thermo_moleFraction(int n,size_t k)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_moleFraction", ...
    "Description", "clib.ctMatlabInterface.thermo_moleFraction Representation of C++ function thermo_moleFraction."); % Modify help description values as needed.
defineArgument(thermo_moleFractionDefinition, "n", "int32");
defineArgument(thermo_moleFractionDefinition, "k", "uint64");
defineOutput(thermo_moleFractionDefinition, "RetVal", "double");
validate(thermo_moleFractionDefinition);

%% C++ function |thermo_massFraction| with MATLAB name |clib.ctMatlabInterface.thermo_massFraction|
% C++ Signature: double thermo_massFraction(int n,size_t k)

thermo_massFractionDefinition = addFunction(libDef, ...
    "double thermo_massFraction(int n,size_t k)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_massFraction", ...
    "Description", "clib.ctMatlabInterface.thermo_massFraction Representation of C++ function thermo_massFraction."); % Modify help description values as needed.
defineArgument(thermo_massFractionDefinition, "n", "int32");
defineArgument(thermo_massFractionDefinition, "k", "uint64");
defineOutput(thermo_massFractionDefinition, "RetVal", "double");
validate(thermo_massFractionDefinition);

%% C++ function |thermo_getMoleFractions| with MATLAB name |clib.ctMatlabInterface.thermo_getMoleFractions|
% C++ Signature: int thermo_getMoleFractions(int n,size_t lenx,double * x)

thermo_getMoleFractionsDefinition = addFunction(libDef, ...
    "int thermo_getMoleFractions(int n,size_t lenx,double * x)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_getMoleFractions", ...
    "Description", "clib.ctMatlabInterface.thermo_getMoleFractions Representation of C++ function thermo_getMoleFractions."); % Modify help description values as needed.
defineArgument(thermo_getMoleFractionsDefinition, "n", "int32");
defineArgument(thermo_getMoleFractionsDefinition, "lenx", "uint64");
defineArgument(thermo_getMoleFractionsDefinition, "x", "clib.array.ctMatlabInterface.Double", "input", "lenx"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_getMoleFractionsDefinition, "RetVal", "int32");
validate(thermo_getMoleFractionsDefinition);

%% C++ function |thermo_getMassFractions| with MATLAB name |clib.ctMatlabInterface.thermo_getMassFractions|
% C++ Signature: int thermo_getMassFractions(int n,size_t leny,double * y)

thermo_getMassFractionsDefinition = addFunction(libDef, ...
    "int thermo_getMassFractions(int n,size_t leny,double * y)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_getMassFractions", ...
    "Description", "clib.ctMatlabInterface.thermo_getMassFractions Representation of C++ function thermo_getMassFractions."); % Modify help description values as needed.
defineArgument(thermo_getMassFractionsDefinition, "n", "int32");
defineArgument(thermo_getMassFractionsDefinition, "leny", "uint64");
defineArgument(thermo_getMassFractionsDefinition, "y", "clib.array.ctMatlabInterface.Double", "input", "leny"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_getMassFractionsDefinition, "RetVal", "int32");
validate(thermo_getMassFractionsDefinition);

%% C++ function |thermo_setMoleFractions| with MATLAB name |clib.ctMatlabInterface.thermo_setMoleFractions|
% C++ Signature: int thermo_setMoleFractions(int n,size_t lenx,double * x,int norm)

thermo_setMoleFractionsDefinition = addFunction(libDef, ...
    "int thermo_setMoleFractions(int n,size_t lenx,double * x,int norm)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_setMoleFractions", ...
    "Description", "clib.ctMatlabInterface.thermo_setMoleFractions Representation of C++ function thermo_setMoleFractions."); % Modify help description values as needed.
defineArgument(thermo_setMoleFractionsDefinition, "n", "int32");
defineArgument(thermo_setMoleFractionsDefinition, "lenx", "uint64");
defineArgument(thermo_setMoleFractionsDefinition, "x", "clib.array.ctMatlabInterface.Double", "input", "lenx"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineArgument(thermo_setMoleFractionsDefinition, "norm", "int32");
defineOutput(thermo_setMoleFractionsDefinition, "RetVal", "int32");
validate(thermo_setMoleFractionsDefinition);

%% C++ function |thermo_setMassFractions| with MATLAB name |clib.ctMatlabInterface.thermo_setMassFractions|
% C++ Signature: int thermo_setMassFractions(int n,size_t leny,double * y,int norm)

thermo_setMassFractionsDefinition = addFunction(libDef, ...
    "int thermo_setMassFractions(int n,size_t leny,double * y,int norm)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_setMassFractions", ...
    "Description", "clib.ctMatlabInterface.thermo_setMassFractions Representation of C++ function thermo_setMassFractions."); % Modify help description values as needed.
defineArgument(thermo_setMassFractionsDefinition, "n", "int32");
defineArgument(thermo_setMassFractionsDefinition, "leny", "uint64");
defineArgument(thermo_setMassFractionsDefinition, "y", "clib.array.ctMatlabInterface.Double", "input", "leny"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineArgument(thermo_setMassFractionsDefinition, "norm", "int32");
defineOutput(thermo_setMassFractionsDefinition, "RetVal", "int32");
validate(thermo_setMassFractionsDefinition);

%% C++ function |thermo_setMoleFractionsByName| with MATLAB name |clib.ctMatlabInterface.thermo_setMoleFractionsByName|
% C++ Signature: int thermo_setMoleFractionsByName(int n,char const * x)

thermo_setMoleFractionsByNameDefinition = addFunction(libDef, ...
    "int thermo_setMoleFractionsByName(int n,char const * x)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_setMoleFractionsByName", ...
    "Description", "clib.ctMatlabInterface.thermo_setMoleFractionsByName Representation of C++ function thermo_setMoleFractionsByName."); % Modify help description values as needed.
defineArgument(thermo_setMoleFractionsByNameDefinition, "n", "int32");
defineArgument(thermo_setMoleFractionsByNameDefinition, "x", "string", "input", "nullTerminated");
defineOutput(thermo_setMoleFractionsByNameDefinition, "RetVal", "int32");
validate(thermo_setMoleFractionsByNameDefinition);

%% C++ function |thermo_setMassFractionsByName| with MATLAB name |clib.ctMatlabInterface.thermo_setMassFractionsByName|
% C++ Signature: int thermo_setMassFractionsByName(int n,char const * y)

thermo_setMassFractionsByNameDefinition = addFunction(libDef, ...
    "int thermo_setMassFractionsByName(int n,char const * y)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_setMassFractionsByName", ...
    "Description", "clib.ctMatlabInterface.thermo_setMassFractionsByName Representation of C++ function thermo_setMassFractionsByName."); % Modify help description values as needed.
defineArgument(thermo_setMassFractionsByNameDefinition, "n", "int32");
defineArgument(thermo_setMassFractionsByNameDefinition, "y", "string", "input", "nullTerminated");
defineOutput(thermo_setMassFractionsByNameDefinition, "RetVal", "int32");
validate(thermo_setMassFractionsByNameDefinition);

%% C++ function |thermo_getAtomicWeights| with MATLAB name |clib.ctMatlabInterface.thermo_getAtomicWeights|
% C++ Signature: int thermo_getAtomicWeights(int n,size_t lenm,double * atw)

thermo_getAtomicWeightsDefinition = addFunction(libDef, ...
    "int thermo_getAtomicWeights(int n,size_t lenm,double * atw)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_getAtomicWeights", ...
    "Description", "clib.ctMatlabInterface.thermo_getAtomicWeights Representation of C++ function thermo_getAtomicWeights."); % Modify help description values as needed.
defineArgument(thermo_getAtomicWeightsDefinition, "n", "int32");
defineArgument(thermo_getAtomicWeightsDefinition, "lenm", "uint64");
defineArgument(thermo_getAtomicWeightsDefinition, "atw", "clib.array.ctMatlabInterface.Double", "input", "lenm"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_getAtomicWeightsDefinition, "RetVal", "int32");
validate(thermo_getAtomicWeightsDefinition);

%% C++ function |thermo_getMolecularWeights| with MATLAB name |clib.ctMatlabInterface.thermo_getMolecularWeights|
% C++ Signature: int thermo_getMolecularWeights(int n,size_t lenm,double * mw)

thermo_getMolecularWeightsDefinition = addFunction(libDef, ...
    "int thermo_getMolecularWeights(int n,size_t lenm,double * mw)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_getMolecularWeights", ...
    "Description", "clib.ctMatlabInterface.thermo_getMolecularWeights Representation of C++ function thermo_getMolecularWeights."); % Modify help description values as needed.
defineArgument(thermo_getMolecularWeightsDefinition, "n", "int32");
defineArgument(thermo_getMolecularWeightsDefinition, "lenm", "uint64");
defineArgument(thermo_getMolecularWeightsDefinition, "mw", "clib.array.ctMatlabInterface.Double", "input", "lenm"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_getMolecularWeightsDefinition, "RetVal", "int32");
validate(thermo_getMolecularWeightsDefinition);

%% C++ function |thermo_getCharges| with MATLAB name |clib.ctMatlabInterface.thermo_getCharges|
% C++ Signature: int thermo_getCharges(int n,size_t lenm,double * sc)

thermo_getChargesDefinition = addFunction(libDef, ...
    "int thermo_getCharges(int n,size_t lenm,double * sc)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_getCharges", ...
    "Description", "clib.ctMatlabInterface.thermo_getCharges Representation of C++ function thermo_getCharges."); % Modify help description values as needed.
defineArgument(thermo_getChargesDefinition, "n", "int32");
defineArgument(thermo_getChargesDefinition, "lenm", "uint64");
defineArgument(thermo_getChargesDefinition, "sc", "clib.array.ctMatlabInterface.Double", "input", "lenm"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_getChargesDefinition, "RetVal", "int32");
validate(thermo_getChargesDefinition);

%% C++ function |thermo_getElementName| with MATLAB name |clib.ctMatlabInterface.thermo_getElementName|
% C++ Signature: int thermo_getElementName(int n,size_t k,size_t lennm,char * nm)

thermo_getElementNameDefinition = addFunction(libDef, ...
    "int thermo_getElementName(int n,size_t k,size_t lennm,char * nm)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_getElementName", ...
    "Description", "clib.ctMatlabInterface.thermo_getElementName Representation of C++ function thermo_getElementName."); % Modify help description values as needed.
defineArgument(thermo_getElementNameDefinition, "n", "int32");
defineArgument(thermo_getElementNameDefinition, "k", "uint64");
defineArgument(thermo_getElementNameDefinition, "lennm", "uint64");
defineArgument(thermo_getElementNameDefinition, "nm", "clib.array.ctMatlabInterface.Char", "input", "lennm"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(thermo_getElementNameDefinition, "RetVal", "int32");
validate(thermo_getElementNameDefinition);

%% C++ function |thermo_getSpeciesName| with MATLAB name |clib.ctMatlabInterface.thermo_getSpeciesName|
% C++ Signature: int thermo_getSpeciesName(int n,size_t m,size_t lennm,char * nm)

thermo_getSpeciesNameDefinition = addFunction(libDef, ...
    "int thermo_getSpeciesName(int n,size_t m,size_t lennm,char * nm)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_getSpeciesName", ...
    "Description", "clib.ctMatlabInterface.thermo_getSpeciesName Representation of C++ function thermo_getSpeciesName."); % Modify help description values as needed.
defineArgument(thermo_getSpeciesNameDefinition, "n", "int32");
defineArgument(thermo_getSpeciesNameDefinition, "m", "uint64");
defineArgument(thermo_getSpeciesNameDefinition, "lennm", "uint64");
defineArgument(thermo_getSpeciesNameDefinition, "nm", "clib.array.ctMatlabInterface.Char", "input", "lennm"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(thermo_getSpeciesNameDefinition, "RetVal", "int32");
validate(thermo_getSpeciesNameDefinition);

%% C++ function |thermo_getName| with MATLAB name |clib.ctMatlabInterface.thermo_getName|
% C++ Signature: int thermo_getName(int n,size_t lennm,char * nm)

thermo_getNameDefinition = addFunction(libDef, ...
    "int thermo_getName(int n,size_t lennm,char * nm)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_getName", ...
    "Description", "clib.ctMatlabInterface.thermo_getName Representation of C++ function thermo_getName."); % Modify help description values as needed.
defineArgument(thermo_getNameDefinition, "n", "int32");
defineArgument(thermo_getNameDefinition, "lennm", "uint64");
defineArgument(thermo_getNameDefinition, "nm", "clib.array.ctMatlabInterface.Char", "input", "lennm"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(thermo_getNameDefinition, "RetVal", "int32");
validate(thermo_getNameDefinition);

%% C++ function |thermo_setName| with MATLAB name |clib.ctMatlabInterface.thermo_setName|
% C++ Signature: int thermo_setName(int n,char const * nm)

thermo_setNameDefinition = addFunction(libDef, ...
    "int thermo_setName(int n,char const * nm)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_setName", ...
    "Description", "clib.ctMatlabInterface.thermo_setName Representation of C++ function thermo_setName."); % Modify help description values as needed.
defineArgument(thermo_setNameDefinition, "n", "int32");
defineArgument(thermo_setNameDefinition, "nm", "string", "input", "nullTerminated");
defineOutput(thermo_setNameDefinition, "RetVal", "int32");
validate(thermo_setNameDefinition);

%% C++ function |thermo_elementIndex| with MATLAB name |clib.ctMatlabInterface.thermo_elementIndex|
% C++ Signature: size_t thermo_elementIndex(int n,char const * nm)

thermo_elementIndexDefinition = addFunction(libDef, ...
    "size_t thermo_elementIndex(int n,char const * nm)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_elementIndex", ...
    "Description", "clib.ctMatlabInterface.thermo_elementIndex Representation of C++ function thermo_elementIndex."); % Modify help description values as needed.
defineArgument(thermo_elementIndexDefinition, "n", "int32");
defineArgument(thermo_elementIndexDefinition, "nm", "string", "input", "nullTerminated");
defineOutput(thermo_elementIndexDefinition, "RetVal", "uint64");
validate(thermo_elementIndexDefinition);

%% C++ function |thermo_speciesIndex| with MATLAB name |clib.ctMatlabInterface.thermo_speciesIndex|
% C++ Signature: size_t thermo_speciesIndex(int n,char const * nm)

thermo_speciesIndexDefinition = addFunction(libDef, ...
    "size_t thermo_speciesIndex(int n,char const * nm)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_speciesIndex", ...
    "Description", "clib.ctMatlabInterface.thermo_speciesIndex Representation of C++ function thermo_speciesIndex."); % Modify help description values as needed.
defineArgument(thermo_speciesIndexDefinition, "n", "int32");
defineArgument(thermo_speciesIndexDefinition, "nm", "string", "input", "nullTerminated");
defineOutput(thermo_speciesIndexDefinition, "RetVal", "uint64");
validate(thermo_speciesIndexDefinition);

%% C++ function |thermo_report| with MATLAB name |clib.ctMatlabInterface.thermo_report|
% C++ Signature: int thermo_report(int nth,int show_thermo,double threshold,int ibuf,char * buf)

thermo_reportDefinition = addFunction(libDef, ...
    "int thermo_report(int nth,int show_thermo,double threshold,int ibuf,char * buf)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_report", ...
    "Description", "clib.ctMatlabInterface.thermo_report Representation of C++ function thermo_report."); % Modify help description values as needed.
defineArgument(thermo_reportDefinition, "nth", "int32");
defineArgument(thermo_reportDefinition, "show_thermo", "int32");
defineArgument(thermo_reportDefinition, "threshold", "double");
defineArgument(thermo_reportDefinition, "ibuf", "int32");
defineArgument(thermo_reportDefinition, "buf", "clib.array.ctMatlabInterface.Char", "input", "ibuf"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(thermo_reportDefinition, "RetVal", "int32");
validate(thermo_reportDefinition);

%% C++ function |thermo_print| with MATLAB name |clib.ctMatlabInterface.thermo_print|
% C++ Signature: int thermo_print(int nth,int show_thermo,double threshold)

thermo_printDefinition = addFunction(libDef, ...
    "int thermo_print(int nth,int show_thermo,double threshold)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_print", ...
    "Description", "clib.ctMatlabInterface.thermo_print Representation of C++ function thermo_print."); % Modify help description values as needed.
defineArgument(thermo_printDefinition, "nth", "int32");
defineArgument(thermo_printDefinition, "show_thermo", "int32");
defineArgument(thermo_printDefinition, "threshold", "double");
defineOutput(thermo_printDefinition, "RetVal", "int32");
validate(thermo_printDefinition);

%% C++ function |thermo_nAtoms| with MATLAB name |clib.ctMatlabInterface.thermo_nAtoms|
% C++ Signature: double thermo_nAtoms(int n,size_t k,size_t m)

thermo_nAtomsDefinition = addFunction(libDef, ...
    "double thermo_nAtoms(int n,size_t k,size_t m)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_nAtoms", ...
    "Description", "clib.ctMatlabInterface.thermo_nAtoms Representation of C++ function thermo_nAtoms."); % Modify help description values as needed.
defineArgument(thermo_nAtomsDefinition, "n", "int32");
defineArgument(thermo_nAtomsDefinition, "k", "uint64");
defineArgument(thermo_nAtomsDefinition, "m", "uint64");
defineOutput(thermo_nAtomsDefinition, "RetVal", "double");
validate(thermo_nAtomsDefinition);

%% C++ function |thermo_addElement| with MATLAB name |clib.ctMatlabInterface.thermo_addElement|
% C++ Signature: int thermo_addElement(int n,char const * name,double weight)

thermo_addElementDefinition = addFunction(libDef, ...
    "int thermo_addElement(int n,char const * name,double weight)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_addElement", ...
    "Description", "clib.ctMatlabInterface.thermo_addElement Representation of C++ function thermo_addElement."); % Modify help description values as needed.
defineArgument(thermo_addElementDefinition, "n", "int32");
defineArgument(thermo_addElementDefinition, "name", "string", "input", "nullTerminated");
defineArgument(thermo_addElementDefinition, "weight", "double");
defineOutput(thermo_addElementDefinition, "RetVal", "int32");
validate(thermo_addElementDefinition);

%% C++ function |thermo_getEosType| with MATLAB name |clib.ctMatlabInterface.thermo_getEosType|
% C++ Signature: int thermo_getEosType(int n,size_t leneos,char * eos)

thermo_getEosTypeDefinition = addFunction(libDef, ...
    "int thermo_getEosType(int n,size_t leneos,char * eos)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_getEosType", ...
    "Description", "clib.ctMatlabInterface.thermo_getEosType Representation of C++ function thermo_getEosType."); % Modify help description values as needed.
defineArgument(thermo_getEosTypeDefinition, "n", "int32");
defineArgument(thermo_getEosTypeDefinition, "leneos", "uint64");
defineArgument(thermo_getEosTypeDefinition, "eos", "clib.array.ctMatlabInterface.Char", "input", "leneos"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(thermo_getEosTypeDefinition, "RetVal", "int32");
validate(thermo_getEosTypeDefinition);

%% C++ function |thermo_refPressure| with MATLAB name |clib.ctMatlabInterface.thermo_refPressure|
% C++ Signature: double thermo_refPressure(int n)

thermo_refPressureDefinition = addFunction(libDef, ...
    "double thermo_refPressure(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_refPressure", ...
    "Description", "clib.ctMatlabInterface.thermo_refPressure Representation of C++ function thermo_refPressure."); % Modify help description values as needed.
defineArgument(thermo_refPressureDefinition, "n", "int32");
defineOutput(thermo_refPressureDefinition, "RetVal", "double");
validate(thermo_refPressureDefinition);

%% C++ function |thermo_minTemp| with MATLAB name |clib.ctMatlabInterface.thermo_minTemp|
% C++ Signature: double thermo_minTemp(int n,int k)

thermo_minTempDefinition = addFunction(libDef, ...
    "double thermo_minTemp(int n,int k)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_minTemp", ...
    "Description", "clib.ctMatlabInterface.thermo_minTemp Representation of C++ function thermo_minTemp."); % Modify help description values as needed.
defineArgument(thermo_minTempDefinition, "n", "int32");
defineArgument(thermo_minTempDefinition, "k", "int32");
defineOutput(thermo_minTempDefinition, "RetVal", "double");
validate(thermo_minTempDefinition);

%% C++ function |thermo_maxTemp| with MATLAB name |clib.ctMatlabInterface.thermo_maxTemp|
% C++ Signature: double thermo_maxTemp(int n,int k)

thermo_maxTempDefinition = addFunction(libDef, ...
    "double thermo_maxTemp(int n,int k)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_maxTemp", ...
    "Description", "clib.ctMatlabInterface.thermo_maxTemp Representation of C++ function thermo_maxTemp."); % Modify help description values as needed.
defineArgument(thermo_maxTempDefinition, "n", "int32");
defineArgument(thermo_maxTempDefinition, "k", "int32");
defineOutput(thermo_maxTempDefinition, "RetVal", "double");
validate(thermo_maxTempDefinition);

%% C++ function |thermo_enthalpy_mole| with MATLAB name |clib.ctMatlabInterface.thermo_enthalpy_mole|
% C++ Signature: double thermo_enthalpy_mole(int n)

thermo_enthalpy_moleDefinition = addFunction(libDef, ...
    "double thermo_enthalpy_mole(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_enthalpy_mole", ...
    "Description", "clib.ctMatlabInterface.thermo_enthalpy_mole Representation of C++ function thermo_enthalpy_mole."); % Modify help description values as needed.
defineArgument(thermo_enthalpy_moleDefinition, "n", "int32");
defineOutput(thermo_enthalpy_moleDefinition, "RetVal", "double");
validate(thermo_enthalpy_moleDefinition);

%% C++ function |thermo_intEnergy_mole| with MATLAB name |clib.ctMatlabInterface.thermo_intEnergy_mole|
% C++ Signature: double thermo_intEnergy_mole(int n)

thermo_intEnergy_moleDefinition = addFunction(libDef, ...
    "double thermo_intEnergy_mole(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_intEnergy_mole", ...
    "Description", "clib.ctMatlabInterface.thermo_intEnergy_mole Representation of C++ function thermo_intEnergy_mole."); % Modify help description values as needed.
defineArgument(thermo_intEnergy_moleDefinition, "n", "int32");
defineOutput(thermo_intEnergy_moleDefinition, "RetVal", "double");
validate(thermo_intEnergy_moleDefinition);

%% C++ function |thermo_entropy_mole| with MATLAB name |clib.ctMatlabInterface.thermo_entropy_mole|
% C++ Signature: double thermo_entropy_mole(int n)

thermo_entropy_moleDefinition = addFunction(libDef, ...
    "double thermo_entropy_mole(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_entropy_mole", ...
    "Description", "clib.ctMatlabInterface.thermo_entropy_mole Representation of C++ function thermo_entropy_mole."); % Modify help description values as needed.
defineArgument(thermo_entropy_moleDefinition, "n", "int32");
defineOutput(thermo_entropy_moleDefinition, "RetVal", "double");
validate(thermo_entropy_moleDefinition);

%% C++ function |thermo_gibbs_mole| with MATLAB name |clib.ctMatlabInterface.thermo_gibbs_mole|
% C++ Signature: double thermo_gibbs_mole(int n)

thermo_gibbs_moleDefinition = addFunction(libDef, ...
    "double thermo_gibbs_mole(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_gibbs_mole", ...
    "Description", "clib.ctMatlabInterface.thermo_gibbs_mole Representation of C++ function thermo_gibbs_mole."); % Modify help description values as needed.
defineArgument(thermo_gibbs_moleDefinition, "n", "int32");
defineOutput(thermo_gibbs_moleDefinition, "RetVal", "double");
validate(thermo_gibbs_moleDefinition);

%% C++ function |thermo_cp_mole| with MATLAB name |clib.ctMatlabInterface.thermo_cp_mole|
% C++ Signature: double thermo_cp_mole(int n)

thermo_cp_moleDefinition = addFunction(libDef, ...
    "double thermo_cp_mole(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_cp_mole", ...
    "Description", "clib.ctMatlabInterface.thermo_cp_mole Representation of C++ function thermo_cp_mole."); % Modify help description values as needed.
defineArgument(thermo_cp_moleDefinition, "n", "int32");
defineOutput(thermo_cp_moleDefinition, "RetVal", "double");
validate(thermo_cp_moleDefinition);

%% C++ function |thermo_cv_mole| with MATLAB name |clib.ctMatlabInterface.thermo_cv_mole|
% C++ Signature: double thermo_cv_mole(int n)

thermo_cv_moleDefinition = addFunction(libDef, ...
    "double thermo_cv_mole(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_cv_mole", ...
    "Description", "clib.ctMatlabInterface.thermo_cv_mole Representation of C++ function thermo_cv_mole."); % Modify help description values as needed.
defineArgument(thermo_cv_moleDefinition, "n", "int32");
defineOutput(thermo_cv_moleDefinition, "RetVal", "double");
validate(thermo_cv_moleDefinition);

%% C++ function |thermo_pressure| with MATLAB name |clib.ctMatlabInterface.thermo_pressure|
% C++ Signature: double thermo_pressure(int n)

thermo_pressureDefinition = addFunction(libDef, ...
    "double thermo_pressure(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_pressure", ...
    "Description", "clib.ctMatlabInterface.thermo_pressure Representation of C++ function thermo_pressure."); % Modify help description values as needed.
defineArgument(thermo_pressureDefinition, "n", "int32");
defineOutput(thermo_pressureDefinition, "RetVal", "double");
validate(thermo_pressureDefinition);

%% C++ function |thermo_setPressure| with MATLAB name |clib.ctMatlabInterface.thermo_setPressure|
% C++ Signature: int thermo_setPressure(int n,double p)

thermo_setPressureDefinition = addFunction(libDef, ...
    "int thermo_setPressure(int n,double p)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_setPressure", ...
    "Description", "clib.ctMatlabInterface.thermo_setPressure Representation of C++ function thermo_setPressure."); % Modify help description values as needed.
defineArgument(thermo_setPressureDefinition, "n", "int32");
defineArgument(thermo_setPressureDefinition, "p", "double");
defineOutput(thermo_setPressureDefinition, "RetVal", "int32");
validate(thermo_setPressureDefinition);

%% C++ function |thermo_enthalpy_mass| with MATLAB name |clib.ctMatlabInterface.thermo_enthalpy_mass|
% C++ Signature: double thermo_enthalpy_mass(int n)

thermo_enthalpy_massDefinition = addFunction(libDef, ...
    "double thermo_enthalpy_mass(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_enthalpy_mass", ...
    "Description", "clib.ctMatlabInterface.thermo_enthalpy_mass Representation of C++ function thermo_enthalpy_mass."); % Modify help description values as needed.
defineArgument(thermo_enthalpy_massDefinition, "n", "int32");
defineOutput(thermo_enthalpy_massDefinition, "RetVal", "double");
validate(thermo_enthalpy_massDefinition);

%% C++ function |thermo_intEnergy_mass| with MATLAB name |clib.ctMatlabInterface.thermo_intEnergy_mass|
% C++ Signature: double thermo_intEnergy_mass(int n)

thermo_intEnergy_massDefinition = addFunction(libDef, ...
    "double thermo_intEnergy_mass(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_intEnergy_mass", ...
    "Description", "clib.ctMatlabInterface.thermo_intEnergy_mass Representation of C++ function thermo_intEnergy_mass."); % Modify help description values as needed.
defineArgument(thermo_intEnergy_massDefinition, "n", "int32");
defineOutput(thermo_intEnergy_massDefinition, "RetVal", "double");
validate(thermo_intEnergy_massDefinition);

%% C++ function |thermo_entropy_mass| with MATLAB name |clib.ctMatlabInterface.thermo_entropy_mass|
% C++ Signature: double thermo_entropy_mass(int n)

thermo_entropy_massDefinition = addFunction(libDef, ...
    "double thermo_entropy_mass(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_entropy_mass", ...
    "Description", "clib.ctMatlabInterface.thermo_entropy_mass Representation of C++ function thermo_entropy_mass."); % Modify help description values as needed.
defineArgument(thermo_entropy_massDefinition, "n", "int32");
defineOutput(thermo_entropy_massDefinition, "RetVal", "double");
validate(thermo_entropy_massDefinition);

%% C++ function |thermo_gibbs_mass| with MATLAB name |clib.ctMatlabInterface.thermo_gibbs_mass|
% C++ Signature: double thermo_gibbs_mass(int n)

thermo_gibbs_massDefinition = addFunction(libDef, ...
    "double thermo_gibbs_mass(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_gibbs_mass", ...
    "Description", "clib.ctMatlabInterface.thermo_gibbs_mass Representation of C++ function thermo_gibbs_mass."); % Modify help description values as needed.
defineArgument(thermo_gibbs_massDefinition, "n", "int32");
defineOutput(thermo_gibbs_massDefinition, "RetVal", "double");
validate(thermo_gibbs_massDefinition);

%% C++ function |thermo_cp_mass| with MATLAB name |clib.ctMatlabInterface.thermo_cp_mass|
% C++ Signature: double thermo_cp_mass(int n)

thermo_cp_massDefinition = addFunction(libDef, ...
    "double thermo_cp_mass(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_cp_mass", ...
    "Description", "clib.ctMatlabInterface.thermo_cp_mass Representation of C++ function thermo_cp_mass."); % Modify help description values as needed.
defineArgument(thermo_cp_massDefinition, "n", "int32");
defineOutput(thermo_cp_massDefinition, "RetVal", "double");
validate(thermo_cp_massDefinition);

%% C++ function |thermo_cv_mass| with MATLAB name |clib.ctMatlabInterface.thermo_cv_mass|
% C++ Signature: double thermo_cv_mass(int n)

thermo_cv_massDefinition = addFunction(libDef, ...
    "double thermo_cv_mass(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_cv_mass", ...
    "Description", "clib.ctMatlabInterface.thermo_cv_mass Representation of C++ function thermo_cv_mass."); % Modify help description values as needed.
defineArgument(thermo_cv_massDefinition, "n", "int32");
defineOutput(thermo_cv_massDefinition, "RetVal", "double");
validate(thermo_cv_massDefinition);

%% C++ function |thermo_electricPotential| with MATLAB name |clib.ctMatlabInterface.thermo_electricPotential|
% C++ Signature: double thermo_electricPotential(int n)

thermo_electricPotentialDefinition = addFunction(libDef, ...
    "double thermo_electricPotential(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_electricPotential", ...
    "Description", "clib.ctMatlabInterface.thermo_electricPotential Representation of C++ function thermo_electricPotential."); % Modify help description values as needed.
defineArgument(thermo_electricPotentialDefinition, "n", "int32");
defineOutput(thermo_electricPotentialDefinition, "RetVal", "double");
validate(thermo_electricPotentialDefinition);

%% C++ function |thermo_thermalExpansionCoeff| with MATLAB name |clib.ctMatlabInterface.thermo_thermalExpansionCoeff|
% C++ Signature: double thermo_thermalExpansionCoeff(int n)

thermo_thermalExpansionCoeffDefinition = addFunction(libDef, ...
    "double thermo_thermalExpansionCoeff(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_thermalExpansionCoeff", ...
    "Description", "clib.ctMatlabInterface.thermo_thermalExpansionCoeff Representation of C++ function thermo_thermalExpansionCoeff."); % Modify help description values as needed.
defineArgument(thermo_thermalExpansionCoeffDefinition, "n", "int32");
defineOutput(thermo_thermalExpansionCoeffDefinition, "RetVal", "double");
validate(thermo_thermalExpansionCoeffDefinition);

%% C++ function |thermo_isothermalCompressibility| with MATLAB name |clib.ctMatlabInterface.thermo_isothermalCompressibility|
% C++ Signature: double thermo_isothermalCompressibility(int n)

thermo_isothermalCompressibilityDefinition = addFunction(libDef, ...
    "double thermo_isothermalCompressibility(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_isothermalCompressibility", ...
    "Description", "clib.ctMatlabInterface.thermo_isothermalCompressibility Representation of C++ function thermo_isothermalCompressibility."); % Modify help description values as needed.
defineArgument(thermo_isothermalCompressibilityDefinition, "n", "int32");
defineOutput(thermo_isothermalCompressibilityDefinition, "RetVal", "double");
validate(thermo_isothermalCompressibilityDefinition);

%% C++ function |thermo_chemPotentials| with MATLAB name |clib.ctMatlabInterface.thermo_chemPotentials|
% C++ Signature: int thermo_chemPotentials(int n,size_t lenm,double * murt)

thermo_chemPotentialsDefinition = addFunction(libDef, ...
    "int thermo_chemPotentials(int n,size_t lenm,double * murt)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_chemPotentials", ...
    "Description", "clib.ctMatlabInterface.thermo_chemPotentials Representation of C++ function thermo_chemPotentials."); % Modify help description values as needed.
defineArgument(thermo_chemPotentialsDefinition, "n", "int32");
defineArgument(thermo_chemPotentialsDefinition, "lenm", "uint64");
defineArgument(thermo_chemPotentialsDefinition, "murt", "clib.array.ctMatlabInterface.Double", "input", "lenm"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_chemPotentialsDefinition, "RetVal", "int32");
validate(thermo_chemPotentialsDefinition);

%% C++ function |thermo_electrochemPotentials| with MATLAB name |clib.ctMatlabInterface.thermo_electrochemPotentials|
% C++ Signature: int thermo_electrochemPotentials(int n,size_t lenm,double * emu)

thermo_electrochemPotentialsDefinition = addFunction(libDef, ...
    "int thermo_electrochemPotentials(int n,size_t lenm,double * emu)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_electrochemPotentials", ...
    "Description", "clib.ctMatlabInterface.thermo_electrochemPotentials Representation of C++ function thermo_electrochemPotentials."); % Modify help description values as needed.
defineArgument(thermo_electrochemPotentialsDefinition, "n", "int32");
defineArgument(thermo_electrochemPotentialsDefinition, "lenm", "uint64");
defineArgument(thermo_electrochemPotentialsDefinition, "emu", "clib.array.ctMatlabInterface.Double", "input", "lenm"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_electrochemPotentialsDefinition, "RetVal", "int32");
validate(thermo_electrochemPotentialsDefinition);

%% C++ function |thermo_getEnthalpies_RT| with MATLAB name |clib.ctMatlabInterface.thermo_getEnthalpies_RT|
% C++ Signature: int thermo_getEnthalpies_RT(int n,size_t lenm,double * h_rt)

thermo_getEnthalpies_RTDefinition = addFunction(libDef, ...
    "int thermo_getEnthalpies_RT(int n,size_t lenm,double * h_rt)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_getEnthalpies_RT", ...
    "Description", "clib.ctMatlabInterface.thermo_getEnthalpies_RT Representation of C++ function thermo_getEnthalpies_RT."); % Modify help description values as needed.
defineArgument(thermo_getEnthalpies_RTDefinition, "n", "int32");
defineArgument(thermo_getEnthalpies_RTDefinition, "lenm", "uint64");
defineArgument(thermo_getEnthalpies_RTDefinition, "h_rt", "clib.array.ctMatlabInterface.Double", "input", "lenm"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_getEnthalpies_RTDefinition, "RetVal", "int32");
validate(thermo_getEnthalpies_RTDefinition);

%% C++ function |thermo_getEntropies_R| with MATLAB name |clib.ctMatlabInterface.thermo_getEntropies_R|
% C++ Signature: int thermo_getEntropies_R(int n,size_t lenm,double * s_r)

thermo_getEntropies_RDefinition = addFunction(libDef, ...
    "int thermo_getEntropies_R(int n,size_t lenm,double * s_r)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_getEntropies_R", ...
    "Description", "clib.ctMatlabInterface.thermo_getEntropies_R Representation of C++ function thermo_getEntropies_R."); % Modify help description values as needed.
defineArgument(thermo_getEntropies_RDefinition, "n", "int32");
defineArgument(thermo_getEntropies_RDefinition, "lenm", "uint64");
defineArgument(thermo_getEntropies_RDefinition, "s_r", "clib.array.ctMatlabInterface.Double", "input", "lenm"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_getEntropies_RDefinition, "RetVal", "int32");
validate(thermo_getEntropies_RDefinition);

%% C++ function |thermo_getCp_R| with MATLAB name |clib.ctMatlabInterface.thermo_getCp_R|
% C++ Signature: int thermo_getCp_R(int n,size_t lenm,double * cp_r)

thermo_getCp_RDefinition = addFunction(libDef, ...
    "int thermo_getCp_R(int n,size_t lenm,double * cp_r)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_getCp_R", ...
    "Description", "clib.ctMatlabInterface.thermo_getCp_R Representation of C++ function thermo_getCp_R."); % Modify help description values as needed.
defineArgument(thermo_getCp_RDefinition, "n", "int32");
defineArgument(thermo_getCp_RDefinition, "lenm", "uint64");
defineArgument(thermo_getCp_RDefinition, "cp_r", "clib.array.ctMatlabInterface.Double", "input", "lenm"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_getCp_RDefinition, "RetVal", "int32");
validate(thermo_getCp_RDefinition);

%% C++ function |thermo_setElectricPotential| with MATLAB name |clib.ctMatlabInterface.thermo_setElectricPotential|
% C++ Signature: int thermo_setElectricPotential(int n,double v)

thermo_setElectricPotentialDefinition = addFunction(libDef, ...
    "int thermo_setElectricPotential(int n,double v)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_setElectricPotential", ...
    "Description", "clib.ctMatlabInterface.thermo_setElectricPotential Representation of C++ function thermo_setElectricPotential."); % Modify help description values as needed.
defineArgument(thermo_setElectricPotentialDefinition, "n", "int32");
defineArgument(thermo_setElectricPotentialDefinition, "v", "double");
defineOutput(thermo_setElectricPotentialDefinition, "RetVal", "int32");
validate(thermo_setElectricPotentialDefinition);

%% C++ function |thermo_getPartialMolarEnthalpies| with MATLAB name |clib.ctMatlabInterface.thermo_getPartialMolarEnthalpies|
% C++ Signature: int thermo_getPartialMolarEnthalpies(int n,size_t lenm,double * pmh)

thermo_getPartialMolarEnthalpiesDefinition = addFunction(libDef, ...
    "int thermo_getPartialMolarEnthalpies(int n,size_t lenm,double * pmh)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_getPartialMolarEnthalpies", ...
    "Description", "clib.ctMatlabInterface.thermo_getPartialMolarEnthalpies Representation of C++ function thermo_getPartialMolarEnthalpies."); % Modify help description values as needed.
defineArgument(thermo_getPartialMolarEnthalpiesDefinition, "n", "int32");
defineArgument(thermo_getPartialMolarEnthalpiesDefinition, "lenm", "uint64");
defineArgument(thermo_getPartialMolarEnthalpiesDefinition, "pmh", "clib.array.ctMatlabInterface.Double", "input", "lenm"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_getPartialMolarEnthalpiesDefinition, "RetVal", "int32");
validate(thermo_getPartialMolarEnthalpiesDefinition);

%% C++ function |thermo_getPartialMolarEntropies| with MATLAB name |clib.ctMatlabInterface.thermo_getPartialMolarEntropies|
% C++ Signature: int thermo_getPartialMolarEntropies(int n,size_t lenm,double * pms)

thermo_getPartialMolarEntropiesDefinition = addFunction(libDef, ...
    "int thermo_getPartialMolarEntropies(int n,size_t lenm,double * pms)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_getPartialMolarEntropies", ...
    "Description", "clib.ctMatlabInterface.thermo_getPartialMolarEntropies Representation of C++ function thermo_getPartialMolarEntropies."); % Modify help description values as needed.
defineArgument(thermo_getPartialMolarEntropiesDefinition, "n", "int32");
defineArgument(thermo_getPartialMolarEntropiesDefinition, "lenm", "uint64");
defineArgument(thermo_getPartialMolarEntropiesDefinition, "pms", "clib.array.ctMatlabInterface.Double", "input", "lenm"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_getPartialMolarEntropiesDefinition, "RetVal", "int32");
validate(thermo_getPartialMolarEntropiesDefinition);

%% C++ function |thermo_getPartialMolarIntEnergies| with MATLAB name |clib.ctMatlabInterface.thermo_getPartialMolarIntEnergies|
% C++ Signature: int thermo_getPartialMolarIntEnergies(int n,size_t lenm,double * pmu)

thermo_getPartialMolarIntEnergiesDefinition = addFunction(libDef, ...
    "int thermo_getPartialMolarIntEnergies(int n,size_t lenm,double * pmu)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_getPartialMolarIntEnergies", ...
    "Description", "clib.ctMatlabInterface.thermo_getPartialMolarIntEnergies Representation of C++ function thermo_getPartialMolarIntEnergies."); % Modify help description values as needed.
defineArgument(thermo_getPartialMolarIntEnergiesDefinition, "n", "int32");
defineArgument(thermo_getPartialMolarIntEnergiesDefinition, "lenm", "uint64");
defineArgument(thermo_getPartialMolarIntEnergiesDefinition, "pmu", "clib.array.ctMatlabInterface.Double", "input", "lenm"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_getPartialMolarIntEnergiesDefinition, "RetVal", "int32");
validate(thermo_getPartialMolarIntEnergiesDefinition);

%% C++ function |thermo_getPartialMolarCp| with MATLAB name |clib.ctMatlabInterface.thermo_getPartialMolarCp|
% C++ Signature: int thermo_getPartialMolarCp(int n,size_t lenm,double * pmcp)

thermo_getPartialMolarCpDefinition = addFunction(libDef, ...
    "int thermo_getPartialMolarCp(int n,size_t lenm,double * pmcp)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_getPartialMolarCp", ...
    "Description", "clib.ctMatlabInterface.thermo_getPartialMolarCp Representation of C++ function thermo_getPartialMolarCp."); % Modify help description values as needed.
defineArgument(thermo_getPartialMolarCpDefinition, "n", "int32");
defineArgument(thermo_getPartialMolarCpDefinition, "lenm", "uint64");
defineArgument(thermo_getPartialMolarCpDefinition, "pmcp", "clib.array.ctMatlabInterface.Double", "input", "lenm"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_getPartialMolarCpDefinition, "RetVal", "int32");
validate(thermo_getPartialMolarCpDefinition);

%% C++ function |thermo_getPartialMolarVolumes| with MATLAB name |clib.ctMatlabInterface.thermo_getPartialMolarVolumes|
% C++ Signature: int thermo_getPartialMolarVolumes(int n,size_t lenm,double * pmv)

thermo_getPartialMolarVolumesDefinition = addFunction(libDef, ...
    "int thermo_getPartialMolarVolumes(int n,size_t lenm,double * pmv)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_getPartialMolarVolumes", ...
    "Description", "clib.ctMatlabInterface.thermo_getPartialMolarVolumes Representation of C++ function thermo_getPartialMolarVolumes."); % Modify help description values as needed.
defineArgument(thermo_getPartialMolarVolumesDefinition, "n", "int32");
defineArgument(thermo_getPartialMolarVolumesDefinition, "lenm", "uint64");
defineArgument(thermo_getPartialMolarVolumesDefinition, "pmv", "clib.array.ctMatlabInterface.Double", "input", "lenm"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_getPartialMolarVolumesDefinition, "RetVal", "int32");
validate(thermo_getPartialMolarVolumesDefinition);

%% C++ function |thermo_set_TP| with MATLAB name |clib.ctMatlabInterface.thermo_set_TP|
% C++ Signature: int thermo_set_TP(int n,double * vals)

thermo_set_TPDefinition = addFunction(libDef, ...
    "int thermo_set_TP(int n,double * vals)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_set_TP", ...
    "Description", "clib.ctMatlabInterface.thermo_set_TP Representation of C++ function thermo_set_TP."); % Modify help description values as needed.
defineArgument(thermo_set_TPDefinition, "n", "int32");
defineArgument(thermo_set_TPDefinition, "vals", "clib.array.ctMatlabInterface.Double", "input", 2); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_set_TPDefinition, "RetVal", "int32");
validate(thermo_set_TPDefinition);

%% C++ function |thermo_set_TD| with MATLAB name |clib.ctMatlabInterface.thermo_set_TD|
% C++ Signature: int thermo_set_TD(int n,double * vals)

thermo_set_TDDefinition = addFunction(libDef, ...
    "int thermo_set_TD(int n,double * vals)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_set_TD", ...
    "Description", "clib.ctMatlabInterface.thermo_set_TD Representation of C++ function thermo_set_TD."); % Modify help description values as needed.
defineArgument(thermo_set_TDDefinition, "n", "int32");
defineArgument(thermo_set_TDDefinition, "vals", "clib.array.ctMatlabInterface.Double", "input", 2); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_set_TDDefinition, "RetVal", "int32");
validate(thermo_set_TDDefinition);

%% C++ function |thermo_set_DP| with MATLAB name |clib.ctMatlabInterface.thermo_set_DP|
% C++ Signature: int thermo_set_DP(int n,double * vals)

thermo_set_DPDefinition = addFunction(libDef, ...
    "int thermo_set_DP(int n,double * vals)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_set_DP", ...
    "Description", "clib.ctMatlabInterface.thermo_set_DP Representation of C++ function thermo_set_DP."); % Modify help description values as needed.
defineArgument(thermo_set_DPDefinition, "n", "int32");
defineArgument(thermo_set_DPDefinition, "vals", "clib.array.ctMatlabInterface.Double", "input", 2); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_set_DPDefinition, "RetVal", "int32");
validate(thermo_set_DPDefinition);

%% C++ function |thermo_set_HP| with MATLAB name |clib.ctMatlabInterface.thermo_set_HP|
% C++ Signature: int thermo_set_HP(int n,double * vals)

thermo_set_HPDefinition = addFunction(libDef, ...
    "int thermo_set_HP(int n,double * vals)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_set_HP", ...
    "Description", "clib.ctMatlabInterface.thermo_set_HP Representation of C++ function thermo_set_HP."); % Modify help description values as needed.
defineArgument(thermo_set_HPDefinition, "n", "int32");
defineArgument(thermo_set_HPDefinition, "vals", "clib.array.ctMatlabInterface.Double", "input", 2); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_set_HPDefinition, "RetVal", "int32");
validate(thermo_set_HPDefinition);

%% C++ function |thermo_set_UV| with MATLAB name |clib.ctMatlabInterface.thermo_set_UV|
% C++ Signature: int thermo_set_UV(int n,double * vals)

thermo_set_UVDefinition = addFunction(libDef, ...
    "int thermo_set_UV(int n,double * vals)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_set_UV", ...
    "Description", "clib.ctMatlabInterface.thermo_set_UV Representation of C++ function thermo_set_UV."); % Modify help description values as needed.
defineArgument(thermo_set_UVDefinition, "n", "int32");
defineArgument(thermo_set_UVDefinition, "vals", "clib.array.ctMatlabInterface.Double", "input", 2); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_set_UVDefinition, "RetVal", "int32");
validate(thermo_set_UVDefinition);

%% C++ function |thermo_set_SV| with MATLAB name |clib.ctMatlabInterface.thermo_set_SV|
% C++ Signature: int thermo_set_SV(int n,double * vals)

thermo_set_SVDefinition = addFunction(libDef, ...
    "int thermo_set_SV(int n,double * vals)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_set_SV", ...
    "Description", "clib.ctMatlabInterface.thermo_set_SV Representation of C++ function thermo_set_SV."); % Modify help description values as needed.
defineArgument(thermo_set_SVDefinition, "n", "int32");
defineArgument(thermo_set_SVDefinition, "vals", "clib.array.ctMatlabInterface.Double", "input", 2); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_set_SVDefinition, "RetVal", "int32");
validate(thermo_set_SVDefinition);

%% C++ function |thermo_set_SP| with MATLAB name |clib.ctMatlabInterface.thermo_set_SP|
% C++ Signature: int thermo_set_SP(int n,double * vals)

thermo_set_SPDefinition = addFunction(libDef, ...
    "int thermo_set_SP(int n,double * vals)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_set_SP", ...
    "Description", "clib.ctMatlabInterface.thermo_set_SP Representation of C++ function thermo_set_SP."); % Modify help description values as needed.
defineArgument(thermo_set_SPDefinition, "n", "int32");
defineArgument(thermo_set_SPDefinition, "vals", "clib.array.ctMatlabInterface.Double", "input", 2); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_set_SPDefinition, "RetVal", "int32");
validate(thermo_set_SPDefinition);

%% C++ function |thermo_set_ST| with MATLAB name |clib.ctMatlabInterface.thermo_set_ST|
% C++ Signature: int thermo_set_ST(int n,double * vals)

thermo_set_STDefinition = addFunction(libDef, ...
    "int thermo_set_ST(int n,double * vals)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_set_ST", ...
    "Description", "clib.ctMatlabInterface.thermo_set_ST Representation of C++ function thermo_set_ST."); % Modify help description values as needed.
defineArgument(thermo_set_STDefinition, "n", "int32");
defineArgument(thermo_set_STDefinition, "vals", "clib.array.ctMatlabInterface.Double", "input", 2); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_set_STDefinition, "RetVal", "int32");
validate(thermo_set_STDefinition);

%% C++ function |thermo_set_TV| with MATLAB name |clib.ctMatlabInterface.thermo_set_TV|
% C++ Signature: int thermo_set_TV(int n,double * vals)

thermo_set_TVDefinition = addFunction(libDef, ...
    "int thermo_set_TV(int n,double * vals)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_set_TV", ...
    "Description", "clib.ctMatlabInterface.thermo_set_TV Representation of C++ function thermo_set_TV."); % Modify help description values as needed.
defineArgument(thermo_set_TVDefinition, "n", "int32");
defineArgument(thermo_set_TVDefinition, "vals", "clib.array.ctMatlabInterface.Double", "input", 2); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_set_TVDefinition, "RetVal", "int32");
validate(thermo_set_TVDefinition);

%% C++ function |thermo_set_PV| with MATLAB name |clib.ctMatlabInterface.thermo_set_PV|
% C++ Signature: int thermo_set_PV(int n,double * vals)

thermo_set_PVDefinition = addFunction(libDef, ...
    "int thermo_set_PV(int n,double * vals)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_set_PV", ...
    "Description", "clib.ctMatlabInterface.thermo_set_PV Representation of C++ function thermo_set_PV."); % Modify help description values as needed.
defineArgument(thermo_set_PVDefinition, "n", "int32");
defineArgument(thermo_set_PVDefinition, "vals", "clib.array.ctMatlabInterface.Double", "input", 2); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_set_PVDefinition, "RetVal", "int32");
validate(thermo_set_PVDefinition);

%% C++ function |thermo_set_UP| with MATLAB name |clib.ctMatlabInterface.thermo_set_UP|
% C++ Signature: int thermo_set_UP(int n,double * vals)

thermo_set_UPDefinition = addFunction(libDef, ...
    "int thermo_set_UP(int n,double * vals)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_set_UP", ...
    "Description", "clib.ctMatlabInterface.thermo_set_UP Representation of C++ function thermo_set_UP."); % Modify help description values as needed.
defineArgument(thermo_set_UPDefinition, "n", "int32");
defineArgument(thermo_set_UPDefinition, "vals", "clib.array.ctMatlabInterface.Double", "input", 2); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_set_UPDefinition, "RetVal", "int32");
validate(thermo_set_UPDefinition);

%% C++ function |thermo_set_VH| with MATLAB name |clib.ctMatlabInterface.thermo_set_VH|
% C++ Signature: int thermo_set_VH(int n,double * vals)

thermo_set_VHDefinition = addFunction(libDef, ...
    "int thermo_set_VH(int n,double * vals)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_set_VH", ...
    "Description", "clib.ctMatlabInterface.thermo_set_VH Representation of C++ function thermo_set_VH."); % Modify help description values as needed.
defineArgument(thermo_set_VHDefinition, "n", "int32");
defineArgument(thermo_set_VHDefinition, "vals", "clib.array.ctMatlabInterface.Double", "input", 2); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_set_VHDefinition, "RetVal", "int32");
validate(thermo_set_VHDefinition);

%% C++ function |thermo_set_TH| with MATLAB name |clib.ctMatlabInterface.thermo_set_TH|
% C++ Signature: int thermo_set_TH(int n,double * vals)

thermo_set_THDefinition = addFunction(libDef, ...
    "int thermo_set_TH(int n,double * vals)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_set_TH", ...
    "Description", "clib.ctMatlabInterface.thermo_set_TH Representation of C++ function thermo_set_TH."); % Modify help description values as needed.
defineArgument(thermo_set_THDefinition, "n", "int32");
defineArgument(thermo_set_THDefinition, "vals", "clib.array.ctMatlabInterface.Double", "input", 2); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_set_THDefinition, "RetVal", "int32");
validate(thermo_set_THDefinition);

%% C++ function |thermo_set_SH| with MATLAB name |clib.ctMatlabInterface.thermo_set_SH|
% C++ Signature: int thermo_set_SH(int n,double * vals)

thermo_set_SHDefinition = addFunction(libDef, ...
    "int thermo_set_SH(int n,double * vals)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_set_SH", ...
    "Description", "clib.ctMatlabInterface.thermo_set_SH Representation of C++ function thermo_set_SH."); % Modify help description values as needed.
defineArgument(thermo_set_SHDefinition, "n", "int32");
defineArgument(thermo_set_SHDefinition, "vals", "clib.array.ctMatlabInterface.Double", "input", 2); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(thermo_set_SHDefinition, "RetVal", "int32");
validate(thermo_set_SHDefinition);

%% C++ function |thermo_equilibrate| with MATLAB name |clib.ctMatlabInterface.thermo_equilibrate|
% C++ Signature: int thermo_equilibrate(int n,char const * XY,char const * solver,double rtol,int maxsteps,int maxiter,int loglevel)

thermo_equilibrateDefinition = addFunction(libDef, ...
    "int thermo_equilibrate(int n,char const * XY,char const * solver,double rtol,int maxsteps,int maxiter,int loglevel)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_equilibrate", ...
    "Description", "clib.ctMatlabInterface.thermo_equilibrate Representation of C++ function thermo_equilibrate."); % Modify help description values as needed.
defineArgument(thermo_equilibrateDefinition, "n", "int32");
defineArgument(thermo_equilibrateDefinition, "XY", "string", "input", "nullTerminated");
defineArgument(thermo_equilibrateDefinition, "solver", "string", "input", "nullTerminated");
defineArgument(thermo_equilibrateDefinition, "rtol", "double");
defineArgument(thermo_equilibrateDefinition, "maxsteps", "int32");
defineArgument(thermo_equilibrateDefinition, "maxiter", "int32");
defineArgument(thermo_equilibrateDefinition, "loglevel", "int32");
defineOutput(thermo_equilibrateDefinition, "RetVal", "int32");
validate(thermo_equilibrateDefinition);

%% C++ function |thermo_critTemperature| with MATLAB name |clib.ctMatlabInterface.thermo_critTemperature|
% C++ Signature: double thermo_critTemperature(int n)

thermo_critTemperatureDefinition = addFunction(libDef, ...
    "double thermo_critTemperature(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_critTemperature", ...
    "Description", "clib.ctMatlabInterface.thermo_critTemperature Representation of C++ function thermo_critTemperature."); % Modify help description values as needed.
defineArgument(thermo_critTemperatureDefinition, "n", "int32");
defineOutput(thermo_critTemperatureDefinition, "RetVal", "double");
validate(thermo_critTemperatureDefinition);

%% C++ function |thermo_critPressure| with MATLAB name |clib.ctMatlabInterface.thermo_critPressure|
% C++ Signature: double thermo_critPressure(int n)

thermo_critPressureDefinition = addFunction(libDef, ...
    "double thermo_critPressure(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_critPressure", ...
    "Description", "clib.ctMatlabInterface.thermo_critPressure Representation of C++ function thermo_critPressure."); % Modify help description values as needed.
defineArgument(thermo_critPressureDefinition, "n", "int32");
defineOutput(thermo_critPressureDefinition, "RetVal", "double");
validate(thermo_critPressureDefinition);

%% C++ function |thermo_critDensity| with MATLAB name |clib.ctMatlabInterface.thermo_critDensity|
% C++ Signature: double thermo_critDensity(int n)

thermo_critDensityDefinition = addFunction(libDef, ...
    "double thermo_critDensity(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_critDensity", ...
    "Description", "clib.ctMatlabInterface.thermo_critDensity Representation of C++ function thermo_critDensity."); % Modify help description values as needed.
defineArgument(thermo_critDensityDefinition, "n", "int32");
defineOutput(thermo_critDensityDefinition, "RetVal", "double");
validate(thermo_critDensityDefinition);

%% C++ function |thermo_vaporFraction| with MATLAB name |clib.ctMatlabInterface.thermo_vaporFraction|
% C++ Signature: double thermo_vaporFraction(int n)

thermo_vaporFractionDefinition = addFunction(libDef, ...
    "double thermo_vaporFraction(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_vaporFraction", ...
    "Description", "clib.ctMatlabInterface.thermo_vaporFraction Representation of C++ function thermo_vaporFraction."); % Modify help description values as needed.
defineArgument(thermo_vaporFractionDefinition, "n", "int32");
defineOutput(thermo_vaporFractionDefinition, "RetVal", "double");
validate(thermo_vaporFractionDefinition);

%% C++ function |thermo_satTemperature| with MATLAB name |clib.ctMatlabInterface.thermo_satTemperature|
% C++ Signature: double thermo_satTemperature(int n,double p)

thermo_satTemperatureDefinition = addFunction(libDef, ...
    "double thermo_satTemperature(int n,double p)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_satTemperature", ...
    "Description", "clib.ctMatlabInterface.thermo_satTemperature Representation of C++ function thermo_satTemperature."); % Modify help description values as needed.
defineArgument(thermo_satTemperatureDefinition, "n", "int32");
defineArgument(thermo_satTemperatureDefinition, "p", "double");
defineOutput(thermo_satTemperatureDefinition, "RetVal", "double");
validate(thermo_satTemperatureDefinition);

%% C++ function |thermo_satPressure| with MATLAB name |clib.ctMatlabInterface.thermo_satPressure|
% C++ Signature: double thermo_satPressure(int n,double t)

thermo_satPressureDefinition = addFunction(libDef, ...
    "double thermo_satPressure(int n,double t)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_satPressure", ...
    "Description", "clib.ctMatlabInterface.thermo_satPressure Representation of C++ function thermo_satPressure."); % Modify help description values as needed.
defineArgument(thermo_satPressureDefinition, "n", "int32");
defineArgument(thermo_satPressureDefinition, "t", "double");
defineOutput(thermo_satPressureDefinition, "RetVal", "double");
validate(thermo_satPressureDefinition);

%% C++ function |thermo_setState_Psat| with MATLAB name |clib.ctMatlabInterface.thermo_setState_Psat|
% C++ Signature: int thermo_setState_Psat(int n,double p,double x)

thermo_setState_PsatDefinition = addFunction(libDef, ...
    "int thermo_setState_Psat(int n,double p,double x)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_setState_Psat", ...
    "Description", "clib.ctMatlabInterface.thermo_setState_Psat Representation of C++ function thermo_setState_Psat."); % Modify help description values as needed.
defineArgument(thermo_setState_PsatDefinition, "n", "int32");
defineArgument(thermo_setState_PsatDefinition, "p", "double");
defineArgument(thermo_setState_PsatDefinition, "x", "double");
defineOutput(thermo_setState_PsatDefinition, "RetVal", "int32");
validate(thermo_setState_PsatDefinition);

%% C++ function |thermo_setState_Tsat| with MATLAB name |clib.ctMatlabInterface.thermo_setState_Tsat|
% C++ Signature: int thermo_setState_Tsat(int n,double t,double x)

thermo_setState_TsatDefinition = addFunction(libDef, ...
    "int thermo_setState_Tsat(int n,double t,double x)", ...
    "MATLABName", "clib.ctMatlabInterface.thermo_setState_Tsat", ...
    "Description", "clib.ctMatlabInterface.thermo_setState_Tsat Representation of C++ function thermo_setState_Tsat."); % Modify help description values as needed.
defineArgument(thermo_setState_TsatDefinition, "n", "int32");
defineArgument(thermo_setState_TsatDefinition, "t", "double");
defineArgument(thermo_setState_TsatDefinition, "x", "double");
defineOutput(thermo_setState_TsatDefinition, "RetVal", "int32");
validate(thermo_setState_TsatDefinition);

%% C++ function |kin_parent| with MATLAB name |clib.ctMatlabInterface.kin_parent|
% C++ Signature: int kin_parent(int n)

kin_parentDefinition = addFunction(libDef, ...
    "int kin_parent(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_parent", ...
    "Description", "clib.ctMatlabInterface.kin_parent Representation of C++ function kin_parent."); % Modify help description values as needed.
defineArgument(kin_parentDefinition, "n", "int32");
defineOutput(kin_parentDefinition, "RetVal", "int32");
validate(kin_parentDefinition);

%% C++ function |kin_del| with MATLAB name |clib.ctMatlabInterface.kin_del|
% C++ Signature: int kin_del(int n)

kin_delDefinition = addFunction(libDef, ...
    "int kin_del(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_del", ...
    "Description", "clib.ctMatlabInterface.kin_del Representation of C++ function kin_del." + newline + ...
    "no-op; object is managed by Solution"); % Modify help description values as needed.
defineArgument(kin_delDefinition, "n", "int32");
defineOutput(kin_delDefinition, "RetVal", "int32");
validate(kin_delDefinition);

%% C++ function |kin_nSpecies| with MATLAB name |clib.ctMatlabInterface.kin_nSpecies|
% C++ Signature: size_t kin_nSpecies(int n)

kin_nSpeciesDefinition = addFunction(libDef, ...
    "size_t kin_nSpecies(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_nSpecies", ...
    "Description", "clib.ctMatlabInterface.kin_nSpecies Representation of C++ function kin_nSpecies."); % Modify help description values as needed.
defineArgument(kin_nSpeciesDefinition, "n", "int32");
defineOutput(kin_nSpeciesDefinition, "RetVal", "uint64");
validate(kin_nSpeciesDefinition);

%% C++ function |kin_nReactions| with MATLAB name |clib.ctMatlabInterface.kin_nReactions|
% C++ Signature: size_t kin_nReactions(int n)

kin_nReactionsDefinition = addFunction(libDef, ...
    "size_t kin_nReactions(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_nReactions", ...
    "Description", "clib.ctMatlabInterface.kin_nReactions Representation of C++ function kin_nReactions."); % Modify help description values as needed.
defineArgument(kin_nReactionsDefinition, "n", "int32");
defineOutput(kin_nReactionsDefinition, "RetVal", "uint64");
validate(kin_nReactionsDefinition);

%% C++ function |kin_nPhases| with MATLAB name |clib.ctMatlabInterface.kin_nPhases|
% C++ Signature: size_t kin_nPhases(int n)

kin_nPhasesDefinition = addFunction(libDef, ...
    "size_t kin_nPhases(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_nPhases", ...
    "Description", "clib.ctMatlabInterface.kin_nPhases Representation of C++ function kin_nPhases."); % Modify help description values as needed.
defineArgument(kin_nPhasesDefinition, "n", "int32");
defineOutput(kin_nPhasesDefinition, "RetVal", "uint64");
validate(kin_nPhasesDefinition);

%% C++ function |kin_phaseIndex| with MATLAB name |clib.ctMatlabInterface.kin_phaseIndex|
% C++ Signature: size_t kin_phaseIndex(int n,char const * ph)

kin_phaseIndexDefinition = addFunction(libDef, ...
    "size_t kin_phaseIndex(int n,char const * ph)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_phaseIndex", ...
    "Description", "clib.ctMatlabInterface.kin_phaseIndex Representation of C++ function kin_phaseIndex."); % Modify help description values as needed.
defineArgument(kin_phaseIndexDefinition, "n", "int32");
defineArgument(kin_phaseIndexDefinition, "ph", "string", "input", "nullTerminated");
defineOutput(kin_phaseIndexDefinition, "RetVal", "uint64");
validate(kin_phaseIndexDefinition);

%% C++ function |kin_reactantStoichCoeff| with MATLAB name |clib.ctMatlabInterface.kin_reactantStoichCoeff|
% C++ Signature: double kin_reactantStoichCoeff(int n,int i,int k)

kin_reactantStoichCoeffDefinition = addFunction(libDef, ...
    "double kin_reactantStoichCoeff(int n,int i,int k)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_reactantStoichCoeff", ...
    "Description", "clib.ctMatlabInterface.kin_reactantStoichCoeff Representation of C++ function kin_reactantStoichCoeff."); % Modify help description values as needed.
defineArgument(kin_reactantStoichCoeffDefinition, "n", "int32");
defineArgument(kin_reactantStoichCoeffDefinition, "i", "int32");
defineArgument(kin_reactantStoichCoeffDefinition, "k", "int32");
defineOutput(kin_reactantStoichCoeffDefinition, "RetVal", "double");
validate(kin_reactantStoichCoeffDefinition);

%% C++ function |kin_productStoichCoeff| with MATLAB name |clib.ctMatlabInterface.kin_productStoichCoeff|
% C++ Signature: double kin_productStoichCoeff(int n,int i,int k)

kin_productStoichCoeffDefinition = addFunction(libDef, ...
    "double kin_productStoichCoeff(int n,int i,int k)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_productStoichCoeff", ...
    "Description", "clib.ctMatlabInterface.kin_productStoichCoeff Representation of C++ function kin_productStoichCoeff."); % Modify help description values as needed.
defineArgument(kin_productStoichCoeffDefinition, "n", "int32");
defineArgument(kin_productStoichCoeffDefinition, "i", "int32");
defineArgument(kin_productStoichCoeffDefinition, "k", "int32");
defineOutput(kin_productStoichCoeffDefinition, "RetVal", "double");
validate(kin_productStoichCoeffDefinition);

%% C++ function |kin_getReactionType| with MATLAB name |clib.ctMatlabInterface.kin_getReactionType|
% C++ Signature: int kin_getReactionType(int n,int i,size_t len,char * name)

kin_getReactionTypeDefinition = addFunction(libDef, ...
    "int kin_getReactionType(int n,int i,size_t len,char * name)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_getReactionType", ...
    "Description", "clib.ctMatlabInterface.kin_getReactionType Representation of C++ function kin_getReactionType."); % Modify help description values as needed.
defineArgument(kin_getReactionTypeDefinition, "n", "int32");
defineArgument(kin_getReactionTypeDefinition, "i", "int32");
defineArgument(kin_getReactionTypeDefinition, "len", "uint64");
defineArgument(kin_getReactionTypeDefinition, "name", "clib.array.ctMatlabInterface.Char", "input", "len"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(kin_getReactionTypeDefinition, "RetVal", "int32");
validate(kin_getReactionTypeDefinition);

%% C++ function |kin_getFwdRatesOfProgress| with MATLAB name |clib.ctMatlabInterface.kin_getFwdRatesOfProgress|
% C++ Signature: int kin_getFwdRatesOfProgress(int n,size_t len,double * fwdROP)

kin_getFwdRatesOfProgressDefinition = addFunction(libDef, ...
    "int kin_getFwdRatesOfProgress(int n,size_t len,double * fwdROP)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_getFwdRatesOfProgress", ...
    "Description", "clib.ctMatlabInterface.kin_getFwdRatesOfProgress Representation of C++ function kin_getFwdRatesOfProgress."); % Modify help description values as needed.
defineArgument(kin_getFwdRatesOfProgressDefinition, "n", "int32");
defineArgument(kin_getFwdRatesOfProgressDefinition, "len", "uint64");
defineArgument(kin_getFwdRatesOfProgressDefinition, "fwdROP", "clib.array.ctMatlabInterface.Double", "input", "len"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(kin_getFwdRatesOfProgressDefinition, "RetVal", "int32");
validate(kin_getFwdRatesOfProgressDefinition);

%% C++ function |kin_getRevRatesOfProgress| with MATLAB name |clib.ctMatlabInterface.kin_getRevRatesOfProgress|
% C++ Signature: int kin_getRevRatesOfProgress(int n,size_t len,double * revROP)

kin_getRevRatesOfProgressDefinition = addFunction(libDef, ...
    "int kin_getRevRatesOfProgress(int n,size_t len,double * revROP)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_getRevRatesOfProgress", ...
    "Description", "clib.ctMatlabInterface.kin_getRevRatesOfProgress Representation of C++ function kin_getRevRatesOfProgress."); % Modify help description values as needed.
defineArgument(kin_getRevRatesOfProgressDefinition, "n", "int32");
defineArgument(kin_getRevRatesOfProgressDefinition, "len", "uint64");
defineArgument(kin_getRevRatesOfProgressDefinition, "revROP", "clib.array.ctMatlabInterface.Double", "input", "len"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(kin_getRevRatesOfProgressDefinition, "RetVal", "int32");
validate(kin_getRevRatesOfProgressDefinition);

%% C++ function |kin_getNetRatesOfProgress| with MATLAB name |clib.ctMatlabInterface.kin_getNetRatesOfProgress|
% C++ Signature: int kin_getNetRatesOfProgress(int n,size_t len,double * netROP)

kin_getNetRatesOfProgressDefinition = addFunction(libDef, ...
    "int kin_getNetRatesOfProgress(int n,size_t len,double * netROP)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_getNetRatesOfProgress", ...
    "Description", "clib.ctMatlabInterface.kin_getNetRatesOfProgress Representation of C++ function kin_getNetRatesOfProgress."); % Modify help description values as needed.
defineArgument(kin_getNetRatesOfProgressDefinition, "n", "int32");
defineArgument(kin_getNetRatesOfProgressDefinition, "len", "uint64");
defineArgument(kin_getNetRatesOfProgressDefinition, "netROP", "clib.array.ctMatlabInterface.Double", "input", "len"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(kin_getNetRatesOfProgressDefinition, "RetVal", "int32");
validate(kin_getNetRatesOfProgressDefinition);

%% C++ function |kin_getEquilibriumConstants| with MATLAB name |clib.ctMatlabInterface.kin_getEquilibriumConstants|
% C++ Signature: int kin_getEquilibriumConstants(int n,size_t len,double * kc)

kin_getEquilibriumConstantsDefinition = addFunction(libDef, ...
    "int kin_getEquilibriumConstants(int n,size_t len,double * kc)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_getEquilibriumConstants", ...
    "Description", "clib.ctMatlabInterface.kin_getEquilibriumConstants Representation of C++ function kin_getEquilibriumConstants."); % Modify help description values as needed.
defineArgument(kin_getEquilibriumConstantsDefinition, "n", "int32");
defineArgument(kin_getEquilibriumConstantsDefinition, "len", "uint64");
defineArgument(kin_getEquilibriumConstantsDefinition, "kc", "clib.array.ctMatlabInterface.Double", "input", "len"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(kin_getEquilibriumConstantsDefinition, "RetVal", "int32");
validate(kin_getEquilibriumConstantsDefinition);

%% C++ function |kin_getFwdRateConstants| with MATLAB name |clib.ctMatlabInterface.kin_getFwdRateConstants|
% C++ Signature: int kin_getFwdRateConstants(int n,size_t len,double * kfwd)

kin_getFwdRateConstantsDefinition = addFunction(libDef, ...
    "int kin_getFwdRateConstants(int n,size_t len,double * kfwd)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_getFwdRateConstants", ...
    "Description", "clib.ctMatlabInterface.kin_getFwdRateConstants Representation of C++ function kin_getFwdRateConstants."); % Modify help description values as needed.
defineArgument(kin_getFwdRateConstantsDefinition, "n", "int32");
defineArgument(kin_getFwdRateConstantsDefinition, "len", "uint64");
defineArgument(kin_getFwdRateConstantsDefinition, "kfwd", "clib.array.ctMatlabInterface.Double", "input", "len"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(kin_getFwdRateConstantsDefinition, "RetVal", "int32");
validate(kin_getFwdRateConstantsDefinition);

%% C++ function |kin_getRevRateConstants| with MATLAB name |clib.ctMatlabInterface.kin_getRevRateConstants|
% C++ Signature: int kin_getRevRateConstants(int n,int doIrreversible,size_t len,double * krev)

kin_getRevRateConstantsDefinition = addFunction(libDef, ...
    "int kin_getRevRateConstants(int n,int doIrreversible,size_t len,double * krev)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_getRevRateConstants", ...
    "Description", "clib.ctMatlabInterface.kin_getRevRateConstants Representation of C++ function kin_getRevRateConstants."); % Modify help description values as needed.
defineArgument(kin_getRevRateConstantsDefinition, "n", "int32");
defineArgument(kin_getRevRateConstantsDefinition, "doIrreversible", "int32");
defineArgument(kin_getRevRateConstantsDefinition, "len", "uint64");
defineArgument(kin_getRevRateConstantsDefinition, "krev", "clib.array.ctMatlabInterface.Double", "input", "len"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(kin_getRevRateConstantsDefinition, "RetVal", "int32");
validate(kin_getRevRateConstantsDefinition);

%% C++ function |kin_getDelta| with MATLAB name |clib.ctMatlabInterface.kin_getDelta|
% C++ Signature: int kin_getDelta(int n,int job,size_t len,double * delta)

kin_getDeltaDefinition = addFunction(libDef, ...
    "int kin_getDelta(int n,int job,size_t len,double * delta)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_getDelta", ...
    "Description", "clib.ctMatlabInterface.kin_getDelta Representation of C++ function kin_getDelta."); % Modify help description values as needed.
defineArgument(kin_getDeltaDefinition, "n", "int32");
defineArgument(kin_getDeltaDefinition, "job", "int32");
defineArgument(kin_getDeltaDefinition, "len", "uint64");
defineArgument(kin_getDeltaDefinition, "delta", "clib.array.ctMatlabInterface.Double", "input", "len"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(kin_getDeltaDefinition, "RetVal", "int32");
validate(kin_getDeltaDefinition);

%% C++ function |kin_getCreationRates| with MATLAB name |clib.ctMatlabInterface.kin_getCreationRates|
% C++ Signature: int kin_getCreationRates(int n,size_t len,double * cdot)

kin_getCreationRatesDefinition = addFunction(libDef, ...
    "int kin_getCreationRates(int n,size_t len,double * cdot)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_getCreationRates", ...
    "Description", "clib.ctMatlabInterface.kin_getCreationRates Representation of C++ function kin_getCreationRates."); % Modify help description values as needed.
defineArgument(kin_getCreationRatesDefinition, "n", "int32");
defineArgument(kin_getCreationRatesDefinition, "len", "uint64");
defineArgument(kin_getCreationRatesDefinition, "cdot", "clib.array.ctMatlabInterface.Double", "input", "len"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(kin_getCreationRatesDefinition, "RetVal", "int32");
validate(kin_getCreationRatesDefinition);

%% C++ function |kin_getDestructionRates| with MATLAB name |clib.ctMatlabInterface.kin_getDestructionRates|
% C++ Signature: int kin_getDestructionRates(int n,size_t len,double * ddot)

kin_getDestructionRatesDefinition = addFunction(libDef, ...
    "int kin_getDestructionRates(int n,size_t len,double * ddot)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_getDestructionRates", ...
    "Description", "clib.ctMatlabInterface.kin_getDestructionRates Representation of C++ function kin_getDestructionRates."); % Modify help description values as needed.
defineArgument(kin_getDestructionRatesDefinition, "n", "int32");
defineArgument(kin_getDestructionRatesDefinition, "len", "uint64");
defineArgument(kin_getDestructionRatesDefinition, "ddot", "clib.array.ctMatlabInterface.Double", "input", "len"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(kin_getDestructionRatesDefinition, "RetVal", "int32");
validate(kin_getDestructionRatesDefinition);

%% C++ function |kin_getNetProductionRates| with MATLAB name |clib.ctMatlabInterface.kin_getNetProductionRates|
% C++ Signature: int kin_getNetProductionRates(int n,size_t len,double * wdot)

kin_getNetProductionRatesDefinition = addFunction(libDef, ...
    "int kin_getNetProductionRates(int n,size_t len,double * wdot)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_getNetProductionRates", ...
    "Description", "clib.ctMatlabInterface.kin_getNetProductionRates Representation of C++ function kin_getNetProductionRates."); % Modify help description values as needed.
defineArgument(kin_getNetProductionRatesDefinition, "n", "int32");
defineArgument(kin_getNetProductionRatesDefinition, "len", "uint64");
defineArgument(kin_getNetProductionRatesDefinition, "wdot", "clib.array.ctMatlabInterface.Double", "input", "len"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(kin_getNetProductionRatesDefinition, "RetVal", "int32");
validate(kin_getNetProductionRatesDefinition);

%% C++ function |kin_getSourceTerms| with MATLAB name |clib.ctMatlabInterface.kin_getSourceTerms|
% C++ Signature: int kin_getSourceTerms(int n,size_t len,double * ydot)

kin_getSourceTermsDefinition = addFunction(libDef, ...
    "int kin_getSourceTerms(int n,size_t len,double * ydot)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_getSourceTerms", ...
    "Description", "clib.ctMatlabInterface.kin_getSourceTerms Representation of C++ function kin_getSourceTerms."); % Modify help description values as needed.
defineArgument(kin_getSourceTermsDefinition, "n", "int32");
defineArgument(kin_getSourceTermsDefinition, "len", "uint64");
defineArgument(kin_getSourceTermsDefinition, "ydot", "clib.array.ctMatlabInterface.Double", "input", "len"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(kin_getSourceTermsDefinition, "RetVal", "int32");
validate(kin_getSourceTermsDefinition);

%% C++ function |kin_multiplier| with MATLAB name |clib.ctMatlabInterface.kin_multiplier|
% C++ Signature: double kin_multiplier(int n,int i)

kin_multiplierDefinition = addFunction(libDef, ...
    "double kin_multiplier(int n,int i)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_multiplier", ...
    "Description", "clib.ctMatlabInterface.kin_multiplier Representation of C++ function kin_multiplier."); % Modify help description values as needed.
defineArgument(kin_multiplierDefinition, "n", "int32");
defineArgument(kin_multiplierDefinition, "i", "int32");
defineOutput(kin_multiplierDefinition, "RetVal", "double");
validate(kin_multiplierDefinition);

%% C++ function |kin_getReactionString| with MATLAB name |clib.ctMatlabInterface.kin_getReactionString|
% C++ Signature: int kin_getReactionString(int n,int i,int len,char * buf)

kin_getReactionStringDefinition = addFunction(libDef, ...
    "int kin_getReactionString(int n,int i,int len,char * buf)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_getReactionString", ...
    "Description", "clib.ctMatlabInterface.kin_getReactionString Representation of C++ function kin_getReactionString."); % Modify help description values as needed.
defineArgument(kin_getReactionStringDefinition, "n", "int32");
defineArgument(kin_getReactionStringDefinition, "i", "int32");
defineArgument(kin_getReactionStringDefinition, "len", "int32");
defineArgument(kin_getReactionStringDefinition, "buf", "clib.array.ctMatlabInterface.Char", "input", "len"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(kin_getReactionStringDefinition, "RetVal", "int32");
validate(kin_getReactionStringDefinition);

%% C++ function |kin_setMultiplier| with MATLAB name |clib.ctMatlabInterface.kin_setMultiplier|
% C++ Signature: int kin_setMultiplier(int n,int i,double v)

kin_setMultiplierDefinition = addFunction(libDef, ...
    "int kin_setMultiplier(int n,int i,double v)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_setMultiplier", ...
    "Description", "clib.ctMatlabInterface.kin_setMultiplier Representation of C++ function kin_setMultiplier."); % Modify help description values as needed.
defineArgument(kin_setMultiplierDefinition, "n", "int32");
defineArgument(kin_setMultiplierDefinition, "i", "int32");
defineArgument(kin_setMultiplierDefinition, "v", "double");
defineOutput(kin_setMultiplierDefinition, "RetVal", "int32");
validate(kin_setMultiplierDefinition);

%% C++ function |kin_isReversible| with MATLAB name |clib.ctMatlabInterface.kin_isReversible|
% C++ Signature: int kin_isReversible(int n,int i)

kin_isReversibleDefinition = addFunction(libDef, ...
    "int kin_isReversible(int n,int i)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_isReversible", ...
    "Description", "clib.ctMatlabInterface.kin_isReversible Representation of C++ function kin_isReversible."); % Modify help description values as needed.
defineArgument(kin_isReversibleDefinition, "n", "int32");
defineArgument(kin_isReversibleDefinition, "i", "int32");
defineOutput(kin_isReversibleDefinition, "RetVal", "int32");
validate(kin_isReversibleDefinition);

%% C++ function |kin_getType| with MATLAB name |clib.ctMatlabInterface.kin_getType|
% C++ Signature: int kin_getType(int n,size_t len,char * name)

kin_getTypeDefinition = addFunction(libDef, ...
    "int kin_getType(int n,size_t len,char * name)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_getType", ...
    "Description", "clib.ctMatlabInterface.kin_getType Representation of C++ function kin_getType."); % Modify help description values as needed.
defineArgument(kin_getTypeDefinition, "n", "int32");
defineArgument(kin_getTypeDefinition, "len", "uint64");
defineArgument(kin_getTypeDefinition, "name", "clib.array.ctMatlabInterface.Char", "input", "len"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(kin_getTypeDefinition, "RetVal", "int32");
validate(kin_getTypeDefinition);

%% C++ function |kin_start| with MATLAB name |clib.ctMatlabInterface.kin_start|
% C++ Signature: size_t kin_start(int n,int p)

kin_startDefinition = addFunction(libDef, ...
    "size_t kin_start(int n,int p)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_start", ...
    "Description", "clib.ctMatlabInterface.kin_start Representation of C++ function kin_start."); % Modify help description values as needed.
defineArgument(kin_startDefinition, "n", "int32");
defineArgument(kin_startDefinition, "p", "int32");
defineOutput(kin_startDefinition, "RetVal", "uint64");
validate(kin_startDefinition);

%% C++ function |kin_speciesIndex| with MATLAB name |clib.ctMatlabInterface.kin_speciesIndex|
% C++ Signature: size_t kin_speciesIndex(int n,char const * nm)

kin_speciesIndexDefinition = addFunction(libDef, ...
    "size_t kin_speciesIndex(int n,char const * nm)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_speciesIndex", ...
    "Description", "clib.ctMatlabInterface.kin_speciesIndex Representation of C++ function kin_speciesIndex."); % Modify help description values as needed.
defineArgument(kin_speciesIndexDefinition, "n", "int32");
defineArgument(kin_speciesIndexDefinition, "nm", "string", "input", "nullTerminated");
defineOutput(kin_speciesIndexDefinition, "RetVal", "uint64");
validate(kin_speciesIndexDefinition);

%% C++ function |kin_advanceCoverages| with MATLAB name |clib.ctMatlabInterface.kin_advanceCoverages|
% C++ Signature: int kin_advanceCoverages(int n,double tstep)

kin_advanceCoveragesDefinition = addFunction(libDef, ...
    "int kin_advanceCoverages(int n,double tstep)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_advanceCoverages", ...
    "Description", "clib.ctMatlabInterface.kin_advanceCoverages Representation of C++ function kin_advanceCoverages."); % Modify help description values as needed.
defineArgument(kin_advanceCoveragesDefinition, "n", "int32");
defineArgument(kin_advanceCoveragesDefinition, "tstep", "double");
defineOutput(kin_advanceCoveragesDefinition, "RetVal", "int32");
validate(kin_advanceCoveragesDefinition);

%% C++ function |kin_phase| with MATLAB name |clib.ctMatlabInterface.kin_phase|
% C++ Signature: size_t kin_phase(int n,size_t i)

kin_phaseDefinition = addFunction(libDef, ...
    "size_t kin_phase(int n,size_t i)", ...
    "MATLABName", "clib.ctMatlabInterface.kin_phase", ...
    "Description", "clib.ctMatlabInterface.kin_phase Representation of C++ function kin_phase."); % Modify help description values as needed.
defineArgument(kin_phaseDefinition, "n", "int32");
defineArgument(kin_phaseDefinition, "i", "uint64");
defineOutput(kin_phaseDefinition, "RetVal", "uint64");
validate(kin_phaseDefinition);

%% C++ function |trans_parent| with MATLAB name |clib.ctMatlabInterface.trans_parent|
% C++ Signature: int trans_parent(int n)

trans_parentDefinition = addFunction(libDef, ...
    "int trans_parent(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.trans_parent", ...
    "Description", "clib.ctMatlabInterface.trans_parent Representation of C++ function trans_parent."); % Modify help description values as needed.
defineArgument(trans_parentDefinition, "n", "int32");
defineOutput(trans_parentDefinition, "RetVal", "int32");
validate(trans_parentDefinition);

%% C++ function |trans_del| with MATLAB name |clib.ctMatlabInterface.trans_del|
% C++ Signature: int trans_del(int n)

trans_delDefinition = addFunction(libDef, ...
    "int trans_del(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.trans_del", ...
    "Description", "clib.ctMatlabInterface.trans_del Representation of C++ function trans_del." + newline + ...
    "no-op; object is managed by Solution"); % Modify help description values as needed.
defineArgument(trans_delDefinition, "n", "int32");
defineOutput(trans_delDefinition, "RetVal", "int32");
validate(trans_delDefinition);

%% C++ function |trans_transportModel| with MATLAB name |clib.ctMatlabInterface.trans_transportModel|
% C++ Signature: int trans_transportModel(int n,int lennm,char * nm)

trans_transportModelDefinition = addFunction(libDef, ...
    "int trans_transportModel(int n,int lennm,char * nm)", ...
    "MATLABName", "clib.ctMatlabInterface.trans_transportModel", ...
    "Description", "clib.ctMatlabInterface.trans_transportModel Representation of C++ function trans_transportModel."); % Modify help description values as needed.
defineArgument(trans_transportModelDefinition, "n", "int32");
defineArgument(trans_transportModelDefinition, "lennm", "int32");
defineArgument(trans_transportModelDefinition, "nm", "clib.array.ctMatlabInterface.Char", "input", "lennm"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(trans_transportModelDefinition, "RetVal", "int32");
validate(trans_transportModelDefinition);

%% C++ function |trans_viscosity| with MATLAB name |clib.ctMatlabInterface.trans_viscosity|
% C++ Signature: double trans_viscosity(int n)

trans_viscosityDefinition = addFunction(libDef, ...
    "double trans_viscosity(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.trans_viscosity", ...
    "Description", "clib.ctMatlabInterface.trans_viscosity Representation of C++ function trans_viscosity."); % Modify help description values as needed.
defineArgument(trans_viscosityDefinition, "n", "int32");
defineOutput(trans_viscosityDefinition, "RetVal", "double");
validate(trans_viscosityDefinition);

%% C++ function |trans_electricalConductivity| with MATLAB name |clib.ctMatlabInterface.trans_electricalConductivity|
% C++ Signature: double trans_electricalConductivity(int n)

trans_electricalConductivityDefinition = addFunction(libDef, ...
    "double trans_electricalConductivity(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.trans_electricalConductivity", ...
    "Description", "clib.ctMatlabInterface.trans_electricalConductivity Representation of C++ function trans_electricalConductivity."); % Modify help description values as needed.
defineArgument(trans_electricalConductivityDefinition, "n", "int32");
defineOutput(trans_electricalConductivityDefinition, "RetVal", "double");
validate(trans_electricalConductivityDefinition);

%% C++ function |trans_thermalConductivity| with MATLAB name |clib.ctMatlabInterface.trans_thermalConductivity|
% C++ Signature: double trans_thermalConductivity(int n)

trans_thermalConductivityDefinition = addFunction(libDef, ...
    "double trans_thermalConductivity(int n)", ...
    "MATLABName", "clib.ctMatlabInterface.trans_thermalConductivity", ...
    "Description", "clib.ctMatlabInterface.trans_thermalConductivity Representation of C++ function trans_thermalConductivity."); % Modify help description values as needed.
defineArgument(trans_thermalConductivityDefinition, "n", "int32");
defineOutput(trans_thermalConductivityDefinition, "RetVal", "double");
validate(trans_thermalConductivityDefinition);

%% C++ function |trans_getThermalDiffCoeffs| with MATLAB name |clib.ctMatlabInterface.trans_getThermalDiffCoeffs|
% C++ Signature: int trans_getThermalDiffCoeffs(int n,int ldt,double * dt)

trans_getThermalDiffCoeffsDefinition = addFunction(libDef, ...
    "int trans_getThermalDiffCoeffs(int n,int ldt,double * dt)", ...
    "MATLABName", "clib.ctMatlabInterface.trans_getThermalDiffCoeffs", ...
    "Description", "clib.ctMatlabInterface.trans_getThermalDiffCoeffs Representation of C++ function trans_getThermalDiffCoeffs."); % Modify help description values as needed.
defineArgument(trans_getThermalDiffCoeffsDefinition, "n", "int32");
defineArgument(trans_getThermalDiffCoeffsDefinition, "ldt", "int32");
defineArgument(trans_getThermalDiffCoeffsDefinition, "dt", "clib.array.ctMatlabInterface.Double", "input", "ldt"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(trans_getThermalDiffCoeffsDefinition, "RetVal", "int32");
validate(trans_getThermalDiffCoeffsDefinition);

%% C++ function |trans_getMixDiffCoeffs| with MATLAB name |clib.ctMatlabInterface.trans_getMixDiffCoeffs|
% C++ Signature: int trans_getMixDiffCoeffs(int n,int ld,double * d)

trans_getMixDiffCoeffsDefinition = addFunction(libDef, ...
    "int trans_getMixDiffCoeffs(int n,int ld,double * d)", ...
    "MATLABName", "clib.ctMatlabInterface.trans_getMixDiffCoeffs", ...
    "Description", "clib.ctMatlabInterface.trans_getMixDiffCoeffs Representation of C++ function trans_getMixDiffCoeffs."); % Modify help description values as needed.
defineArgument(trans_getMixDiffCoeffsDefinition, "n", "int32");
defineArgument(trans_getMixDiffCoeffsDefinition, "ld", "int32");
defineArgument(trans_getMixDiffCoeffsDefinition, "d", "clib.array.ctMatlabInterface.Double", "input", "ld"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(trans_getMixDiffCoeffsDefinition, "RetVal", "int32");
validate(trans_getMixDiffCoeffsDefinition);

%% C++ function |trans_getBinDiffCoeffs| with MATLAB name |clib.ctMatlabInterface.trans_getBinDiffCoeffs|
% C++ Signature: int trans_getBinDiffCoeffs(int n,int ld,double * d)

trans_getBinDiffCoeffsDefinition = addFunction(libDef, ...
    "int trans_getBinDiffCoeffs(int n,int ld,double * d)", ...
    "MATLABName", "clib.ctMatlabInterface.trans_getBinDiffCoeffs", ...
    "Description", "clib.ctMatlabInterface.trans_getBinDiffCoeffs Representation of C++ function trans_getBinDiffCoeffs."); % Modify help description values as needed.
defineArgument(trans_getBinDiffCoeffsDefinition, "n", "int32");
defineArgument(trans_getBinDiffCoeffsDefinition, "ld", "int32");
defineArgument(trans_getBinDiffCoeffsDefinition, "d", "clib.array.ctMatlabInterface.Double", "input", "ld"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(trans_getBinDiffCoeffsDefinition, "RetVal", "int32");
validate(trans_getBinDiffCoeffsDefinition);

%% C++ function |trans_getMultiDiffCoeffs| with MATLAB name |clib.ctMatlabInterface.trans_getMultiDiffCoeffs|
% C++ Signature: int trans_getMultiDiffCoeffs(int n,int ld,double * d)

trans_getMultiDiffCoeffsDefinition = addFunction(libDef, ...
    "int trans_getMultiDiffCoeffs(int n,int ld,double * d)", ...
    "MATLABName", "clib.ctMatlabInterface.trans_getMultiDiffCoeffs", ...
    "Description", "clib.ctMatlabInterface.trans_getMultiDiffCoeffs Representation of C++ function trans_getMultiDiffCoeffs."); % Modify help description values as needed.
defineArgument(trans_getMultiDiffCoeffsDefinition, "n", "int32");
defineArgument(trans_getMultiDiffCoeffsDefinition, "ld", "int32");
defineArgument(trans_getMultiDiffCoeffsDefinition, "d", "clib.array.ctMatlabInterface.Double", "input", "ld"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(trans_getMultiDiffCoeffsDefinition, "RetVal", "int32");
validate(trans_getMultiDiffCoeffsDefinition);

%% C++ function |trans_getMolarFluxes| with MATLAB name |clib.ctMatlabInterface.trans_getMolarFluxes|
% C++ Signature: int trans_getMolarFluxes(int n,double const * state1,double const * state2,double delta,double * fluxes)

% trans_getMolarFluxesDefinition = addFunction(libDef, ...
%     "int trans_getMolarFluxes(int n,double const * state1,double const * state2,double delta,double * fluxes)", ...
%     "MATLABName", "clib.ctMatlabInterface.trans_getMolarFluxes", ...
%     "Description", "clib.ctMatlabInterface.trans_getMolarFluxes Representation of C++ function trans_getMolarFluxes."); % Modify help description values as needed.
% defineArgument(trans_getMolarFluxesDefinition, "n", "int32");
% defineArgument(trans_getMolarFluxesDefinition, "state1", "clib.array.ctMatlabInterface.Double", "input", 2); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
% defineArgument(trans_getMolarFluxesDefinition, "state2", "clib.array.ctMatlabInterface.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
% defineArgument(trans_getMolarFluxesDefinition, "delta", "double");
% defineArgument(trans_getMolarFluxesDefinition, "fluxes", "clib.array.ctMatlabInterface.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
% defineOutput(trans_getMolarFluxesDefinition, "RetVal", "int32");
% validate(trans_getMolarFluxesDefinition);

%% C++ function |trans_getMassFluxes| with MATLAB name |clib.ctMatlabInterface.trans_getMassFluxes|
% C++ Signature: int trans_getMassFluxes(int n,double const * state1,double const * state2,double delta,double * fluxes)

% trans_getMassFluxesDefinition = addFunction(libDef, ...
%     "int trans_getMassFluxes(int n,double const * state1,double const * state2,double delta,double * fluxes)", ...
%     "MATLABName", "clib.ctMatlabInterface.trans_getMassFluxes", ...
%     "Description", "clib.ctMatlabInterface.trans_getMassFluxes Representation of C++ function trans_getMassFluxes."); % Modify help description values as needed.
% defineArgument(trans_getMassFluxesDefinition, "n", "int32");
% defineArgument(trans_getMassFluxesDefinition, "state1", "clib.array.ctMatlabInterface.Double", "input", 2); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
% defineArgument(trans_getMassFluxesDefinition, "state2", "clib.array.ctMatlabInterface.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
% defineArgument(trans_getMassFluxesDefinition, "delta", "double");
% defineArgument(trans_getMassFluxesDefinition, "fluxes", "clib.array.ctMatlabInterface.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
% defineOutput(trans_getMassFluxesDefinition, "RetVal", "int32");
% validate(trans_getMassFluxesDefinition);

%% C++ function |ct_getCanteraError| with MATLAB name |clib.ctMatlabInterface.ct_getCanteraError|
% C++ Signature: int ct_getCanteraError(int buflen,char * buf)

ct_getCanteraErrorDefinition = addFunction(libDef, ...
    "int ct_getCanteraError(int buflen,char * buf)", ...
    "MATLABName", "clib.ctMatlabInterface.ct_getCanteraError", ...
    "Description", "clib.ctMatlabInterface.ct_getCanteraError Representation of C++ function ct_getCanteraError."); % Modify help description values as needed.
defineArgument(ct_getCanteraErrorDefinition, "buflen", "int32");
defineArgument(ct_getCanteraErrorDefinition, "buf", "clib.array.ctMatlabInterface.Char", "input", "buflen"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(ct_getCanteraErrorDefinition, "RetVal", "int32");
validate(ct_getCanteraErrorDefinition);

%% C++ function |ct_setLogWriter| with MATLAB name |clib.ctMatlabInterface.ct_setLogWriter|
% C++ Signature: int ct_setLogWriter(void * logger)

%ct_setLogWriterDefinition = addFunction(libDef, ...
%    "int ct_setLogWriter(void * logger)", ...
%    "MATLABName", "clib.ctMatlabInterface.ct_setLogWriter", ...
%    "Description", "clib.ctMatlabInterface.ct_setLogWriter Representation of C++ function ct_setLogWriter."); % Modify help description values as needed.
%defineArgument(ct_setLogWriterDefinition, "logger", <MLTYPE>, <DIRECTION>, 1); % <MLTYPE> can be primitive type, user-defined type, clib.array type, or a list of existing typedef names for void*.
%defineOutput(ct_setLogWriterDefinition, "RetVal", "int32");
%validate(ct_setLogWriterDefinition);

%% C++ function |ct_setLogCallback| with MATLAB name |clib.ctMatlabInterface.ct_setLogCallback|
% C++ Signature: int ct_setLogCallback(LogCallback writer)

ct_setLogCallbackDefinition = addFunction(libDef, ...
    "int ct_setLogCallback(LogCallback writer)", ...
    "MATLABName", "clib.ctMatlabInterface.ct_setLogCallback", ...
    "Description", "clib.ctMatlabInterface.ct_setLogCallback Representation of C++ function ct_setLogCallback."); % Modify help description values as needed.
defineArgument(ct_setLogCallbackDefinition, "writer", "clib.ctMatlabInterface.LogCallback");
defineOutput(ct_setLogCallbackDefinition, "RetVal", "int32");
validate(ct_setLogCallbackDefinition);

%% C++ function |ct_addCanteraDirectory| with MATLAB name |clib.ctMatlabInterface.ct_addCanteraDirectory|
% C++ Signature: int ct_addCanteraDirectory(size_t buflen,char const * buf)

ct_addCanteraDirectoryDefinition = addFunction(libDef, ...
    "int ct_addCanteraDirectory(size_t buflen,char const * buf)", ...
    "MATLABName", "clib.ctMatlabInterface.ct_addCanteraDirectory", ...
    "Description", "clib.ctMatlabInterface.ct_addCanteraDirectory Representation of C++ function ct_addCanteraDirectory."); % Modify help description values as needed.
defineArgument(ct_addCanteraDirectoryDefinition, "buflen", "uint64");
defineArgument(ct_addCanteraDirectoryDefinition, "buf", "string", "input", "nullTerminated");
defineOutput(ct_addCanteraDirectoryDefinition, "RetVal", "int32");
validate(ct_addCanteraDirectoryDefinition);

%% C++ function |ct_getDataDirectories| with MATLAB name |clib.ctMatlabInterface.ct_getDataDirectories|
% C++ Signature: int ct_getDataDirectories(char const * sep,int buflen,char * buf)

ct_getDataDirectoriesDefinition = addFunction(libDef, ...
    "int ct_getDataDirectories(char const * sep,int buflen,char * buf)", ...
    "MATLABName", "clib.ctMatlabInterface.ct_getDataDirectories", ...
    "Description", "clib.ctMatlabInterface.ct_getDataDirectories Representation of C++ function ct_getDataDirectories."); % Modify help description values as needed.
defineArgument(ct_getDataDirectoriesDefinition, "sep", "string", "input", "nullTerminated");
defineArgument(ct_getDataDirectoriesDefinition, "buflen", "int32");
defineArgument(ct_getDataDirectoriesDefinition, "buf", "clib.array.ctMatlabInterface.Char", "input", "buflen"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(ct_getDataDirectoriesDefinition, "RetVal", "int32");
validate(ct_getDataDirectoriesDefinition);

%% C++ function |ct_getCanteraVersion| with MATLAB name |clib.ctMatlabInterface.ct_getCanteraVersion|
% C++ Signature: int ct_getCanteraVersion(int buflen,char * buf)

ct_getCanteraVersionDefinition = addFunction(libDef, ...
    "int ct_getCanteraVersion(int buflen,char * buf)", ...
    "MATLABName", "clib.ctMatlabInterface.ct_getCanteraVersion", ...
    "Description", "clib.ctMatlabInterface.ct_getCanteraVersion Representation of C++ function ct_getCanteraVersion."); % Modify help description values as needed.
defineArgument(ct_getCanteraVersionDefinition, "buflen", "int32");
defineArgument(ct_getCanteraVersionDefinition, "buf", "clib.array.ctMatlabInterface.Char", "input", "buflen"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(ct_getCanteraVersionDefinition, "RetVal", "int32");
validate(ct_getCanteraVersionDefinition);

%% C++ function |ct_getGitCommit| with MATLAB name |clib.ctMatlabInterface.ct_getGitCommit|
% C++ Signature: int ct_getGitCommit(int buflen,char * buf)

ct_getGitCommitDefinition = addFunction(libDef, ...
    "int ct_getGitCommit(int buflen,char * buf)", ...
    "MATLABName", "clib.ctMatlabInterface.ct_getGitCommit", ...
    "Description", "clib.ctMatlabInterface.ct_getGitCommit Representation of C++ function ct_getGitCommit."); % Modify help description values as needed.
defineArgument(ct_getGitCommitDefinition, "buflen", "int32");
defineArgument(ct_getGitCommitDefinition, "buf", "clib.array.ctMatlabInterface.Char", "input", "buflen"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(ct_getGitCommitDefinition, "RetVal", "int32");
validate(ct_getGitCommitDefinition);

%% C++ function |ct_suppress_thermo_warnings| with MATLAB name |clib.ctMatlabInterface.ct_suppress_thermo_warnings|
% C++ Signature: int ct_suppress_thermo_warnings(int suppress)

ct_suppress_thermo_warningsDefinition = addFunction(libDef, ...
    "int ct_suppress_thermo_warnings(int suppress)", ...
    "MATLABName", "clib.ctMatlabInterface.ct_suppress_thermo_warnings", ...
    "Description", "clib.ctMatlabInterface.ct_suppress_thermo_warnings Representation of C++ function ct_suppress_thermo_warnings."); % Modify help description values as needed.
defineArgument(ct_suppress_thermo_warningsDefinition, "suppress", "int32");
defineOutput(ct_suppress_thermo_warningsDefinition, "RetVal", "int32");
validate(ct_suppress_thermo_warningsDefinition);

%% C++ function |ct_use_legacy_rate_constants| with MATLAB name |clib.ctMatlabInterface.ct_use_legacy_rate_constants|
% C++ Signature: int ct_use_legacy_rate_constants(int legacy)

ct_use_legacy_rate_constantsDefinition = addFunction(libDef, ...
    "int ct_use_legacy_rate_constants(int legacy)", ...
    "MATLABName", "clib.ctMatlabInterface.ct_use_legacy_rate_constants", ...
    "Description", "clib.ctMatlabInterface.ct_use_legacy_rate_constants Representation of C++ function ct_use_legacy_rate_constants."); % Modify help description values as needed.
defineArgument(ct_use_legacy_rate_constantsDefinition, "legacy", "int32");
defineOutput(ct_use_legacy_rate_constantsDefinition, "RetVal", "int32");
validate(ct_use_legacy_rate_constantsDefinition);

%% C++ function |ct_clearStorage| with MATLAB name |clib.ctMatlabInterface.ct_clearStorage|
% C++ Signature: int ct_clearStorage()

ct_clearStorageDefinition = addFunction(libDef, ...
    "int ct_clearStorage()", ...
    "MATLABName", "clib.ctMatlabInterface.ct_clearStorage", ...
    "Description", "clib.ctMatlabInterface.ct_clearStorage Representation of C++ function ct_clearStorage."); % Modify help description values as needed.
defineOutput(ct_clearStorageDefinition, "RetVal", "int32");
validate(ct_clearStorageDefinition);

%% C++ function |ct_resetStorage| with MATLAB name |clib.ctMatlabInterface.ct_resetStorage|
% C++ Signature: int ct_resetStorage()

ct_resetStorageDefinition = addFunction(libDef, ...
    "int ct_resetStorage()", ...
    "MATLABName", "clib.ctMatlabInterface.ct_resetStorage", ...
    "Description", "clib.ctMatlabInterface.ct_resetStorage Representation of C++ function ct_resetStorage."); % Modify help description values as needed.
defineOutput(ct_resetStorageDefinition, "RetVal", "int32");
validate(ct_resetStorageDefinition);

%% C++ function |func_check| with MATLAB name |clib.ctMatlabInterface.func_check|
% C++ Signature: int func_check(char const * type,size_t len,char * buf)

func_checkDefinition = addFunction(libDef, ...
    "int func_check(char const * type,size_t len,char * buf)", ...
    "MATLABName", "clib.ctMatlabInterface.func_check", ...
    "Description", "clib.ctMatlabInterface.func_check Representation of C++ function func_check." + newline + ...
    "@since New in %Cantera 3.1"); % Modify help description values as needed.
defineArgument(func_checkDefinition, "type", "string", "input", "nullTerminated");
defineArgument(func_checkDefinition, "len", "uint64");
defineArgument(func_checkDefinition, "buf", "clib.array.ctMatlabInterface.Char", "input", "len"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(func_checkDefinition, "RetVal", "int32");
validate(func_checkDefinition);

%% C++ function |func_new_basic| with MATLAB name |clib.ctMatlabInterface.func_new_basic|
% C++ Signature: int func_new_basic(char const * type,double c)

func_new_basicDefinition = addFunction(libDef, ...
    "int func_new_basic(char const * type,double c)", ...
    "MATLABName", "clib.ctMatlabInterface.func_new_basic", ...
    "Description", "clib.ctMatlabInterface.func_new_basic Representation of C++ function func_new_basic."); % Modify help description values as needed.
defineArgument(func_new_basicDefinition, "type", "string", "input", "nullTerminated");
defineArgument(func_new_basicDefinition, "c", "double");
defineOutput(func_new_basicDefinition, "RetVal", "int32");
validate(func_new_basicDefinition);

%% C++ function |func_new_advanced| with MATLAB name |clib.ctMatlabInterface.func_new_advanced|
% C++ Signature: int func_new_advanced(char const * type,size_t lenp,double const * p)

func_new_advancedDefinition = addFunction(libDef, ...
    "int func_new_advanced(char const * type,size_t lenp,double const * p)", ...
    "MATLABName", "clib.ctMatlabInterface.func_new_advanced", ...
    "Description", "clib.ctMatlabInterface.func_new_advanced Representation of C++ function func_new_advanced."); % Modify help description values as needed.
defineArgument(func_new_advancedDefinition, "type", "string", "input", "nullTerminated");
defineArgument(func_new_advancedDefinition, "lenp", "uint64");
defineArgument(func_new_advancedDefinition, "p", "clib.array.ctMatlabInterface.Double", "input", "lenp"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(func_new_advancedDefinition, "RetVal", "int32");
validate(func_new_advancedDefinition);

%% C++ function |func_new_compound| with MATLAB name |clib.ctMatlabInterface.func_new_compound|
% C++ Signature: int func_new_compound(char const * type,int a,int b)

func_new_compoundDefinition = addFunction(libDef, ...
    "int func_new_compound(char const * type,int a,int b)", ...
    "MATLABName", "clib.ctMatlabInterface.func_new_compound", ...
    "Description", "clib.ctMatlabInterface.func_new_compound Representation of C++ function func_new_compound."); % Modify help description values as needed.
defineArgument(func_new_compoundDefinition, "type", "string", "input", "nullTerminated");
defineArgument(func_new_compoundDefinition, "a", "int32");
defineArgument(func_new_compoundDefinition, "b", "int32");
defineOutput(func_new_compoundDefinition, "RetVal", "int32");
validate(func_new_compoundDefinition);

%% C++ function |func_new_modified| with MATLAB name |clib.ctMatlabInterface.func_new_modified|
% C++ Signature: int func_new_modified(char const * type,int a,double c)

func_new_modifiedDefinition = addFunction(libDef, ...
    "int func_new_modified(char const * type,int a,double c)", ...
    "MATLABName", "clib.ctMatlabInterface.func_new_modified", ...
    "Description", "clib.ctMatlabInterface.func_new_modified Representation of C++ function func_new_modified."); % Modify help description values as needed.
defineArgument(func_new_modifiedDefinition, "type", "string", "input", "nullTerminated");
defineArgument(func_new_modifiedDefinition, "a", "int32");
defineArgument(func_new_modifiedDefinition, "c", "double");
defineOutput(func_new_modifiedDefinition, "RetVal", "int32");
validate(func_new_modifiedDefinition);

%% C++ function |func_new_sum| with MATLAB name |clib.ctMatlabInterface.func_new_sum|
% C++ Signature: int func_new_sum(int a,int b)

func_new_sumDefinition = addFunction(libDef, ...
    "int func_new_sum(int a,int b)", ...
    "MATLABName", "clib.ctMatlabInterface.func_new_sum", ...
    "Description", "clib.ctMatlabInterface.func_new_sum Representation of C++ function func_new_sum."); % Modify help description values as needed.
defineArgument(func_new_sumDefinition, "a", "int32");
defineArgument(func_new_sumDefinition, "b", "int32");
defineOutput(func_new_sumDefinition, "RetVal", "int32");
validate(func_new_sumDefinition);

%% C++ function |func_new_diff| with MATLAB name |clib.ctMatlabInterface.func_new_diff|
% C++ Signature: int func_new_diff(int a,int b)

func_new_diffDefinition = addFunction(libDef, ...
    "int func_new_diff(int a,int b)", ...
    "MATLABName", "clib.ctMatlabInterface.func_new_diff", ...
    "Description", "clib.ctMatlabInterface.func_new_diff Representation of C++ function func_new_diff."); % Modify help description values as needed.
defineArgument(func_new_diffDefinition, "a", "int32");
defineArgument(func_new_diffDefinition, "b", "int32");
defineOutput(func_new_diffDefinition, "RetVal", "int32");
validate(func_new_diffDefinition);

%% C++ function |func_new_prod| with MATLAB name |clib.ctMatlabInterface.func_new_prod|
% C++ Signature: int func_new_prod(int a,int b)

func_new_prodDefinition = addFunction(libDef, ...
    "int func_new_prod(int a,int b)", ...
    "MATLABName", "clib.ctMatlabInterface.func_new_prod", ...
    "Description", "clib.ctMatlabInterface.func_new_prod Representation of C++ function func_new_prod."); % Modify help description values as needed.
defineArgument(func_new_prodDefinition, "a", "int32");
defineArgument(func_new_prodDefinition, "b", "int32");
defineOutput(func_new_prodDefinition, "RetVal", "int32");
validate(func_new_prodDefinition);

%% C++ function |func_new_ratio| with MATLAB name |clib.ctMatlabInterface.func_new_ratio|
% C++ Signature: int func_new_ratio(int a,int b)

func_new_ratioDefinition = addFunction(libDef, ...
    "int func_new_ratio(int a,int b)", ...
    "MATLABName", "clib.ctMatlabInterface.func_new_ratio", ...
    "Description", "clib.ctMatlabInterface.func_new_ratio Representation of C++ function func_new_ratio."); % Modify help description values as needed.
defineArgument(func_new_ratioDefinition, "a", "int32");
defineArgument(func_new_ratioDefinition, "b", "int32");
defineOutput(func_new_ratioDefinition, "RetVal", "int32");
validate(func_new_ratioDefinition);

%% C++ function |func_del| with MATLAB name |clib.ctMatlabInterface.func_del|
% C++ Signature: int func_del(int i)

func_delDefinition = addFunction(libDef, ...
    "int func_del(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.func_del", ...
    "Description", "clib.ctMatlabInterface.func_del Representation of C++ function func_del."); % Modify help description values as needed.
defineArgument(func_delDefinition, "i", "int32");
defineOutput(func_delDefinition, "RetVal", "int32");
validate(func_delDefinition);

%% C++ function |func_type| with MATLAB name |clib.ctMatlabInterface.func_type|
% C++ Signature: int func_type(int i,size_t lennm,char * nm)

func_typeDefinition = addFunction(libDef, ...
    "int func_type(int i,size_t lennm,char * nm)", ...
    "MATLABName", "clib.ctMatlabInterface.func_type", ...
    "Description", "clib.ctMatlabInterface.func_type Representation of C++ function func_type."); % Modify help description values as needed.
defineArgument(func_typeDefinition, "i", "int32");
defineArgument(func_typeDefinition, "lennm", "uint64");
defineArgument(func_typeDefinition, "nm", "clib.array.ctMatlabInterface.Char", "input", "lennm"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(func_typeDefinition, "RetVal", "int32");
validate(func_typeDefinition);

%% C++ function |func_value| with MATLAB name |clib.ctMatlabInterface.func_value|
% C++ Signature: double func_value(int i,double t)

func_valueDefinition = addFunction(libDef, ...
    "double func_value(int i,double t)", ...
    "MATLABName", "clib.ctMatlabInterface.func_value", ...
    "Description", "clib.ctMatlabInterface.func_value Representation of C++ function func_value."); % Modify help description values as needed.
defineArgument(func_valueDefinition, "i", "int32");
defineArgument(func_valueDefinition, "t", "double");
defineOutput(func_valueDefinition, "RetVal", "double");
validate(func_valueDefinition);

%% C++ function |func_derivative| with MATLAB name |clib.ctMatlabInterface.func_derivative|
% C++ Signature: int func_derivative(int i)

func_derivativeDefinition = addFunction(libDef, ...
    "int func_derivative(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.func_derivative", ...
    "Description", "clib.ctMatlabInterface.func_derivative Representation of C++ function func_derivative."); % Modify help description values as needed.
defineArgument(func_derivativeDefinition, "i", "int32");
defineOutput(func_derivativeDefinition, "RetVal", "int32");
validate(func_derivativeDefinition);

%% C++ function |func_duplicate| with MATLAB name |clib.ctMatlabInterface.func_duplicate|
% C++ Signature: int func_duplicate(int i)

func_duplicateDefinition = addFunction(libDef, ...
    "int func_duplicate(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.func_duplicate", ...
    "Description", "clib.ctMatlabInterface.func_duplicate Representation of C++ function func_duplicate."); % Modify help description values as needed.
defineArgument(func_duplicateDefinition, "i", "int32");
defineOutput(func_duplicateDefinition, "RetVal", "int32");
validate(func_duplicateDefinition);

%% C++ function |func_write| with MATLAB name |clib.ctMatlabInterface.func_write|
% C++ Signature: int func_write(int i,char const * arg,size_t lennm,char * nm)

func_writeDefinition = addFunction(libDef, ...
    "int func_write(int i,char const * arg,size_t lennm,char * nm)", ...
    "MATLABName", "clib.ctMatlabInterface.func_write", ...
    "Description", "clib.ctMatlabInterface.func_write Representation of C++ function func_write."); % Modify help description values as needed.
defineArgument(func_writeDefinition, "i", "int32");
defineArgument(func_writeDefinition, "arg", "string", "input", "nullTerminated");
defineArgument(func_writeDefinition, "lennm", "uint64");
defineArgument(func_writeDefinition, "nm", "clib.array.ctMatlabInterface.Char", "input", "lennm"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(func_writeDefinition, "RetVal", "int32");
validate(func_writeDefinition);

%% C++ function |ct_clearFunc| with MATLAB name |clib.ctMatlabInterface.ct_clearFunc|
% C++ Signature: int ct_clearFunc()

ct_clearFuncDefinition = addFunction(libDef, ...
    "int ct_clearFunc()", ...
    "MATLABName", "clib.ctMatlabInterface.ct_clearFunc", ...
    "Description", "clib.ctMatlabInterface.ct_clearFunc Representation of C++ function ct_clearFunc."); % Modify help description values as needed.
defineOutput(ct_clearFuncDefinition, "RetVal", "int32");
validate(ct_clearFuncDefinition);

%% C++ function |mix_new| with MATLAB name |clib.ctMatlabInterface.mix_new|
% C++ Signature: int mix_new()

mix_newDefinition = addFunction(libDef, ...
    "int mix_new()", ...
    "MATLABName", "clib.ctMatlabInterface.mix_new", ...
    "Description", "clib.ctMatlabInterface.mix_new Representation of C++ function mix_new."); % Modify help description values as needed.
defineOutput(mix_newDefinition, "RetVal", "int32");
validate(mix_newDefinition);

%% C++ function |mix_del| with MATLAB name |clib.ctMatlabInterface.mix_del|
% C++ Signature: int mix_del(int i)

mix_delDefinition = addFunction(libDef, ...
    "int mix_del(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_del", ...
    "Description", "clib.ctMatlabInterface.mix_del Representation of C++ function mix_del."); % Modify help description values as needed.
defineArgument(mix_delDefinition, "i", "int32");
defineOutput(mix_delDefinition, "RetVal", "int32");
validate(mix_delDefinition);

%% C++ function |ct_clearMix| with MATLAB name |clib.ctMatlabInterface.ct_clearMix|
% C++ Signature: int ct_clearMix()

ct_clearMixDefinition = addFunction(libDef, ...
    "int ct_clearMix()", ...
    "MATLABName", "clib.ctMatlabInterface.ct_clearMix", ...
    "Description", "clib.ctMatlabInterface.ct_clearMix Representation of C++ function ct_clearMix."); % Modify help description values as needed.
defineOutput(ct_clearMixDefinition, "RetVal", "int32");
validate(ct_clearMixDefinition);

%% C++ function |mix_addPhase| with MATLAB name |clib.ctMatlabInterface.mix_addPhase|
% C++ Signature: int mix_addPhase(int i,int j,double moles)

mix_addPhaseDefinition = addFunction(libDef, ...
    "int mix_addPhase(int i,int j,double moles)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_addPhase", ...
    "Description", "clib.ctMatlabInterface.mix_addPhase Representation of C++ function mix_addPhase."); % Modify help description values as needed.
defineArgument(mix_addPhaseDefinition, "i", "int32");
defineArgument(mix_addPhaseDefinition, "j", "int32");
defineArgument(mix_addPhaseDefinition, "moles", "double");
defineOutput(mix_addPhaseDefinition, "RetVal", "int32");
validate(mix_addPhaseDefinition);

%% C++ function |mix_init| with MATLAB name |clib.ctMatlabInterface.mix_init|
% C++ Signature: int mix_init(int i)

mix_initDefinition = addFunction(libDef, ...
    "int mix_init(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_init", ...
    "Description", "clib.ctMatlabInterface.mix_init Representation of C++ function mix_init."); % Modify help description values as needed.
defineArgument(mix_initDefinition, "i", "int32");
defineOutput(mix_initDefinition, "RetVal", "int32");
validate(mix_initDefinition);

%% C++ function |mix_updatePhases| with MATLAB name |clib.ctMatlabInterface.mix_updatePhases|
% C++ Signature: int mix_updatePhases(int i)

mix_updatePhasesDefinition = addFunction(libDef, ...
    "int mix_updatePhases(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_updatePhases", ...
    "Description", "clib.ctMatlabInterface.mix_updatePhases Representation of C++ function mix_updatePhases."); % Modify help description values as needed.
defineArgument(mix_updatePhasesDefinition, "i", "int32");
defineOutput(mix_updatePhasesDefinition, "RetVal", "int32");
validate(mix_updatePhasesDefinition);

%% C++ function |mix_nElements| with MATLAB name |clib.ctMatlabInterface.mix_nElements|
% C++ Signature: size_t mix_nElements(int i)

mix_nElementsDefinition = addFunction(libDef, ...
    "size_t mix_nElements(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_nElements", ...
    "Description", "clib.ctMatlabInterface.mix_nElements Representation of C++ function mix_nElements."); % Modify help description values as needed.
defineArgument(mix_nElementsDefinition, "i", "int32");
defineOutput(mix_nElementsDefinition, "RetVal", "uint64");
validate(mix_nElementsDefinition);

%% C++ function |mix_elementIndex| with MATLAB name |clib.ctMatlabInterface.mix_elementIndex|
% C++ Signature: size_t mix_elementIndex(int i,char const * name)

mix_elementIndexDefinition = addFunction(libDef, ...
    "size_t mix_elementIndex(int i,char const * name)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_elementIndex", ...
    "Description", "clib.ctMatlabInterface.mix_elementIndex Representation of C++ function mix_elementIndex."); % Modify help description values as needed.
defineArgument(mix_elementIndexDefinition, "i", "int32");
defineArgument(mix_elementIndexDefinition, "name", "string", "input", "nullTerminated");
defineOutput(mix_elementIndexDefinition, "RetVal", "uint64");
validate(mix_elementIndexDefinition);

%% C++ function |mix_speciesIndex| with MATLAB name |clib.ctMatlabInterface.mix_speciesIndex|
% C++ Signature: size_t mix_speciesIndex(int i,int k,int p)

mix_speciesIndexDefinition = addFunction(libDef, ...
    "size_t mix_speciesIndex(int i,int k,int p)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_speciesIndex", ...
    "Description", "clib.ctMatlabInterface.mix_speciesIndex Representation of C++ function mix_speciesIndex."); % Modify help description values as needed.
defineArgument(mix_speciesIndexDefinition, "i", "int32");
defineArgument(mix_speciesIndexDefinition, "k", "int32");
defineArgument(mix_speciesIndexDefinition, "p", "int32");
defineOutput(mix_speciesIndexDefinition, "RetVal", "uint64");
validate(mix_speciesIndexDefinition);

%% C++ function |mix_nSpecies| with MATLAB name |clib.ctMatlabInterface.mix_nSpecies|
% C++ Signature: size_t mix_nSpecies(int i)

mix_nSpeciesDefinition = addFunction(libDef, ...
    "size_t mix_nSpecies(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_nSpecies", ...
    "Description", "clib.ctMatlabInterface.mix_nSpecies Representation of C++ function mix_nSpecies."); % Modify help description values as needed.
defineArgument(mix_nSpeciesDefinition, "i", "int32");
defineOutput(mix_nSpeciesDefinition, "RetVal", "uint64");
validate(mix_nSpeciesDefinition);

%% C++ function |mix_setTemperature| with MATLAB name |clib.ctMatlabInterface.mix_setTemperature|
% C++ Signature: int mix_setTemperature(int i,double t)

mix_setTemperatureDefinition = addFunction(libDef, ...
    "int mix_setTemperature(int i,double t)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_setTemperature", ...
    "Description", "clib.ctMatlabInterface.mix_setTemperature Representation of C++ function mix_setTemperature."); % Modify help description values as needed.
defineArgument(mix_setTemperatureDefinition, "i", "int32");
defineArgument(mix_setTemperatureDefinition, "t", "double");
defineOutput(mix_setTemperatureDefinition, "RetVal", "int32");
validate(mix_setTemperatureDefinition);

%% C++ function |mix_temperature| with MATLAB name |clib.ctMatlabInterface.mix_temperature|
% C++ Signature: double mix_temperature(int i)

mix_temperatureDefinition = addFunction(libDef, ...
    "double mix_temperature(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_temperature", ...
    "Description", "clib.ctMatlabInterface.mix_temperature Representation of C++ function mix_temperature."); % Modify help description values as needed.
defineArgument(mix_temperatureDefinition, "i", "int32");
defineOutput(mix_temperatureDefinition, "RetVal", "double");
validate(mix_temperatureDefinition);

%% C++ function |mix_minTemp| with MATLAB name |clib.ctMatlabInterface.mix_minTemp|
% C++ Signature: double mix_minTemp(int i)

mix_minTempDefinition = addFunction(libDef, ...
    "double mix_minTemp(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_minTemp", ...
    "Description", "clib.ctMatlabInterface.mix_minTemp Representation of C++ function mix_minTemp."); % Modify help description values as needed.
defineArgument(mix_minTempDefinition, "i", "int32");
defineOutput(mix_minTempDefinition, "RetVal", "double");
validate(mix_minTempDefinition);

%% C++ function |mix_maxTemp| with MATLAB name |clib.ctMatlabInterface.mix_maxTemp|
% C++ Signature: double mix_maxTemp(int i)

mix_maxTempDefinition = addFunction(libDef, ...
    "double mix_maxTemp(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_maxTemp", ...
    "Description", "clib.ctMatlabInterface.mix_maxTemp Representation of C++ function mix_maxTemp."); % Modify help description values as needed.
defineArgument(mix_maxTempDefinition, "i", "int32");
defineOutput(mix_maxTempDefinition, "RetVal", "double");
validate(mix_maxTempDefinition);

%% C++ function |mix_charge| with MATLAB name |clib.ctMatlabInterface.mix_charge|
% C++ Signature: double mix_charge(int i)

mix_chargeDefinition = addFunction(libDef, ...
    "double mix_charge(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_charge", ...
    "Description", "clib.ctMatlabInterface.mix_charge Representation of C++ function mix_charge."); % Modify help description values as needed.
defineArgument(mix_chargeDefinition, "i", "int32");
defineOutput(mix_chargeDefinition, "RetVal", "double");
validate(mix_chargeDefinition);

%% C++ function |mix_phaseCharge| with MATLAB name |clib.ctMatlabInterface.mix_phaseCharge|
% C++ Signature: double mix_phaseCharge(int i,int p)

mix_phaseChargeDefinition = addFunction(libDef, ...
    "double mix_phaseCharge(int i,int p)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_phaseCharge", ...
    "Description", "clib.ctMatlabInterface.mix_phaseCharge Representation of C++ function mix_phaseCharge."); % Modify help description values as needed.
defineArgument(mix_phaseChargeDefinition, "i", "int32");
defineArgument(mix_phaseChargeDefinition, "p", "int32");
defineOutput(mix_phaseChargeDefinition, "RetVal", "double");
validate(mix_phaseChargeDefinition);

%% C++ function |mix_setPressure| with MATLAB name |clib.ctMatlabInterface.mix_setPressure|
% C++ Signature: int mix_setPressure(int i,double p)

mix_setPressureDefinition = addFunction(libDef, ...
    "int mix_setPressure(int i,double p)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_setPressure", ...
    "Description", "clib.ctMatlabInterface.mix_setPressure Representation of C++ function mix_setPressure."); % Modify help description values as needed.
defineArgument(mix_setPressureDefinition, "i", "int32");
defineArgument(mix_setPressureDefinition, "p", "double");
defineOutput(mix_setPressureDefinition, "RetVal", "int32");
validate(mix_setPressureDefinition);

%% C++ function |mix_pressure| with MATLAB name |clib.ctMatlabInterface.mix_pressure|
% C++ Signature: double mix_pressure(int i)

mix_pressureDefinition = addFunction(libDef, ...
    "double mix_pressure(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_pressure", ...
    "Description", "clib.ctMatlabInterface.mix_pressure Representation of C++ function mix_pressure."); % Modify help description values as needed.
defineArgument(mix_pressureDefinition, "i", "int32");
defineOutput(mix_pressureDefinition, "RetVal", "double");
validate(mix_pressureDefinition);

%% C++ function |mix_nAtoms| with MATLAB name |clib.ctMatlabInterface.mix_nAtoms|
% C++ Signature: double mix_nAtoms(int i,int k,int m)

mix_nAtomsDefinition = addFunction(libDef, ...
    "double mix_nAtoms(int i,int k,int m)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_nAtoms", ...
    "Description", "clib.ctMatlabInterface.mix_nAtoms Representation of C++ function mix_nAtoms."); % Modify help description values as needed.
defineArgument(mix_nAtomsDefinition, "i", "int32");
defineArgument(mix_nAtomsDefinition, "k", "int32");
defineArgument(mix_nAtomsDefinition, "m", "int32");
defineOutput(mix_nAtomsDefinition, "RetVal", "double");
validate(mix_nAtomsDefinition);

%% C++ function |mix_nPhases| with MATLAB name |clib.ctMatlabInterface.mix_nPhases|
% C++ Signature: size_t mix_nPhases(int i)

mix_nPhasesDefinition = addFunction(libDef, ...
    "size_t mix_nPhases(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_nPhases", ...
    "Description", "clib.ctMatlabInterface.mix_nPhases Representation of C++ function mix_nPhases."); % Modify help description values as needed.
defineArgument(mix_nPhasesDefinition, "i", "int32");
defineOutput(mix_nPhasesDefinition, "RetVal", "uint64");
validate(mix_nPhasesDefinition);

%% C++ function |mix_phaseMoles| with MATLAB name |clib.ctMatlabInterface.mix_phaseMoles|
% C++ Signature: double mix_phaseMoles(int i,int n)

mix_phaseMolesDefinition = addFunction(libDef, ...
    "double mix_phaseMoles(int i,int n)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_phaseMoles", ...
    "Description", "clib.ctMatlabInterface.mix_phaseMoles Representation of C++ function mix_phaseMoles."); % Modify help description values as needed.
defineArgument(mix_phaseMolesDefinition, "i", "int32");
defineArgument(mix_phaseMolesDefinition, "n", "int32");
defineOutput(mix_phaseMolesDefinition, "RetVal", "double");
validate(mix_phaseMolesDefinition);

%% C++ function |mix_setPhaseMoles| with MATLAB name |clib.ctMatlabInterface.mix_setPhaseMoles|
% C++ Signature: int mix_setPhaseMoles(int i,int n,double v)

mix_setPhaseMolesDefinition = addFunction(libDef, ...
    "int mix_setPhaseMoles(int i,int n,double v)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_setPhaseMoles", ...
    "Description", "clib.ctMatlabInterface.mix_setPhaseMoles Representation of C++ function mix_setPhaseMoles."); % Modify help description values as needed.
defineArgument(mix_setPhaseMolesDefinition, "i", "int32");
defineArgument(mix_setPhaseMolesDefinition, "n", "int32");
defineArgument(mix_setPhaseMolesDefinition, "v", "double");
defineOutput(mix_setPhaseMolesDefinition, "RetVal", "int32");
validate(mix_setPhaseMolesDefinition);

%% C++ function |mix_setMoles| with MATLAB name |clib.ctMatlabInterface.mix_setMoles|
% C++ Signature: int mix_setMoles(int i,size_t nlen,double const * n)

mix_setMolesDefinition = addFunction(libDef, ...
    "int mix_setMoles(int i,size_t nlen,double const * n)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_setMoles", ...
    "Description", "clib.ctMatlabInterface.mix_setMoles Representation of C++ function mix_setMoles."); % Modify help description values as needed.
defineArgument(mix_setMolesDefinition, "i", "int32");
defineArgument(mix_setMolesDefinition, "nlen", "uint64");
defineArgument(mix_setMolesDefinition, "n", "clib.array.ctMatlabInterface.Double", "input", "nlen"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(mix_setMolesDefinition, "RetVal", "int32");
validate(mix_setMolesDefinition);

%% C++ function |mix_setMolesByName| with MATLAB name |clib.ctMatlabInterface.mix_setMolesByName|
% C++ Signature: int mix_setMolesByName(int i,char const * n)

mix_setMolesByNameDefinition = addFunction(libDef, ...
    "int mix_setMolesByName(int i,char const * n)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_setMolesByName", ...
    "Description", "clib.ctMatlabInterface.mix_setMolesByName Representation of C++ function mix_setMolesByName."); % Modify help description values as needed.
defineArgument(mix_setMolesByNameDefinition, "i", "int32");
defineArgument(mix_setMolesByNameDefinition, "n", "string", "input", "nullTerminated");
defineOutput(mix_setMolesByNameDefinition, "RetVal", "int32");
validate(mix_setMolesByNameDefinition);

%% C++ function |mix_speciesMoles| with MATLAB name |clib.ctMatlabInterface.mix_speciesMoles|
% C++ Signature: double mix_speciesMoles(int i,int k)

mix_speciesMolesDefinition = addFunction(libDef, ...
    "double mix_speciesMoles(int i,int k)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_speciesMoles", ...
    "Description", "clib.ctMatlabInterface.mix_speciesMoles Representation of C++ function mix_speciesMoles."); % Modify help description values as needed.
defineArgument(mix_speciesMolesDefinition, "i", "int32");
defineArgument(mix_speciesMolesDefinition, "k", "int32");
defineOutput(mix_speciesMolesDefinition, "RetVal", "double");
validate(mix_speciesMolesDefinition);

%% C++ function |mix_elementMoles| with MATLAB name |clib.ctMatlabInterface.mix_elementMoles|
% C++ Signature: double mix_elementMoles(int i,int m)

mix_elementMolesDefinition = addFunction(libDef, ...
    "double mix_elementMoles(int i,int m)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_elementMoles", ...
    "Description", "clib.ctMatlabInterface.mix_elementMoles Representation of C++ function mix_elementMoles."); % Modify help description values as needed.
defineArgument(mix_elementMolesDefinition, "i", "int32");
defineArgument(mix_elementMolesDefinition, "m", "int32");
defineOutput(mix_elementMolesDefinition, "RetVal", "double");
validate(mix_elementMolesDefinition);

%% C++ function |mix_equilibrate| with MATLAB name |clib.ctMatlabInterface.mix_equilibrate|
% C++ Signature: double mix_equilibrate(int i,char const * XY,double err,int maxsteps,int maxiter,int loglevel)

mix_equilibrateDefinition = addFunction(libDef, ...
    "double mix_equilibrate(int i,char const * XY,double err,int maxsteps,int maxiter,int loglevel)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_equilibrate", ...
    "Description", "clib.ctMatlabInterface.mix_equilibrate Representation of C++ function mix_equilibrate."); % Modify help description values as needed.
defineArgument(mix_equilibrateDefinition, "i", "int32");
defineArgument(mix_equilibrateDefinition, "XY", "string", "input", "nullTerminated");
defineArgument(mix_equilibrateDefinition, "err", "double");
defineArgument(mix_equilibrateDefinition, "maxsteps", "int32");
defineArgument(mix_equilibrateDefinition, "maxiter", "int32");
defineArgument(mix_equilibrateDefinition, "loglevel", "int32");
defineOutput(mix_equilibrateDefinition, "RetVal", "double");
validate(mix_equilibrateDefinition);

%% C++ function |mix_getChemPotentials| with MATLAB name |clib.ctMatlabInterface.mix_getChemPotentials|
% C++ Signature: int mix_getChemPotentials(int i,size_t lenmu,double * mu)

mix_getChemPotentialsDefinition = addFunction(libDef, ...
    "int mix_getChemPotentials(int i,size_t lenmu,double * mu)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_getChemPotentials", ...
    "Description", "clib.ctMatlabInterface.mix_getChemPotentials Representation of C++ function mix_getChemPotentials."); % Modify help description values as needed.
defineArgument(mix_getChemPotentialsDefinition, "i", "int32");
defineArgument(mix_getChemPotentialsDefinition, "lenmu", "uint64");
defineArgument(mix_getChemPotentialsDefinition, "mu", "clib.array.ctMatlabInterface.Double", "input", "lenmu"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(mix_getChemPotentialsDefinition, "RetVal", "int32");
validate(mix_getChemPotentialsDefinition);

%% C++ function |mix_enthalpy| with MATLAB name |clib.ctMatlabInterface.mix_enthalpy|
% C++ Signature: double mix_enthalpy(int i)

mix_enthalpyDefinition = addFunction(libDef, ...
    "double mix_enthalpy(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_enthalpy", ...
    "Description", "clib.ctMatlabInterface.mix_enthalpy Representation of C++ function mix_enthalpy."); % Modify help description values as needed.
defineArgument(mix_enthalpyDefinition, "i", "int32");
defineOutput(mix_enthalpyDefinition, "RetVal", "double");
validate(mix_enthalpyDefinition);

%% C++ function |mix_entropy| with MATLAB name |clib.ctMatlabInterface.mix_entropy|
% C++ Signature: double mix_entropy(int i)

mix_entropyDefinition = addFunction(libDef, ...
    "double mix_entropy(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_entropy", ...
    "Description", "clib.ctMatlabInterface.mix_entropy Representation of C++ function mix_entropy."); % Modify help description values as needed.
defineArgument(mix_entropyDefinition, "i", "int32");
defineOutput(mix_entropyDefinition, "RetVal", "double");
validate(mix_entropyDefinition);

%% C++ function |mix_gibbs| with MATLAB name |clib.ctMatlabInterface.mix_gibbs|
% C++ Signature: double mix_gibbs(int i)

mix_gibbsDefinition = addFunction(libDef, ...
    "double mix_gibbs(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_gibbs", ...
    "Description", "clib.ctMatlabInterface.mix_gibbs Representation of C++ function mix_gibbs."); % Modify help description values as needed.
defineArgument(mix_gibbsDefinition, "i", "int32");
defineOutput(mix_gibbsDefinition, "RetVal", "double");
validate(mix_gibbsDefinition);

%% C++ function |mix_cp| with MATLAB name |clib.ctMatlabInterface.mix_cp|
% C++ Signature: double mix_cp(int i)

mix_cpDefinition = addFunction(libDef, ...
    "double mix_cp(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_cp", ...
    "Description", "clib.ctMatlabInterface.mix_cp Representation of C++ function mix_cp."); % Modify help description values as needed.
defineArgument(mix_cpDefinition, "i", "int32");
defineOutput(mix_cpDefinition, "RetVal", "double");
validate(mix_cpDefinition);

%% C++ function |mix_volume| with MATLAB name |clib.ctMatlabInterface.mix_volume|
% C++ Signature: double mix_volume(int i)

mix_volumeDefinition = addFunction(libDef, ...
    "double mix_volume(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_volume", ...
    "Description", "clib.ctMatlabInterface.mix_volume Representation of C++ function mix_volume."); % Modify help description values as needed.
defineArgument(mix_volumeDefinition, "i", "int32");
defineOutput(mix_volumeDefinition, "RetVal", "double");
validate(mix_volumeDefinition);

%% C++ function |mix_speciesPhaseIndex| with MATLAB name |clib.ctMatlabInterface.mix_speciesPhaseIndex|
% C++ Signature: size_t mix_speciesPhaseIndex(int i,int k)

mix_speciesPhaseIndexDefinition = addFunction(libDef, ...
    "size_t mix_speciesPhaseIndex(int i,int k)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_speciesPhaseIndex", ...
    "Description", "clib.ctMatlabInterface.mix_speciesPhaseIndex Representation of C++ function mix_speciesPhaseIndex."); % Modify help description values as needed.
defineArgument(mix_speciesPhaseIndexDefinition, "i", "int32");
defineArgument(mix_speciesPhaseIndexDefinition, "k", "int32");
defineOutput(mix_speciesPhaseIndexDefinition, "RetVal", "uint64");
validate(mix_speciesPhaseIndexDefinition);

%% C++ function |mix_moleFraction| with MATLAB name |clib.ctMatlabInterface.mix_moleFraction|
% C++ Signature: double mix_moleFraction(int i,int k)

mix_moleFractionDefinition = addFunction(libDef, ...
    "double mix_moleFraction(int i,int k)", ...
    "MATLABName", "clib.ctMatlabInterface.mix_moleFraction", ...
    "Description", "clib.ctMatlabInterface.mix_moleFraction Representation of C++ function mix_moleFraction."); % Modify help description values as needed.
defineArgument(mix_moleFractionDefinition, "i", "int32");
defineArgument(mix_moleFractionDefinition, "k", "int32");
defineOutput(mix_moleFractionDefinition, "RetVal", "double");
validate(mix_moleFractionDefinition);

%% C++ function |ct_clearOneDim| with MATLAB name |clib.ctMatlabInterface.ct_clearOneDim|
% C++ Signature: int ct_clearOneDim()

ct_clearOneDimDefinition = addFunction(libDef, ...
    "int ct_clearOneDim()", ...
    "MATLABName", "clib.ctMatlabInterface.ct_clearOneDim", ...
    "Description", "clib.ctMatlabInterface.ct_clearOneDim Representation of C++ function ct_clearOneDim."); % Modify help description values as needed.
defineOutput(ct_clearOneDimDefinition, "RetVal", "int32");
validate(ct_clearOneDimDefinition);

%% C++ function |domain_new| with MATLAB name |clib.ctMatlabInterface.domain_new|
% C++ Signature: int domain_new(char const * type,int i,char const * id)

domain_newDefinition = addFunction(libDef, ...
    "int domain_new(char const * type,int i,char const * id)", ...
    "MATLABName", "clib.ctMatlabInterface.domain_new", ...
    "Description", "clib.ctMatlabInterface.domain_new Representation of C++ function domain_new."); % Modify help description values as needed.
defineArgument(domain_newDefinition, "type", "string", "input", "nullTerminated");
defineArgument(domain_newDefinition, "i", "int32");
defineArgument(domain_newDefinition, "id", "string", "input", "nullTerminated");
defineOutput(domain_newDefinition, "RetVal", "int32");
validate(domain_newDefinition);

%% C++ function |domain_del| with MATLAB name |clib.ctMatlabInterface.domain_del|
% C++ Signature: int domain_del(int i)

domain_delDefinition = addFunction(libDef, ...
    "int domain_del(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.domain_del", ...
    "Description", "clib.ctMatlabInterface.domain_del Representation of C++ function domain_del."); % Modify help description values as needed.
defineArgument(domain_delDefinition, "i", "int32");
defineOutput(domain_delDefinition, "RetVal", "int32");
validate(domain_delDefinition);

%% C++ function |domain_type| with MATLAB name |clib.ctMatlabInterface.domain_type|
% C++ Signature: int domain_type(int i,size_t lennm,char * nm)

domain_typeDefinition = addFunction(libDef, ...
    "int domain_type(int i,size_t lennm,char * nm)", ...
    "MATLABName", "clib.ctMatlabInterface.domain_type", ...
    "Description", "clib.ctMatlabInterface.domain_type Representation of C++ function domain_type."); % Modify help description values as needed.
defineArgument(domain_typeDefinition, "i", "int32");
defineArgument(domain_typeDefinition, "lennm", "uint64");
defineArgument(domain_typeDefinition, "nm", "clib.array.ctMatlabInterface.Char", "input", "lennm"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(domain_typeDefinition, "RetVal", "int32");
validate(domain_typeDefinition);

%% C++ function |domain_index| with MATLAB name |clib.ctMatlabInterface.domain_index|
% C++ Signature: size_t domain_index(int i)

domain_indexDefinition = addFunction(libDef, ...
    "size_t domain_index(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.domain_index", ...
    "Description", "clib.ctMatlabInterface.domain_index Representation of C++ function domain_index."); % Modify help description values as needed.
defineArgument(domain_indexDefinition, "i", "int32");
defineOutput(domain_indexDefinition, "RetVal", "uint64");
validate(domain_indexDefinition);

%% C++ function |domain_nComponents| with MATLAB name |clib.ctMatlabInterface.domain_nComponents|
% C++ Signature: size_t domain_nComponents(int i)

domain_nComponentsDefinition = addFunction(libDef, ...
    "size_t domain_nComponents(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.domain_nComponents", ...
    "Description", "clib.ctMatlabInterface.domain_nComponents Representation of C++ function domain_nComponents."); % Modify help description values as needed.
defineArgument(domain_nComponentsDefinition, "i", "int32");
defineOutput(domain_nComponentsDefinition, "RetVal", "uint64");
validate(domain_nComponentsDefinition);

%% C++ function |domain_nPoints| with MATLAB name |clib.ctMatlabInterface.domain_nPoints|
% C++ Signature: size_t domain_nPoints(int i)

domain_nPointsDefinition = addFunction(libDef, ...
    "size_t domain_nPoints(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.domain_nPoints", ...
    "Description", "clib.ctMatlabInterface.domain_nPoints Representation of C++ function domain_nPoints."); % Modify help description values as needed.
defineArgument(domain_nPointsDefinition, "i", "int32");
defineOutput(domain_nPointsDefinition, "RetVal", "uint64");
validate(domain_nPointsDefinition);

%% C++ function |domain_componentName| with MATLAB name |clib.ctMatlabInterface.domain_componentName|
% C++ Signature: int domain_componentName(int i,int n,int sz,char * nameout)

domain_componentNameDefinition = addFunction(libDef, ...
    "int domain_componentName(int i,int n,int sz,char * nameout)", ...
    "MATLABName", "clib.ctMatlabInterface.domain_componentName", ...
    "Description", "clib.ctMatlabInterface.domain_componentName Representation of C++ function domain_componentName."); % Modify help description values as needed.
defineArgument(domain_componentNameDefinition, "i", "int32");
defineArgument(domain_componentNameDefinition, "n", "int32");
defineArgument(domain_componentNameDefinition, "sz", "int32");
defineArgument(domain_componentNameDefinition, "nameout", "clib.array.ctMatlabInterface.Char", "input", "sz"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(domain_componentNameDefinition, "RetVal", "int32");
validate(domain_componentNameDefinition);

%% C++ function |domain_componentIndex| with MATLAB name |clib.ctMatlabInterface.domain_componentIndex|
% C++ Signature: size_t domain_componentIndex(int i,char const * name)

domain_componentIndexDefinition = addFunction(libDef, ...
    "size_t domain_componentIndex(int i,char const * name)", ...
    "MATLABName", "clib.ctMatlabInterface.domain_componentIndex", ...
    "Description", "clib.ctMatlabInterface.domain_componentIndex Representation of C++ function domain_componentIndex."); % Modify help description values as needed.
defineArgument(domain_componentIndexDefinition, "i", "int32");
defineArgument(domain_componentIndexDefinition, "name", "string", "input", "nullTerminated");
defineOutput(domain_componentIndexDefinition, "RetVal", "uint64");
validate(domain_componentIndexDefinition);

%% C++ function |domain_setBounds| with MATLAB name |clib.ctMatlabInterface.domain_setBounds|
% C++ Signature: int domain_setBounds(int i,int n,double lower,double upper)

domain_setBoundsDefinition = addFunction(libDef, ...
    "int domain_setBounds(int i,int n,double lower,double upper)", ...
    "MATLABName", "clib.ctMatlabInterface.domain_setBounds", ...
    "Description", "clib.ctMatlabInterface.domain_setBounds Representation of C++ function domain_setBounds."); % Modify help description values as needed.
defineArgument(domain_setBoundsDefinition, "i", "int32");
defineArgument(domain_setBoundsDefinition, "n", "int32");
defineArgument(domain_setBoundsDefinition, "lower", "double");
defineArgument(domain_setBoundsDefinition, "upper", "double");
defineOutput(domain_setBoundsDefinition, "RetVal", "int32");
validate(domain_setBoundsDefinition);

%% C++ function |domain_lowerBound| with MATLAB name |clib.ctMatlabInterface.domain_lowerBound|
% C++ Signature: double domain_lowerBound(int i,int n)

domain_lowerBoundDefinition = addFunction(libDef, ...
    "double domain_lowerBound(int i,int n)", ...
    "MATLABName", "clib.ctMatlabInterface.domain_lowerBound", ...
    "Description", "clib.ctMatlabInterface.domain_lowerBound Representation of C++ function domain_lowerBound."); % Modify help description values as needed.
defineArgument(domain_lowerBoundDefinition, "i", "int32");
defineArgument(domain_lowerBoundDefinition, "n", "int32");
defineOutput(domain_lowerBoundDefinition, "RetVal", "double");
validate(domain_lowerBoundDefinition);

%% C++ function |domain_upperBound| with MATLAB name |clib.ctMatlabInterface.domain_upperBound|
% C++ Signature: double domain_upperBound(int i,int n)

domain_upperBoundDefinition = addFunction(libDef, ...
    "double domain_upperBound(int i,int n)", ...
    "MATLABName", "clib.ctMatlabInterface.domain_upperBound", ...
    "Description", "clib.ctMatlabInterface.domain_upperBound Representation of C++ function domain_upperBound."); % Modify help description values as needed.
defineArgument(domain_upperBoundDefinition, "i", "int32");
defineArgument(domain_upperBoundDefinition, "n", "int32");
defineOutput(domain_upperBoundDefinition, "RetVal", "double");
validate(domain_upperBoundDefinition);

%% C++ function |domain_setSteadyTolerances| with MATLAB name |clib.ctMatlabInterface.domain_setSteadyTolerances|
% C++ Signature: int domain_setSteadyTolerances(int i,int n,double rtol,double atol)

domain_setSteadyTolerancesDefinition = addFunction(libDef, ...
    "int domain_setSteadyTolerances(int i,int n,double rtol,double atol)", ...
    "MATLABName", "clib.ctMatlabInterface.domain_setSteadyTolerances", ...
    "Description", "clib.ctMatlabInterface.domain_setSteadyTolerances Representation of C++ function domain_setSteadyTolerances."); % Modify help description values as needed.
defineArgument(domain_setSteadyTolerancesDefinition, "i", "int32");
defineArgument(domain_setSteadyTolerancesDefinition, "n", "int32");
defineArgument(domain_setSteadyTolerancesDefinition, "rtol", "double");
defineArgument(domain_setSteadyTolerancesDefinition, "atol", "double");
defineOutput(domain_setSteadyTolerancesDefinition, "RetVal", "int32");
validate(domain_setSteadyTolerancesDefinition);

%% C++ function |domain_setTransientTolerances| with MATLAB name |clib.ctMatlabInterface.domain_setTransientTolerances|
% C++ Signature: int domain_setTransientTolerances(int i,int n,double rtol,double atol)

domain_setTransientTolerancesDefinition = addFunction(libDef, ...
    "int domain_setTransientTolerances(int i,int n,double rtol,double atol)", ...
    "MATLABName", "clib.ctMatlabInterface.domain_setTransientTolerances", ...
    "Description", "clib.ctMatlabInterface.domain_setTransientTolerances Representation of C++ function domain_setTransientTolerances."); % Modify help description values as needed.
defineArgument(domain_setTransientTolerancesDefinition, "i", "int32");
defineArgument(domain_setTransientTolerancesDefinition, "n", "int32");
defineArgument(domain_setTransientTolerancesDefinition, "rtol", "double");
defineArgument(domain_setTransientTolerancesDefinition, "atol", "double");
defineOutput(domain_setTransientTolerancesDefinition, "RetVal", "int32");
validate(domain_setTransientTolerancesDefinition);

%% C++ function |domain_rtol| with MATLAB name |clib.ctMatlabInterface.domain_rtol|
% C++ Signature: double domain_rtol(int i,int n)

domain_rtolDefinition = addFunction(libDef, ...
    "double domain_rtol(int i,int n)", ...
    "MATLABName", "clib.ctMatlabInterface.domain_rtol", ...
    "Description", "clib.ctMatlabInterface.domain_rtol Representation of C++ function domain_rtol."); % Modify help description values as needed.
defineArgument(domain_rtolDefinition, "i", "int32");
defineArgument(domain_rtolDefinition, "n", "int32");
defineOutput(domain_rtolDefinition, "RetVal", "double");
validate(domain_rtolDefinition);

%% C++ function |domain_atol| with MATLAB name |clib.ctMatlabInterface.domain_atol|
% C++ Signature: double domain_atol(int i,int n)

domain_atolDefinition = addFunction(libDef, ...
    "double domain_atol(int i,int n)", ...
    "MATLABName", "clib.ctMatlabInterface.domain_atol", ...
    "Description", "clib.ctMatlabInterface.domain_atol Representation of C++ function domain_atol."); % Modify help description values as needed.
defineArgument(domain_atolDefinition, "i", "int32");
defineArgument(domain_atolDefinition, "n", "int32");
defineOutput(domain_atolDefinition, "RetVal", "double");
validate(domain_atolDefinition);

%% C++ function |domain_setupGrid| with MATLAB name |clib.ctMatlabInterface.domain_setupGrid|
% C++ Signature: int domain_setupGrid(int i,size_t npts,double const * grid)

domain_setupGridDefinition = addFunction(libDef, ...
    "int domain_setupGrid(int i,size_t npts,double const * grid)", ...
    "MATLABName", "clib.ctMatlabInterface.domain_setupGrid", ...
    "Description", "clib.ctMatlabInterface.domain_setupGrid Representation of C++ function domain_setupGrid."); % Modify help description values as needed.
defineArgument(domain_setupGridDefinition, "i", "int32");
defineArgument(domain_setupGridDefinition, "npts", "uint64");
defineArgument(domain_setupGridDefinition, "grid", "clib.array.ctMatlabInterface.Double", "input", "npts"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(domain_setupGridDefinition, "RetVal", "int32");
validate(domain_setupGridDefinition);

%% C++ function |domain_setID| with MATLAB name |clib.ctMatlabInterface.domain_setID|
% C++ Signature: int domain_setID(int i,char const * id)

domain_setIDDefinition = addFunction(libDef, ...
    "int domain_setID(int i,char const * id)", ...
    "MATLABName", "clib.ctMatlabInterface.domain_setID", ...
    "Description", "clib.ctMatlabInterface.domain_setID Representation of C++ function domain_setID."); % Modify help description values as needed.
defineArgument(domain_setIDDefinition, "i", "int32");
defineArgument(domain_setIDDefinition, "id", "string", "input", "nullTerminated");
defineOutput(domain_setIDDefinition, "RetVal", "int32");
validate(domain_setIDDefinition);

%% C++ function |domain_grid| with MATLAB name |clib.ctMatlabInterface.domain_grid|
% C++ Signature: double domain_grid(int i,int n)

domain_gridDefinition = addFunction(libDef, ...
    "double domain_grid(int i,int n)", ...
    "MATLABName", "clib.ctMatlabInterface.domain_grid", ...
    "Description", "clib.ctMatlabInterface.domain_grid Representation of C++ function domain_grid."); % Modify help description values as needed.
defineArgument(domain_gridDefinition, "i", "int32");
defineArgument(domain_gridDefinition, "n", "int32");
defineOutput(domain_gridDefinition, "RetVal", "double");
validate(domain_gridDefinition);

%% C++ function |bdry_setMdot| with MATLAB name |clib.ctMatlabInterface.bdry_setMdot|
% C++ Signature: int bdry_setMdot(int i,double mdot)

bdry_setMdotDefinition = addFunction(libDef, ...
    "int bdry_setMdot(int i,double mdot)", ...
    "MATLABName", "clib.ctMatlabInterface.bdry_setMdot", ...
    "Description", "clib.ctMatlabInterface.bdry_setMdot Representation of C++ function bdry_setMdot."); % Modify help description values as needed.
defineArgument(bdry_setMdotDefinition, "i", "int32");
defineArgument(bdry_setMdotDefinition, "mdot", "double");
defineOutput(bdry_setMdotDefinition, "RetVal", "int32");
validate(bdry_setMdotDefinition);

%% C++ function |bdry_setTemperature| with MATLAB name |clib.ctMatlabInterface.bdry_setTemperature|
% C++ Signature: int bdry_setTemperature(int i,double t)

bdry_setTemperatureDefinition = addFunction(libDef, ...
    "int bdry_setTemperature(int i,double t)", ...
    "MATLABName", "clib.ctMatlabInterface.bdry_setTemperature", ...
    "Description", "clib.ctMatlabInterface.bdry_setTemperature Representation of C++ function bdry_setTemperature."); % Modify help description values as needed.
defineArgument(bdry_setTemperatureDefinition, "i", "int32");
defineArgument(bdry_setTemperatureDefinition, "t", "double");
defineOutput(bdry_setTemperatureDefinition, "RetVal", "int32");
validate(bdry_setTemperatureDefinition);

%% C++ function |bdry_setSpreadRate| with MATLAB name |clib.ctMatlabInterface.bdry_setSpreadRate|
% C++ Signature: int bdry_setSpreadRate(int i,double v)

bdry_setSpreadRateDefinition = addFunction(libDef, ...
    "int bdry_setSpreadRate(int i,double v)", ...
    "MATLABName", "clib.ctMatlabInterface.bdry_setSpreadRate", ...
    "Description", "clib.ctMatlabInterface.bdry_setSpreadRate Representation of C++ function bdry_setSpreadRate."); % Modify help description values as needed.
defineArgument(bdry_setSpreadRateDefinition, "i", "int32");
defineArgument(bdry_setSpreadRateDefinition, "v", "double");
defineOutput(bdry_setSpreadRateDefinition, "RetVal", "int32");
validate(bdry_setSpreadRateDefinition);

%% C++ function |bdry_setMoleFractions| with MATLAB name |clib.ctMatlabInterface.bdry_setMoleFractions|
% C++ Signature: int bdry_setMoleFractions(int i,char const * x)

bdry_setMoleFractionsDefinition = addFunction(libDef, ...
    "int bdry_setMoleFractions(int i,char const * x)", ...
    "MATLABName", "clib.ctMatlabInterface.bdry_setMoleFractions", ...
    "Description", "clib.ctMatlabInterface.bdry_setMoleFractions Representation of C++ function bdry_setMoleFractions."); % Modify help description values as needed.
defineArgument(bdry_setMoleFractionsDefinition, "i", "int32");
defineArgument(bdry_setMoleFractionsDefinition, "x", "string", "input", "nullTerminated");
defineOutput(bdry_setMoleFractionsDefinition, "RetVal", "int32");
validate(bdry_setMoleFractionsDefinition);

%% C++ function |bdry_temperature| with MATLAB name |clib.ctMatlabInterface.bdry_temperature|
% C++ Signature: double bdry_temperature(int i)

bdry_temperatureDefinition = addFunction(libDef, ...
    "double bdry_temperature(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.bdry_temperature", ...
    "Description", "clib.ctMatlabInterface.bdry_temperature Representation of C++ function bdry_temperature."); % Modify help description values as needed.
defineArgument(bdry_temperatureDefinition, "i", "int32");
defineOutput(bdry_temperatureDefinition, "RetVal", "double");
validate(bdry_temperatureDefinition);

%% C++ function |bdry_spreadRate| with MATLAB name |clib.ctMatlabInterface.bdry_spreadRate|
% C++ Signature: double bdry_spreadRate(int i)

bdry_spreadRateDefinition = addFunction(libDef, ...
    "double bdry_spreadRate(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.bdry_spreadRate", ...
    "Description", "clib.ctMatlabInterface.bdry_spreadRate Representation of C++ function bdry_spreadRate."); % Modify help description values as needed.
defineArgument(bdry_spreadRateDefinition, "i", "int32");
defineOutput(bdry_spreadRateDefinition, "RetVal", "double");
validate(bdry_spreadRateDefinition);

%% C++ function |bdry_massFraction| with MATLAB name |clib.ctMatlabInterface.bdry_massFraction|
% C++ Signature: double bdry_massFraction(int i,int k)

bdry_massFractionDefinition = addFunction(libDef, ...
    "double bdry_massFraction(int i,int k)", ...
    "MATLABName", "clib.ctMatlabInterface.bdry_massFraction", ...
    "Description", "clib.ctMatlabInterface.bdry_massFraction Representation of C++ function bdry_massFraction."); % Modify help description values as needed.
defineArgument(bdry_massFractionDefinition, "i", "int32");
defineArgument(bdry_massFractionDefinition, "k", "int32");
defineOutput(bdry_massFractionDefinition, "RetVal", "double");
validate(bdry_massFractionDefinition);

%% C++ function |bdry_mdot| with MATLAB name |clib.ctMatlabInterface.bdry_mdot|
% C++ Signature: double bdry_mdot(int i)

bdry_mdotDefinition = addFunction(libDef, ...
    "double bdry_mdot(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.bdry_mdot", ...
    "Description", "clib.ctMatlabInterface.bdry_mdot Representation of C++ function bdry_mdot."); % Modify help description values as needed.
defineArgument(bdry_mdotDefinition, "i", "int32");
defineOutput(bdry_mdotDefinition, "RetVal", "double");
validate(bdry_mdotDefinition);

%% C++ function |reactingsurf_enableCoverageEqs| with MATLAB name |clib.ctMatlabInterface.reactingsurf_enableCoverageEqs|
% C++ Signature: int reactingsurf_enableCoverageEqs(int i,int onoff)

reactingsurf_enableCoverageEqsDefinition = addFunction(libDef, ...
    "int reactingsurf_enableCoverageEqs(int i,int onoff)", ...
    "MATLABName", "clib.ctMatlabInterface.reactingsurf_enableCoverageEqs", ...
    "Description", "clib.ctMatlabInterface.reactingsurf_enableCoverageEqs Representation of C++ function reactingsurf_enableCoverageEqs."); % Modify help description values as needed.
defineArgument(reactingsurf_enableCoverageEqsDefinition, "i", "int32");
defineArgument(reactingsurf_enableCoverageEqsDefinition, "onoff", "int32");
defineOutput(reactingsurf_enableCoverageEqsDefinition, "RetVal", "int32");
validate(reactingsurf_enableCoverageEqsDefinition);

%% C++ function |flow1D_setTransport| with MATLAB name |clib.ctMatlabInterface.flow1D_setTransport|
% C++ Signature: int flow1D_setTransport(int i,int itr)

flow1D_setTransportDefinition = addFunction(libDef, ...
    "int flow1D_setTransport(int i,int itr)", ...
    "MATLABName", "clib.ctMatlabInterface.flow1D_setTransport", ...
    "Description", "clib.ctMatlabInterface.flow1D_setTransport Representation of C++ function flow1D_setTransport."); % Modify help description values as needed.
defineArgument(flow1D_setTransportDefinition, "i", "int32");
defineArgument(flow1D_setTransportDefinition, "itr", "int32");
defineOutput(flow1D_setTransportDefinition, "RetVal", "int32");
validate(flow1D_setTransportDefinition);

%% C++ function |flow1D_enableSoret| with MATLAB name |clib.ctMatlabInterface.flow1D_enableSoret|
% C++ Signature: int flow1D_enableSoret(int i,int iSoret)

flow1D_enableSoretDefinition = addFunction(libDef, ...
    "int flow1D_enableSoret(int i,int iSoret)", ...
    "MATLABName", "clib.ctMatlabInterface.flow1D_enableSoret", ...
    "Description", "clib.ctMatlabInterface.flow1D_enableSoret Representation of C++ function flow1D_enableSoret."); % Modify help description values as needed.
defineArgument(flow1D_enableSoretDefinition, "i", "int32");
defineArgument(flow1D_enableSoretDefinition, "iSoret", "int32");
defineOutput(flow1D_enableSoretDefinition, "RetVal", "int32");
validate(flow1D_enableSoretDefinition);

%% C++ function |flow1D_setPressure| with MATLAB name |clib.ctMatlabInterface.flow1D_setPressure|
% C++ Signature: int flow1D_setPressure(int i,double p)

flow1D_setPressureDefinition = addFunction(libDef, ...
    "int flow1D_setPressure(int i,double p)", ...
    "MATLABName", "clib.ctMatlabInterface.flow1D_setPressure", ...
    "Description", "clib.ctMatlabInterface.flow1D_setPressure Representation of C++ function flow1D_setPressure."); % Modify help description values as needed.
defineArgument(flow1D_setPressureDefinition, "i", "int32");
defineArgument(flow1D_setPressureDefinition, "p", "double");
defineOutput(flow1D_setPressureDefinition, "RetVal", "int32");
validate(flow1D_setPressureDefinition);

%% C++ function |flow1D_pressure| with MATLAB name |clib.ctMatlabInterface.flow1D_pressure|
% C++ Signature: double flow1D_pressure(int i)

flow1D_pressureDefinition = addFunction(libDef, ...
    "double flow1D_pressure(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.flow1D_pressure", ...
    "Description", "clib.ctMatlabInterface.flow1D_pressure Representation of C++ function flow1D_pressure."); % Modify help description values as needed.
defineArgument(flow1D_pressureDefinition, "i", "int32");
defineOutput(flow1D_pressureDefinition, "RetVal", "double");
validate(flow1D_pressureDefinition);

%% C++ function |flow1D_setFixedTempProfile| with MATLAB name |clib.ctMatlabInterface.flow1D_setFixedTempProfile|
% C++ Signature: int flow1D_setFixedTempProfile(int i,size_t n,double const * pos,size_t m,double const * temp)

flow1D_setFixedTempProfileDefinition = addFunction(libDef, ...
    "int flow1D_setFixedTempProfile(int i,size_t n,double const * pos,size_t m,double const * temp)", ...
    "MATLABName", "clib.ctMatlabInterface.flow1D_setFixedTempProfile", ...
    "Description", "clib.ctMatlabInterface.flow1D_setFixedTempProfile Representation of C++ function flow1D_setFixedTempProfile."); % Modify help description values as needed.
defineArgument(flow1D_setFixedTempProfileDefinition, "i", "int32");
defineArgument(flow1D_setFixedTempProfileDefinition, "n", "uint64");
defineArgument(flow1D_setFixedTempProfileDefinition, "pos", "clib.array.ctMatlabInterface.Double", "input", "n"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineArgument(flow1D_setFixedTempProfileDefinition, "m", "uint64");
defineArgument(flow1D_setFixedTempProfileDefinition, "temp", "clib.array.ctMatlabInterface.Double", "input", "m"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(flow1D_setFixedTempProfileDefinition, "RetVal", "int32");
validate(flow1D_setFixedTempProfileDefinition);

%% C++ function |flow1D_solveEnergyEqn| with MATLAB name |clib.ctMatlabInterface.flow1D_solveEnergyEqn|
% C++ Signature: int flow1D_solveEnergyEqn(int i,int flag)

flow1D_solveEnergyEqnDefinition = addFunction(libDef, ...
    "int flow1D_solveEnergyEqn(int i,int flag)", ...
    "MATLABName", "clib.ctMatlabInterface.flow1D_solveEnergyEqn", ...
    "Description", "clib.ctMatlabInterface.flow1D_solveEnergyEqn Representation of C++ function flow1D_solveEnergyEqn."); % Modify help description values as needed.
defineArgument(flow1D_solveEnergyEqnDefinition, "i", "int32");
defineArgument(flow1D_solveEnergyEqnDefinition, "flag", "int32");
defineOutput(flow1D_solveEnergyEqnDefinition, "RetVal", "int32");
validate(flow1D_solveEnergyEqnDefinition);

%% C++ function |sim1D_new| with MATLAB name |clib.ctMatlabInterface.sim1D_new|
% C++ Signature: int sim1D_new(size_t nd,int const * domains)

sim1D_newDefinition = addFunction(libDef, ...
    "int sim1D_new(size_t nd,int const * domains)", ...
    "MATLABName", "clib.ctMatlabInterface.sim1D_new", ...
    "Description", "clib.ctMatlabInterface.sim1D_new Representation of C++ function sim1D_new."); % Modify help description values as needed.
defineArgument(sim1D_newDefinition, "nd", "uint64");
defineArgument(sim1D_newDefinition, "domains", "clib.array.ctMatlabInterface.Int", "input", "nd"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Int", or "int32"
defineOutput(sim1D_newDefinition, "RetVal", "int32");
validate(sim1D_newDefinition);

%% C++ function |sim1D_del| with MATLAB name |clib.ctMatlabInterface.sim1D_del|
% C++ Signature: int sim1D_del(int i)

sim1D_delDefinition = addFunction(libDef, ...
    "int sim1D_del(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.sim1D_del", ...
    "Description", "clib.ctMatlabInterface.sim1D_del Representation of C++ function sim1D_del."); % Modify help description values as needed.
defineArgument(sim1D_delDefinition, "i", "int32");
defineOutput(sim1D_delDefinition, "RetVal", "int32");
validate(sim1D_delDefinition);

%% C++ function |sim1D_setValue| with MATLAB name |clib.ctMatlabInterface.sim1D_setValue|
% C++ Signature: int sim1D_setValue(int i,int dom,int comp,int localPoint,double value)

sim1D_setValueDefinition = addFunction(libDef, ...
    "int sim1D_setValue(int i,int dom,int comp,int localPoint,double value)", ...
    "MATLABName", "clib.ctMatlabInterface.sim1D_setValue", ...
    "Description", "clib.ctMatlabInterface.sim1D_setValue Representation of C++ function sim1D_setValue."); % Modify help description values as needed.
defineArgument(sim1D_setValueDefinition, "i", "int32");
defineArgument(sim1D_setValueDefinition, "dom", "int32");
defineArgument(sim1D_setValueDefinition, "comp", "int32");
defineArgument(sim1D_setValueDefinition, "localPoint", "int32");
defineArgument(sim1D_setValueDefinition, "value", "double");
defineOutput(sim1D_setValueDefinition, "RetVal", "int32");
validate(sim1D_setValueDefinition);

%% C++ function |sim1D_setProfile| with MATLAB name |clib.ctMatlabInterface.sim1D_setProfile|
% C++ Signature: int sim1D_setProfile(int i,int dom,int comp,size_t np,double const * pos,size_t nv,double const * v)

sim1D_setProfileDefinition = addFunction(libDef, ...
    "int sim1D_setProfile(int i,int dom,int comp,size_t np,double const * pos,size_t nv,double const * v)", ...
    "MATLABName", "clib.ctMatlabInterface.sim1D_setProfile", ...
    "Description", "clib.ctMatlabInterface.sim1D_setProfile Representation of C++ function sim1D_setProfile."); % Modify help description values as needed.
defineArgument(sim1D_setProfileDefinition, "i", "int32");
defineArgument(sim1D_setProfileDefinition, "dom", "int32");
defineArgument(sim1D_setProfileDefinition, "comp", "int32");
defineArgument(sim1D_setProfileDefinition, "np", "uint64");
defineArgument(sim1D_setProfileDefinition, "pos", "clib.array.ctMatlabInterface.Double", "input", "np"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineArgument(sim1D_setProfileDefinition, "nv", "uint64");
defineArgument(sim1D_setProfileDefinition, "v", "clib.array.ctMatlabInterface.Double", "input", "nv"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(sim1D_setProfileDefinition, "RetVal", "int32");
validate(sim1D_setProfileDefinition);

%% C++ function |sim1D_setFlatProfile| with MATLAB name |clib.ctMatlabInterface.sim1D_setFlatProfile|
% C++ Signature: int sim1D_setFlatProfile(int i,int dom,int comp,double v)

sim1D_setFlatProfileDefinition = addFunction(libDef, ...
    "int sim1D_setFlatProfile(int i,int dom,int comp,double v)", ...
    "MATLABName", "clib.ctMatlabInterface.sim1D_setFlatProfile", ...
    "Description", "clib.ctMatlabInterface.sim1D_setFlatProfile Representation of C++ function sim1D_setFlatProfile."); % Modify help description values as needed.
defineArgument(sim1D_setFlatProfileDefinition, "i", "int32");
defineArgument(sim1D_setFlatProfileDefinition, "dom", "int32");
defineArgument(sim1D_setFlatProfileDefinition, "comp", "int32");
defineArgument(sim1D_setFlatProfileDefinition, "v", "double");
defineOutput(sim1D_setFlatProfileDefinition, "RetVal", "int32");
validate(sim1D_setFlatProfileDefinition);

%% C++ function |sim1D_show| with MATLAB name |clib.ctMatlabInterface.sim1D_show|
% C++ Signature: int sim1D_show(int i)

sim1D_showDefinition = addFunction(libDef, ...
    "int sim1D_show(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.sim1D_show", ...
    "Description", "clib.ctMatlabInterface.sim1D_show Representation of C++ function sim1D_show."); % Modify help description values as needed.
defineArgument(sim1D_showDefinition, "i", "int32");
defineOutput(sim1D_showDefinition, "RetVal", "int32");
validate(sim1D_showDefinition);

%% C++ function |sim1D_setTimeStep| with MATLAB name |clib.ctMatlabInterface.sim1D_setTimeStep|
% C++ Signature: int sim1D_setTimeStep(int i,double stepsize,size_t ns,int const * nsteps)

sim1D_setTimeStepDefinition = addFunction(libDef, ...
    "int sim1D_setTimeStep(int i,double stepsize,size_t ns,int const * nsteps)", ...
    "MATLABName", "clib.ctMatlabInterface.sim1D_setTimeStep", ...
    "Description", "clib.ctMatlabInterface.sim1D_setTimeStep Representation of C++ function sim1D_setTimeStep."); % Modify help description values as needed.
defineArgument(sim1D_setTimeStepDefinition, "i", "int32");
defineArgument(sim1D_setTimeStepDefinition, "stepsize", "double");
defineArgument(sim1D_setTimeStepDefinition, "ns", "uint64");
defineArgument(sim1D_setTimeStepDefinition, "nsteps", "clib.array.ctMatlabInterface.Int", "input", "ns"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Int", or "int32"
defineOutput(sim1D_setTimeStepDefinition, "RetVal", "int32");
validate(sim1D_setTimeStepDefinition);

%% C++ function |sim1D_getInitialSoln| with MATLAB name |clib.ctMatlabInterface.sim1D_getInitialSoln|
% C++ Signature: int sim1D_getInitialSoln(int i)

sim1D_getInitialSolnDefinition = addFunction(libDef, ...
    "int sim1D_getInitialSoln(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.sim1D_getInitialSoln", ...
    "Description", "clib.ctMatlabInterface.sim1D_getInitialSoln Representation of C++ function sim1D_getInitialSoln."); % Modify help description values as needed.
defineArgument(sim1D_getInitialSolnDefinition, "i", "int32");
defineOutput(sim1D_getInitialSolnDefinition, "RetVal", "int32");
validate(sim1D_getInitialSolnDefinition);

%% C++ function |sim1D_solve| with MATLAB name |clib.ctMatlabInterface.sim1D_solve|
% C++ Signature: int sim1D_solve(int i,int loglevel,int refine_grid)

sim1D_solveDefinition = addFunction(libDef, ...
    "int sim1D_solve(int i,int loglevel,int refine_grid)", ...
    "MATLABName", "clib.ctMatlabInterface.sim1D_solve", ...
    "Description", "clib.ctMatlabInterface.sim1D_solve Representation of C++ function sim1D_solve."); % Modify help description values as needed.
defineArgument(sim1D_solveDefinition, "i", "int32");
defineArgument(sim1D_solveDefinition, "loglevel", "int32");
defineArgument(sim1D_solveDefinition, "refine_grid", "int32");
defineOutput(sim1D_solveDefinition, "RetVal", "int32");
validate(sim1D_solveDefinition);

%% C++ function |sim1D_refine| with MATLAB name |clib.ctMatlabInterface.sim1D_refine|
% C++ Signature: int sim1D_refine(int i,int loglevel)

sim1D_refineDefinition = addFunction(libDef, ...
    "int sim1D_refine(int i,int loglevel)", ...
    "MATLABName", "clib.ctMatlabInterface.sim1D_refine", ...
    "Description", "clib.ctMatlabInterface.sim1D_refine Representation of C++ function sim1D_refine."); % Modify help description values as needed.
defineArgument(sim1D_refineDefinition, "i", "int32");
defineArgument(sim1D_refineDefinition, "loglevel", "int32");
defineOutput(sim1D_refineDefinition, "RetVal", "int32");
validate(sim1D_refineDefinition);

%% C++ function |sim1D_setRefineCriteria| with MATLAB name |clib.ctMatlabInterface.sim1D_setRefineCriteria|
% C++ Signature: int sim1D_setRefineCriteria(int i,int dom,double ratio,double slope,double curve,double prune)

sim1D_setRefineCriteriaDefinition = addFunction(libDef, ...
    "int sim1D_setRefineCriteria(int i,int dom,double ratio,double slope,double curve,double prune)", ...
    "MATLABName", "clib.ctMatlabInterface.sim1D_setRefineCriteria", ...
    "Description", "clib.ctMatlabInterface.sim1D_setRefineCriteria Representation of C++ function sim1D_setRefineCriteria."); % Modify help description values as needed.
defineArgument(sim1D_setRefineCriteriaDefinition, "i", "int32");
defineArgument(sim1D_setRefineCriteriaDefinition, "dom", "int32");
defineArgument(sim1D_setRefineCriteriaDefinition, "ratio", "double");
defineArgument(sim1D_setRefineCriteriaDefinition, "slope", "double");
defineArgument(sim1D_setRefineCriteriaDefinition, "curve", "double");
defineArgument(sim1D_setRefineCriteriaDefinition, "prune", "double");
defineOutput(sim1D_setRefineCriteriaDefinition, "RetVal", "int32");
validate(sim1D_setRefineCriteriaDefinition);

%% C++ function |sim1D_setGridMin| with MATLAB name |clib.ctMatlabInterface.sim1D_setGridMin|
% C++ Signature: int sim1D_setGridMin(int i,int dom,double gridmin)

sim1D_setGridMinDefinition = addFunction(libDef, ...
    "int sim1D_setGridMin(int i,int dom,double gridmin)", ...
    "MATLABName", "clib.ctMatlabInterface.sim1D_setGridMin", ...
    "Description", "clib.ctMatlabInterface.sim1D_setGridMin Representation of C++ function sim1D_setGridMin."); % Modify help description values as needed.
defineArgument(sim1D_setGridMinDefinition, "i", "int32");
defineArgument(sim1D_setGridMinDefinition, "dom", "int32");
defineArgument(sim1D_setGridMinDefinition, "gridmin", "double");
defineOutput(sim1D_setGridMinDefinition, "RetVal", "int32");
validate(sim1D_setGridMinDefinition);

%% C++ function |sim1D_save| with MATLAB name |clib.ctMatlabInterface.sim1D_save|
% C++ Signature: int sim1D_save(int i,char const * fname,char const * id,char const * desc)

sim1D_saveDefinition = addFunction(libDef, ...
    "int sim1D_save(int i,char const * fname,char const * id,char const * desc)", ...
    "MATLABName", "clib.ctMatlabInterface.sim1D_save", ...
    "Description", "clib.ctMatlabInterface.sim1D_save Representation of C++ function sim1D_save."); % Modify help description values as needed.
defineArgument(sim1D_saveDefinition, "i", "int32");
defineArgument(sim1D_saveDefinition, "fname", "string", "input", "nullTerminated");
defineArgument(sim1D_saveDefinition, "id", "string", "input", "nullTerminated");
defineArgument(sim1D_saveDefinition, "desc", "string", "input", "nullTerminated");
defineOutput(sim1D_saveDefinition, "RetVal", "int32");
validate(sim1D_saveDefinition);

%% C++ function |sim1D_restore| with MATLAB name |clib.ctMatlabInterface.sim1D_restore|
% C++ Signature: int sim1D_restore(int i,char const * fname,char const * id)

sim1D_restoreDefinition = addFunction(libDef, ...
    "int sim1D_restore(int i,char const * fname,char const * id)", ...
    "MATLABName", "clib.ctMatlabInterface.sim1D_restore", ...
    "Description", "clib.ctMatlabInterface.sim1D_restore Representation of C++ function sim1D_restore."); % Modify help description values as needed.
defineArgument(sim1D_restoreDefinition, "i", "int32");
defineArgument(sim1D_restoreDefinition, "fname", "string", "input", "nullTerminated");
defineArgument(sim1D_restoreDefinition, "id", "string", "input", "nullTerminated");
defineOutput(sim1D_restoreDefinition, "RetVal", "int32");
validate(sim1D_restoreDefinition);

%% C++ function |sim1D_writeStats| with MATLAB name |clib.ctMatlabInterface.sim1D_writeStats|
% C++ Signature: int sim1D_writeStats(int i,int printTime)

sim1D_writeStatsDefinition = addFunction(libDef, ...
    "int sim1D_writeStats(int i,int printTime)", ...
    "MATLABName", "clib.ctMatlabInterface.sim1D_writeStats", ...
    "Description", "clib.ctMatlabInterface.sim1D_writeStats Representation of C++ function sim1D_writeStats."); % Modify help description values as needed.
defineArgument(sim1D_writeStatsDefinition, "i", "int32");
defineArgument(sim1D_writeStatsDefinition, "printTime", "int32");
defineOutput(sim1D_writeStatsDefinition, "RetVal", "int32");
validate(sim1D_writeStatsDefinition);

%% C++ function |sim1D_domainIndex| with MATLAB name |clib.ctMatlabInterface.sim1D_domainIndex|
% C++ Signature: int sim1D_domainIndex(int i,char const * name)

sim1D_domainIndexDefinition = addFunction(libDef, ...
    "int sim1D_domainIndex(int i,char const * name)", ...
    "MATLABName", "clib.ctMatlabInterface.sim1D_domainIndex", ...
    "Description", "clib.ctMatlabInterface.sim1D_domainIndex Representation of C++ function sim1D_domainIndex."); % Modify help description values as needed.
defineArgument(sim1D_domainIndexDefinition, "i", "int32");
defineArgument(sim1D_domainIndexDefinition, "name", "string", "input", "nullTerminated");
defineOutput(sim1D_domainIndexDefinition, "RetVal", "int32");
validate(sim1D_domainIndexDefinition);

%% C++ function |sim1D_value| with MATLAB name |clib.ctMatlabInterface.sim1D_value|
% C++ Signature: double sim1D_value(int i,int idom,int icomp,int localPoint)

sim1D_valueDefinition = addFunction(libDef, ...
    "double sim1D_value(int i,int idom,int icomp,int localPoint)", ...
    "MATLABName", "clib.ctMatlabInterface.sim1D_value", ...
    "Description", "clib.ctMatlabInterface.sim1D_value Representation of C++ function sim1D_value."); % Modify help description values as needed.
defineArgument(sim1D_valueDefinition, "i", "int32");
defineArgument(sim1D_valueDefinition, "idom", "int32");
defineArgument(sim1D_valueDefinition, "icomp", "int32");
defineArgument(sim1D_valueDefinition, "localPoint", "int32");
defineOutput(sim1D_valueDefinition, "RetVal", "double");
validate(sim1D_valueDefinition);

%% C++ function |sim1D_workValue| with MATLAB name |clib.ctMatlabInterface.sim1D_workValue|
% C++ Signature: double sim1D_workValue(int i,int idom,int icomp,int localPoint)

sim1D_workValueDefinition = addFunction(libDef, ...
    "double sim1D_workValue(int i,int idom,int icomp,int localPoint)", ...
    "MATLABName", "clib.ctMatlabInterface.sim1D_workValue", ...
    "Description", "clib.ctMatlabInterface.sim1D_workValue Representation of C++ function sim1D_workValue."); % Modify help description values as needed.
defineArgument(sim1D_workValueDefinition, "i", "int32");
defineArgument(sim1D_workValueDefinition, "idom", "int32");
defineArgument(sim1D_workValueDefinition, "icomp", "int32");
defineArgument(sim1D_workValueDefinition, "localPoint", "int32");
defineOutput(sim1D_workValueDefinition, "RetVal", "double");
validate(sim1D_workValueDefinition);

%% C++ function |sim1D_eval| with MATLAB name |clib.ctMatlabInterface.sim1D_eval|
% C++ Signature: int sim1D_eval(int i,double rdt,int count)

sim1D_evalDefinition = addFunction(libDef, ...
    "int sim1D_eval(int i,double rdt,int count)", ...
    "MATLABName", "clib.ctMatlabInterface.sim1D_eval", ...
    "Description", "clib.ctMatlabInterface.sim1D_eval Representation of C++ function sim1D_eval."); % Modify help description values as needed.
defineArgument(sim1D_evalDefinition, "i", "int32");
defineArgument(sim1D_evalDefinition, "rdt", "double");
defineArgument(sim1D_evalDefinition, "count", "int32");
defineOutput(sim1D_evalDefinition, "RetVal", "int32");
validate(sim1D_evalDefinition);

%% C++ function |sim1D_setMaxJacAge| with MATLAB name |clib.ctMatlabInterface.sim1D_setMaxJacAge|
% C++ Signature: int sim1D_setMaxJacAge(int i,int ss_age,int ts_age)

sim1D_setMaxJacAgeDefinition = addFunction(libDef, ...
    "int sim1D_setMaxJacAge(int i,int ss_age,int ts_age)", ...
    "MATLABName", "clib.ctMatlabInterface.sim1D_setMaxJacAge", ...
    "Description", "clib.ctMatlabInterface.sim1D_setMaxJacAge Representation of C++ function sim1D_setMaxJacAge."); % Modify help description values as needed.
defineArgument(sim1D_setMaxJacAgeDefinition, "i", "int32");
defineArgument(sim1D_setMaxJacAgeDefinition, "ss_age", "int32");
defineArgument(sim1D_setMaxJacAgeDefinition, "ts_age", "int32");
defineOutput(sim1D_setMaxJacAgeDefinition, "RetVal", "int32");
validate(sim1D_setMaxJacAgeDefinition);

%% C++ function |sim1D_setFixedTemperature| with MATLAB name |clib.ctMatlabInterface.sim1D_setFixedTemperature|
% C++ Signature: int sim1D_setFixedTemperature(int i,double temp)

sim1D_setFixedTemperatureDefinition = addFunction(libDef, ...
    "int sim1D_setFixedTemperature(int i,double temp)", ...
    "MATLABName", "clib.ctMatlabInterface.sim1D_setFixedTemperature", ...
    "Description", "clib.ctMatlabInterface.sim1D_setFixedTemperature Representation of C++ function sim1D_setFixedTemperature."); % Modify help description values as needed.
defineArgument(sim1D_setFixedTemperatureDefinition, "i", "int32");
defineArgument(sim1D_setFixedTemperatureDefinition, "temp", "double");
defineOutput(sim1D_setFixedTemperatureDefinition, "RetVal", "int32");
validate(sim1D_setFixedTemperatureDefinition);

%% C++ function |reactor_new| with MATLAB name |clib.ctMatlabInterface.reactor_new|
% C++ Signature: int reactor_new(char const * type,int n,char const * name)

reactor_newDefinition = addFunction(libDef, ...
    "int reactor_new(char const * type,int n,char const * name)", ...
    "MATLABName", "clib.ctMatlabInterface.reactor_new", ...
    "Description", "clib.ctMatlabInterface.reactor_new Representation of C++ function reactor_new."); % Modify help description values as needed.
defineArgument(reactor_newDefinition, "type", "string", "input", "nullTerminated");
defineArgument(reactor_newDefinition, "n", "int32");
defineArgument(reactor_newDefinition, "name", "string", "input", "nullTerminated");
defineOutput(reactor_newDefinition, "RetVal", "int32");
validate(reactor_newDefinition);

%% C++ function |reactor_del| with MATLAB name |clib.ctMatlabInterface.reactor_del|
% C++ Signature: int reactor_del(int i)

reactor_delDefinition = addFunction(libDef, ...
    "int reactor_del(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.reactor_del", ...
    "Description", "clib.ctMatlabInterface.reactor_del Representation of C++ function reactor_del."); % Modify help description values as needed.
defineArgument(reactor_delDefinition, "i", "int32");
defineOutput(reactor_delDefinition, "RetVal", "int32");
validate(reactor_delDefinition);

%% C++ function |reactor_type| with MATLAB name |clib.ctMatlabInterface.reactor_type|
% C++ Signature: int reactor_type(int i,int len,char * nbuf)

reactor_typeDefinition = addFunction(libDef, ...
    "int reactor_type(int i,int len,char * nbuf)", ...
    "MATLABName", "clib.ctMatlabInterface.reactor_type", ...
    "Description", "clib.ctMatlabInterface.reactor_type Representation of C++ function reactor_type."); % Modify help description values as needed.
defineArgument(reactor_typeDefinition, "i", "int32");
defineArgument(reactor_typeDefinition, "len", "int32");
defineArgument(reactor_typeDefinition, "nbuf", "clib.array.ctMatlabInterface.Char", "input", "len"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(reactor_typeDefinition, "RetVal", "int32");
validate(reactor_typeDefinition);

%% C++ function |reactor_name| with MATLAB name |clib.ctMatlabInterface.reactor_name|
% C++ Signature: int reactor_name(int i,int len,char * nbuf)

reactor_nameDefinition = addFunction(libDef, ...
    "int reactor_name(int i,int len,char * nbuf)", ...
    "MATLABName", "clib.ctMatlabInterface.reactor_name", ...
    "Description", "clib.ctMatlabInterface.reactor_name Representation of C++ function reactor_name."); % Modify help description values as needed.
defineArgument(reactor_nameDefinition, "i", "int32");
defineArgument(reactor_nameDefinition, "len", "int32");
defineArgument(reactor_nameDefinition, "nbuf", "clib.array.ctMatlabInterface.Char", "input", "len"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(reactor_nameDefinition, "RetVal", "int32");
validate(reactor_nameDefinition);

%% C++ function |reactor_setName| with MATLAB name |clib.ctMatlabInterface.reactor_setName|
% C++ Signature: int reactor_setName(int i,char const * name)

reactor_setNameDefinition = addFunction(libDef, ...
    "int reactor_setName(int i,char const * name)", ...
    "MATLABName", "clib.ctMatlabInterface.reactor_setName", ...
    "Description", "clib.ctMatlabInterface.reactor_setName Representation of C++ function reactor_setName."); % Modify help description values as needed.
defineArgument(reactor_setNameDefinition, "i", "int32");
defineArgument(reactor_setNameDefinition, "name", "string", "input", "nullTerminated");
defineOutput(reactor_setNameDefinition, "RetVal", "int32");
validate(reactor_setNameDefinition);

%% C++ function |reactor_setInitialVolume| with MATLAB name |clib.ctMatlabInterface.reactor_setInitialVolume|
% C++ Signature: int reactor_setInitialVolume(int i,double v)

reactor_setInitialVolumeDefinition = addFunction(libDef, ...
    "int reactor_setInitialVolume(int i,double v)", ...
    "MATLABName", "clib.ctMatlabInterface.reactor_setInitialVolume", ...
    "Description", "clib.ctMatlabInterface.reactor_setInitialVolume Representation of C++ function reactor_setInitialVolume."); % Modify help description values as needed.
defineArgument(reactor_setInitialVolumeDefinition, "i", "int32");
defineArgument(reactor_setInitialVolumeDefinition, "v", "double");
defineOutput(reactor_setInitialVolumeDefinition, "RetVal", "int32");
validate(reactor_setInitialVolumeDefinition);

%% C++ function |reactor_setChemistry| with MATLAB name |clib.ctMatlabInterface.reactor_setChemistry|
% C++ Signature: int reactor_setChemistry(int i,int cflag)

reactor_setChemistryDefinition = addFunction(libDef, ...
    "int reactor_setChemistry(int i,int cflag)", ...
    "MATLABName", "clib.ctMatlabInterface.reactor_setChemistry", ...
    "Description", "clib.ctMatlabInterface.reactor_setChemistry Representation of C++ function reactor_setChemistry."); % Modify help description values as needed.
defineArgument(reactor_setChemistryDefinition, "i", "int32");
defineArgument(reactor_setChemistryDefinition, "cflag", "int32");
defineOutput(reactor_setChemistryDefinition, "RetVal", "int32");
validate(reactor_setChemistryDefinition);

%% C++ function |reactor_setEnergy| with MATLAB name |clib.ctMatlabInterface.reactor_setEnergy|
% C++ Signature: int reactor_setEnergy(int i,int eflag)

reactor_setEnergyDefinition = addFunction(libDef, ...
    "int reactor_setEnergy(int i,int eflag)", ...
    "MATLABName", "clib.ctMatlabInterface.reactor_setEnergy", ...
    "Description", "clib.ctMatlabInterface.reactor_setEnergy Representation of C++ function reactor_setEnergy."); % Modify help description values as needed.
defineArgument(reactor_setEnergyDefinition, "i", "int32");
defineArgument(reactor_setEnergyDefinition, "eflag", "int32");
defineOutput(reactor_setEnergyDefinition, "RetVal", "int32");
validate(reactor_setEnergyDefinition);

%% C++ function |reactor_mass| with MATLAB name |clib.ctMatlabInterface.reactor_mass|
% C++ Signature: double reactor_mass(int i)

reactor_massDefinition = addFunction(libDef, ...
    "double reactor_mass(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.reactor_mass", ...
    "Description", "clib.ctMatlabInterface.reactor_mass Representation of C++ function reactor_mass."); % Modify help description values as needed.
defineArgument(reactor_massDefinition, "i", "int32");
defineOutput(reactor_massDefinition, "RetVal", "double");
validate(reactor_massDefinition);

%% C++ function |reactor_volume| with MATLAB name |clib.ctMatlabInterface.reactor_volume|
% C++ Signature: double reactor_volume(int i)

reactor_volumeDefinition = addFunction(libDef, ...
    "double reactor_volume(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.reactor_volume", ...
    "Description", "clib.ctMatlabInterface.reactor_volume Representation of C++ function reactor_volume."); % Modify help description values as needed.
defineArgument(reactor_volumeDefinition, "i", "int32");
defineOutput(reactor_volumeDefinition, "RetVal", "double");
validate(reactor_volumeDefinition);

%% C++ function |reactor_density| with MATLAB name |clib.ctMatlabInterface.reactor_density|
% C++ Signature: double reactor_density(int i)

reactor_densityDefinition = addFunction(libDef, ...
    "double reactor_density(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.reactor_density", ...
    "Description", "clib.ctMatlabInterface.reactor_density Representation of C++ function reactor_density."); % Modify help description values as needed.
defineArgument(reactor_densityDefinition, "i", "int32");
defineOutput(reactor_densityDefinition, "RetVal", "double");
validate(reactor_densityDefinition);

%% C++ function |reactor_temperature| with MATLAB name |clib.ctMatlabInterface.reactor_temperature|
% C++ Signature: double reactor_temperature(int i)

reactor_temperatureDefinition = addFunction(libDef, ...
    "double reactor_temperature(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.reactor_temperature", ...
    "Description", "clib.ctMatlabInterface.reactor_temperature Representation of C++ function reactor_temperature."); % Modify help description values as needed.
defineArgument(reactor_temperatureDefinition, "i", "int32");
defineOutput(reactor_temperatureDefinition, "RetVal", "double");
validate(reactor_temperatureDefinition);

%% C++ function |reactor_enthalpy_mass| with MATLAB name |clib.ctMatlabInterface.reactor_enthalpy_mass|
% C++ Signature: double reactor_enthalpy_mass(int i)

reactor_enthalpy_massDefinition = addFunction(libDef, ...
    "double reactor_enthalpy_mass(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.reactor_enthalpy_mass", ...
    "Description", "clib.ctMatlabInterface.reactor_enthalpy_mass Representation of C++ function reactor_enthalpy_mass."); % Modify help description values as needed.
defineArgument(reactor_enthalpy_massDefinition, "i", "int32");
defineOutput(reactor_enthalpy_massDefinition, "RetVal", "double");
validate(reactor_enthalpy_massDefinition);

%% C++ function |reactor_intEnergy_mass| with MATLAB name |clib.ctMatlabInterface.reactor_intEnergy_mass|
% C++ Signature: double reactor_intEnergy_mass(int i)

reactor_intEnergy_massDefinition = addFunction(libDef, ...
    "double reactor_intEnergy_mass(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.reactor_intEnergy_mass", ...
    "Description", "clib.ctMatlabInterface.reactor_intEnergy_mass Representation of C++ function reactor_intEnergy_mass."); % Modify help description values as needed.
defineArgument(reactor_intEnergy_massDefinition, "i", "int32");
defineOutput(reactor_intEnergy_massDefinition, "RetVal", "double");
validate(reactor_intEnergy_massDefinition);

%% C++ function |reactor_pressure| with MATLAB name |clib.ctMatlabInterface.reactor_pressure|
% C++ Signature: double reactor_pressure(int i)

reactor_pressureDefinition = addFunction(libDef, ...
    "double reactor_pressure(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.reactor_pressure", ...
    "Description", "clib.ctMatlabInterface.reactor_pressure Representation of C++ function reactor_pressure."); % Modify help description values as needed.
defineArgument(reactor_pressureDefinition, "i", "int32");
defineOutput(reactor_pressureDefinition, "RetVal", "double");
validate(reactor_pressureDefinition);

%% C++ function |reactor_massFraction| with MATLAB name |clib.ctMatlabInterface.reactor_massFraction|
% C++ Signature: double reactor_massFraction(int i,int k)

reactor_massFractionDefinition = addFunction(libDef, ...
    "double reactor_massFraction(int i,int k)", ...
    "MATLABName", "clib.ctMatlabInterface.reactor_massFraction", ...
    "Description", "clib.ctMatlabInterface.reactor_massFraction Representation of C++ function reactor_massFraction."); % Modify help description values as needed.
defineArgument(reactor_massFractionDefinition, "i", "int32");
defineArgument(reactor_massFractionDefinition, "k", "int32");
defineOutput(reactor_massFractionDefinition, "RetVal", "double");
validate(reactor_massFractionDefinition);

%% C++ function |reactor_nSensParams| with MATLAB name |clib.ctMatlabInterface.reactor_nSensParams|
% C++ Signature: size_t reactor_nSensParams(int i)

reactor_nSensParamsDefinition = addFunction(libDef, ...
    "size_t reactor_nSensParams(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.reactor_nSensParams", ...
    "Description", "clib.ctMatlabInterface.reactor_nSensParams Representation of C++ function reactor_nSensParams."); % Modify help description values as needed.
defineArgument(reactor_nSensParamsDefinition, "i", "int32");
defineOutput(reactor_nSensParamsDefinition, "RetVal", "uint64");
validate(reactor_nSensParamsDefinition);

%% C++ function |reactor_addSensitivityReaction| with MATLAB name |clib.ctMatlabInterface.reactor_addSensitivityReaction|
% C++ Signature: int reactor_addSensitivityReaction(int i,int rxn)

reactor_addSensitivityReactionDefinition = addFunction(libDef, ...
    "int reactor_addSensitivityReaction(int i,int rxn)", ...
    "MATLABName", "clib.ctMatlabInterface.reactor_addSensitivityReaction", ...
    "Description", "clib.ctMatlabInterface.reactor_addSensitivityReaction Representation of C++ function reactor_addSensitivityReaction."); % Modify help description values as needed.
defineArgument(reactor_addSensitivityReactionDefinition, "i", "int32");
defineArgument(reactor_addSensitivityReactionDefinition, "rxn", "int32");
defineOutput(reactor_addSensitivityReactionDefinition, "RetVal", "int32");
validate(reactor_addSensitivityReactionDefinition);

%% C++ function |reactor_setMassFlowRate| with MATLAB name |clib.ctMatlabInterface.reactor_setMassFlowRate|
% C++ Signature: int reactor_setMassFlowRate(int i,double mdot)

reactor_setMassFlowRateDefinition = addFunction(libDef, ...
    "int reactor_setMassFlowRate(int i,double mdot)", ...
    "MATLABName", "clib.ctMatlabInterface.reactor_setMassFlowRate", ...
    "Description", "clib.ctMatlabInterface.reactor_setMassFlowRate Representation of C++ function reactor_setMassFlowRate."); % Modify help description values as needed.
defineArgument(reactor_setMassFlowRateDefinition, "i", "int32");
defineArgument(reactor_setMassFlowRateDefinition, "mdot", "double");
defineOutput(reactor_setMassFlowRateDefinition, "RetVal", "int32");
validate(reactor_setMassFlowRateDefinition);

%% C++ function |reactornet_new| with MATLAB name |clib.ctMatlabInterface.reactornet_new|
% C++ Signature: int reactornet_new()

reactornet_newDefinition = addFunction(libDef, ...
    "int reactornet_new()", ...
    "MATLABName", "clib.ctMatlabInterface.reactornet_new", ...
    "Description", "clib.ctMatlabInterface.reactornet_new Representation of C++ function reactornet_new."); % Modify help description values as needed.
defineOutput(reactornet_newDefinition, "RetVal", "int32");
validate(reactornet_newDefinition);

%% C++ function |reactornet_del| with MATLAB name |clib.ctMatlabInterface.reactornet_del|
% C++ Signature: int reactornet_del(int i)

reactornet_delDefinition = addFunction(libDef, ...
    "int reactornet_del(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.reactornet_del", ...
    "Description", "clib.ctMatlabInterface.reactornet_del Representation of C++ function reactornet_del."); % Modify help description values as needed.
defineArgument(reactornet_delDefinition, "i", "int32");
defineOutput(reactornet_delDefinition, "RetVal", "int32");
validate(reactornet_delDefinition);

%% C++ function |reactornet_setInitialTime| with MATLAB name |clib.ctMatlabInterface.reactornet_setInitialTime|
% C++ Signature: int reactornet_setInitialTime(int i,double t)

reactornet_setInitialTimeDefinition = addFunction(libDef, ...
    "int reactornet_setInitialTime(int i,double t)", ...
    "MATLABName", "clib.ctMatlabInterface.reactornet_setInitialTime", ...
    "Description", "clib.ctMatlabInterface.reactornet_setInitialTime Representation of C++ function reactornet_setInitialTime."); % Modify help description values as needed.
defineArgument(reactornet_setInitialTimeDefinition, "i", "int32");
defineArgument(reactornet_setInitialTimeDefinition, "t", "double");
defineOutput(reactornet_setInitialTimeDefinition, "RetVal", "int32");
validate(reactornet_setInitialTimeDefinition);

%% C++ function |reactornet_setMaxTimeStep| with MATLAB name |clib.ctMatlabInterface.reactornet_setMaxTimeStep|
% C++ Signature: int reactornet_setMaxTimeStep(int i,double maxstep)

reactornet_setMaxTimeStepDefinition = addFunction(libDef, ...
    "int reactornet_setMaxTimeStep(int i,double maxstep)", ...
    "MATLABName", "clib.ctMatlabInterface.reactornet_setMaxTimeStep", ...
    "Description", "clib.ctMatlabInterface.reactornet_setMaxTimeStep Representation of C++ function reactornet_setMaxTimeStep."); % Modify help description values as needed.
defineArgument(reactornet_setMaxTimeStepDefinition, "i", "int32");
defineArgument(reactornet_setMaxTimeStepDefinition, "maxstep", "double");
defineOutput(reactornet_setMaxTimeStepDefinition, "RetVal", "int32");
validate(reactornet_setMaxTimeStepDefinition);

%% C++ function |reactornet_setTolerances| with MATLAB name |clib.ctMatlabInterface.reactornet_setTolerances|
% C++ Signature: int reactornet_setTolerances(int i,double rtol,double atol)

reactornet_setTolerancesDefinition = addFunction(libDef, ...
    "int reactornet_setTolerances(int i,double rtol,double atol)", ...
    "MATLABName", "clib.ctMatlabInterface.reactornet_setTolerances", ...
    "Description", "clib.ctMatlabInterface.reactornet_setTolerances Representation of C++ function reactornet_setTolerances."); % Modify help description values as needed.
defineArgument(reactornet_setTolerancesDefinition, "i", "int32");
defineArgument(reactornet_setTolerancesDefinition, "rtol", "double");
defineArgument(reactornet_setTolerancesDefinition, "atol", "double");
defineOutput(reactornet_setTolerancesDefinition, "RetVal", "int32");
validate(reactornet_setTolerancesDefinition);

%% C++ function |reactornet_setSensitivityTolerances| with MATLAB name |clib.ctMatlabInterface.reactornet_setSensitivityTolerances|
% C++ Signature: int reactornet_setSensitivityTolerances(int i,double rtol,double atol)

reactornet_setSensitivityTolerancesDefinition = addFunction(libDef, ...
    "int reactornet_setSensitivityTolerances(int i,double rtol,double atol)", ...
    "MATLABName", "clib.ctMatlabInterface.reactornet_setSensitivityTolerances", ...
    "Description", "clib.ctMatlabInterface.reactornet_setSensitivityTolerances Representation of C++ function reactornet_setSensitivityTolerances."); % Modify help description values as needed.
defineArgument(reactornet_setSensitivityTolerancesDefinition, "i", "int32");
defineArgument(reactornet_setSensitivityTolerancesDefinition, "rtol", "double");
defineArgument(reactornet_setSensitivityTolerancesDefinition, "atol", "double");
defineOutput(reactornet_setSensitivityTolerancesDefinition, "RetVal", "int32");
validate(reactornet_setSensitivityTolerancesDefinition);

%% C++ function |reactornet_addreactor| with MATLAB name |clib.ctMatlabInterface.reactornet_addreactor|
% C++ Signature: int reactornet_addreactor(int i,int n)

reactornet_addreactorDefinition = addFunction(libDef, ...
    "int reactornet_addreactor(int i,int n)", ...
    "MATLABName", "clib.ctMatlabInterface.reactornet_addreactor", ...
    "Description", "clib.ctMatlabInterface.reactornet_addreactor Representation of C++ function reactornet_addreactor."); % Modify help description values as needed.
defineArgument(reactornet_addreactorDefinition, "i", "int32");
defineArgument(reactornet_addreactorDefinition, "n", "int32");
defineOutput(reactornet_addreactorDefinition, "RetVal", "int32");
validate(reactornet_addreactorDefinition);

%% C++ function |reactornet_advance| with MATLAB name |clib.ctMatlabInterface.reactornet_advance|
% C++ Signature: int reactornet_advance(int i,double t)

reactornet_advanceDefinition = addFunction(libDef, ...
    "int reactornet_advance(int i,double t)", ...
    "MATLABName", "clib.ctMatlabInterface.reactornet_advance", ...
    "Description", "clib.ctMatlabInterface.reactornet_advance Representation of C++ function reactornet_advance."); % Modify help description values as needed.
defineArgument(reactornet_advanceDefinition, "i", "int32");
defineArgument(reactornet_advanceDefinition, "t", "double");
defineOutput(reactornet_advanceDefinition, "RetVal", "int32");
validate(reactornet_advanceDefinition);

%% C++ function |reactornet_step| with MATLAB name |clib.ctMatlabInterface.reactornet_step|
% C++ Signature: double reactornet_step(int i)

reactornet_stepDefinition = addFunction(libDef, ...
    "double reactornet_step(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.reactornet_step", ...
    "Description", "clib.ctMatlabInterface.reactornet_step Representation of C++ function reactornet_step."); % Modify help description values as needed.
defineArgument(reactornet_stepDefinition, "i", "int32");
defineOutput(reactornet_stepDefinition, "RetVal", "double");
validate(reactornet_stepDefinition);

%% C++ function |reactornet_time| with MATLAB name |clib.ctMatlabInterface.reactornet_time|
% C++ Signature: double reactornet_time(int i)

reactornet_timeDefinition = addFunction(libDef, ...
    "double reactornet_time(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.reactornet_time", ...
    "Description", "clib.ctMatlabInterface.reactornet_time Representation of C++ function reactornet_time."); % Modify help description values as needed.
defineArgument(reactornet_timeDefinition, "i", "int32");
defineOutput(reactornet_timeDefinition, "RetVal", "double");
validate(reactornet_timeDefinition);

%% C++ function |reactornet_rtol| with MATLAB name |clib.ctMatlabInterface.reactornet_rtol|
% C++ Signature: double reactornet_rtol(int i)

reactornet_rtolDefinition = addFunction(libDef, ...
    "double reactornet_rtol(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.reactornet_rtol", ...
    "Description", "clib.ctMatlabInterface.reactornet_rtol Representation of C++ function reactornet_rtol."); % Modify help description values as needed.
defineArgument(reactornet_rtolDefinition, "i", "int32");
defineOutput(reactornet_rtolDefinition, "RetVal", "double");
validate(reactornet_rtolDefinition);

%% C++ function |reactornet_atol| with MATLAB name |clib.ctMatlabInterface.reactornet_atol|
% C++ Signature: double reactornet_atol(int i)

reactornet_atolDefinition = addFunction(libDef, ...
    "double reactornet_atol(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.reactornet_atol", ...
    "Description", "clib.ctMatlabInterface.reactornet_atol Representation of C++ function reactornet_atol."); % Modify help description values as needed.
defineArgument(reactornet_atolDefinition, "i", "int32");
defineOutput(reactornet_atolDefinition, "RetVal", "double");
validate(reactornet_atolDefinition);

%% C++ function |reactornet_sensitivity| with MATLAB name |clib.ctMatlabInterface.reactornet_sensitivity|
% C++ Signature: double reactornet_sensitivity(int i,char const * v,int p,int r)

reactornet_sensitivityDefinition = addFunction(libDef, ...
    "double reactornet_sensitivity(int i,char const * v,int p,int r)", ...
    "MATLABName", "clib.ctMatlabInterface.reactornet_sensitivity", ...
    "Description", "clib.ctMatlabInterface.reactornet_sensitivity Representation of C++ function reactornet_sensitivity."); % Modify help description values as needed.
defineArgument(reactornet_sensitivityDefinition, "i", "int32");
defineArgument(reactornet_sensitivityDefinition, "v", "string", "input", "nullTerminated");
defineArgument(reactornet_sensitivityDefinition, "p", "int32");
defineArgument(reactornet_sensitivityDefinition, "r", "int32");
defineOutput(reactornet_sensitivityDefinition, "RetVal", "double");
validate(reactornet_sensitivityDefinition);

%% C++ function |connector_new| with MATLAB name |clib.ctMatlabInterface.connector_new|
% C++ Signature: int connector_new(char const * type,int r0,int r1,char const * name)

connector_newDefinition = addFunction(libDef, ...
    "int connector_new(char const * type,int r0,int r1,char const * name)", ...
    "MATLABName", "clib.ctMatlabInterface.connector_new", ...
    "Description", "clib.ctMatlabInterface.connector_new Representation of C++ function connector_new."); % Modify help description values as needed.
defineArgument(connector_newDefinition, "type", "string", "input", "nullTerminated");
defineArgument(connector_newDefinition, "r0", "int32");
defineArgument(connector_newDefinition, "r1", "int32");
defineArgument(connector_newDefinition, "name", "string", "input", "nullTerminated");
defineOutput(connector_newDefinition, "RetVal", "int32");
validate(connector_newDefinition);

%% C++ function |connector_del| with MATLAB name |clib.ctMatlabInterface.connector_del|
% C++ Signature: int connector_del(int i)

connector_delDefinition = addFunction(libDef, ...
    "int connector_del(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.connector_del", ...
    "Description", "clib.ctMatlabInterface.connector_del Representation of C++ function connector_del."); % Modify help description values as needed.
defineArgument(connector_delDefinition, "i", "int32");
defineOutput(connector_delDefinition, "RetVal", "int32");
validate(connector_delDefinition);

%% C++ function |connector_type| with MATLAB name |clib.ctMatlabInterface.connector_type|
% C++ Signature: int connector_type(int i,int len,char * nbuf)

connector_typeDefinition = addFunction(libDef, ...
    "int connector_type(int i,int len,char * nbuf)", ...
    "MATLABName", "clib.ctMatlabInterface.connector_type", ...
    "Description", "clib.ctMatlabInterface.connector_type Representation of C++ function connector_type."); % Modify help description values as needed.
defineArgument(connector_typeDefinition, "i", "int32");
defineArgument(connector_typeDefinition, "len", "int32");
defineArgument(connector_typeDefinition, "nbuf", "clib.array.ctMatlabInterface.Char", "input", "len"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(connector_typeDefinition, "RetVal", "int32");
validate(connector_typeDefinition);

%% C++ function |connector_name| with MATLAB name |clib.ctMatlabInterface.connector_name|
% C++ Signature: int connector_name(int i,int len,char * nbuf)

connector_nameDefinition = addFunction(libDef, ...
    "int connector_name(int i,int len,char * nbuf)", ...
    "MATLABName", "clib.ctMatlabInterface.connector_name", ...
    "Description", "clib.ctMatlabInterface.connector_name Representation of C++ function connector_name."); % Modify help description values as needed.
defineArgument(connector_nameDefinition, "i", "int32");
defineArgument(connector_nameDefinition, "len", "int32");
defineArgument(connector_nameDefinition, "nbuf", "clib.array.ctMatlabInterface.Char", "input", "len"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(connector_nameDefinition, "RetVal", "int32");
validate(connector_nameDefinition);

%% C++ function |connector_setName| with MATLAB name |clib.ctMatlabInterface.connector_setName|
% C++ Signature: int connector_setName(int i,char const * name)

connector_setNameDefinition = addFunction(libDef, ...
    "int connector_setName(int i,char const * name)", ...
    "MATLABName", "clib.ctMatlabInterface.connector_setName", ...
    "Description", "clib.ctMatlabInterface.connector_setName Representation of C++ function connector_setName."); % Modify help description values as needed.
defineArgument(connector_setNameDefinition, "i", "int32");
defineArgument(connector_setNameDefinition, "name", "string", "input", "nullTerminated");
defineOutput(connector_setNameDefinition, "RetVal", "int32");
validate(connector_setNameDefinition);

%% C++ function |flowdev_setPrimary| with MATLAB name |clib.ctMatlabInterface.flowdev_setPrimary|
% C++ Signature: int flowdev_setPrimary(int i,int n)

flowdev_setPrimaryDefinition = addFunction(libDef, ...
    "int flowdev_setPrimary(int i,int n)", ...
    "MATLABName", "clib.ctMatlabInterface.flowdev_setPrimary", ...
    "Description", "clib.ctMatlabInterface.flowdev_setPrimary Representation of C++ function flowdev_setPrimary."); % Modify help description values as needed.
defineArgument(flowdev_setPrimaryDefinition, "i", "int32");
defineArgument(flowdev_setPrimaryDefinition, "n", "int32");
defineOutput(flowdev_setPrimaryDefinition, "RetVal", "int32");
validate(flowdev_setPrimaryDefinition);

%% C++ function |flowdev_massFlowRate| with MATLAB name |clib.ctMatlabInterface.flowdev_massFlowRate|
% C++ Signature: double flowdev_massFlowRate(int i)

flowdev_massFlowRateDefinition = addFunction(libDef, ...
    "double flowdev_massFlowRate(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.flowdev_massFlowRate", ...
    "Description", "clib.ctMatlabInterface.flowdev_massFlowRate Representation of C++ function flowdev_massFlowRate."); % Modify help description values as needed.
defineArgument(flowdev_massFlowRateDefinition, "i", "int32");
defineOutput(flowdev_massFlowRateDefinition, "RetVal", "double");
validate(flowdev_massFlowRateDefinition);

%% C++ function |flowdev_setMassFlowCoeff| with MATLAB name |clib.ctMatlabInterface.flowdev_setMassFlowCoeff|
% C++ Signature: int flowdev_setMassFlowCoeff(int i,double v)

flowdev_setMassFlowCoeffDefinition = addFunction(libDef, ...
    "int flowdev_setMassFlowCoeff(int i,double v)", ...
    "MATLABName", "clib.ctMatlabInterface.flowdev_setMassFlowCoeff", ...
    "Description", "clib.ctMatlabInterface.flowdev_setMassFlowCoeff Representation of C++ function flowdev_setMassFlowCoeff."); % Modify help description values as needed.
defineArgument(flowdev_setMassFlowCoeffDefinition, "i", "int32");
defineArgument(flowdev_setMassFlowCoeffDefinition, "v", "double");
defineOutput(flowdev_setMassFlowCoeffDefinition, "RetVal", "int32");
validate(flowdev_setMassFlowCoeffDefinition);

%% C++ function |flowdev_setValveCoeff| with MATLAB name |clib.ctMatlabInterface.flowdev_setValveCoeff|
% C++ Signature: int flowdev_setValveCoeff(int i,double v)

flowdev_setValveCoeffDefinition = addFunction(libDef, ...
    "int flowdev_setValveCoeff(int i,double v)", ...
    "MATLABName", "clib.ctMatlabInterface.flowdev_setValveCoeff", ...
    "Description", "clib.ctMatlabInterface.flowdev_setValveCoeff Representation of C++ function flowdev_setValveCoeff."); % Modify help description values as needed.
defineArgument(flowdev_setValveCoeffDefinition, "i", "int32");
defineArgument(flowdev_setValveCoeffDefinition, "v", "double");
defineOutput(flowdev_setValveCoeffDefinition, "RetVal", "int32");
validate(flowdev_setValveCoeffDefinition);

%% C++ function |flowdev_setPressureCoeff| with MATLAB name |clib.ctMatlabInterface.flowdev_setPressureCoeff|
% C++ Signature: int flowdev_setPressureCoeff(int i,double v)

flowdev_setPressureCoeffDefinition = addFunction(libDef, ...
    "int flowdev_setPressureCoeff(int i,double v)", ...
    "MATLABName", "clib.ctMatlabInterface.flowdev_setPressureCoeff", ...
    "Description", "clib.ctMatlabInterface.flowdev_setPressureCoeff Representation of C++ function flowdev_setPressureCoeff."); % Modify help description values as needed.
defineArgument(flowdev_setPressureCoeffDefinition, "i", "int32");
defineArgument(flowdev_setPressureCoeffDefinition, "v", "double");
defineOutput(flowdev_setPressureCoeffDefinition, "RetVal", "int32");
validate(flowdev_setPressureCoeffDefinition);

%% C++ function |flowdev_setPressureFunction| with MATLAB name |clib.ctMatlabInterface.flowdev_setPressureFunction|
% C++ Signature: int flowdev_setPressureFunction(int i,int n)

flowdev_setPressureFunctionDefinition = addFunction(libDef, ...
    "int flowdev_setPressureFunction(int i,int n)", ...
    "MATLABName", "clib.ctMatlabInterface.flowdev_setPressureFunction", ...
    "Description", "clib.ctMatlabInterface.flowdev_setPressureFunction Representation of C++ function flowdev_setPressureFunction."); % Modify help description values as needed.
defineArgument(flowdev_setPressureFunctionDefinition, "i", "int32");
defineArgument(flowdev_setPressureFunctionDefinition, "n", "int32");
defineOutput(flowdev_setPressureFunctionDefinition, "RetVal", "int32");
validate(flowdev_setPressureFunctionDefinition);

%% C++ function |flowdev_setTimeFunction| with MATLAB name |clib.ctMatlabInterface.flowdev_setTimeFunction|
% C++ Signature: int flowdev_setTimeFunction(int i,int n)

flowdev_setTimeFunctionDefinition = addFunction(libDef, ...
    "int flowdev_setTimeFunction(int i,int n)", ...
    "MATLABName", "clib.ctMatlabInterface.flowdev_setTimeFunction", ...
    "Description", "clib.ctMatlabInterface.flowdev_setTimeFunction Representation of C++ function flowdev_setTimeFunction."); % Modify help description values as needed.
defineArgument(flowdev_setTimeFunctionDefinition, "i", "int32");
defineArgument(flowdev_setTimeFunctionDefinition, "n", "int32");
defineOutput(flowdev_setTimeFunctionDefinition, "RetVal", "int32");
validate(flowdev_setTimeFunctionDefinition);

%% C++ function |wall_expansionRate| with MATLAB name |clib.ctMatlabInterface.wall_expansionRate|
% C++ Signature: double wall_expansionRate(int i)

wall_expansionRateDefinition = addFunction(libDef, ...
    "double wall_expansionRate(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.wall_expansionRate", ...
    "Description", "clib.ctMatlabInterface.wall_expansionRate Representation of C++ function wall_expansionRate."); % Modify help description values as needed.
defineArgument(wall_expansionRateDefinition, "i", "int32");
defineOutput(wall_expansionRateDefinition, "RetVal", "double");
validate(wall_expansionRateDefinition);

%% C++ function |wall_heatRate| with MATLAB name |clib.ctMatlabInterface.wall_heatRate|
% C++ Signature: double wall_heatRate(int i)

wall_heatRateDefinition = addFunction(libDef, ...
    "double wall_heatRate(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.wall_heatRate", ...
    "Description", "clib.ctMatlabInterface.wall_heatRate Representation of C++ function wall_heatRate."); % Modify help description values as needed.
defineArgument(wall_heatRateDefinition, "i", "int32");
defineOutput(wall_heatRateDefinition, "RetVal", "double");
validate(wall_heatRateDefinition);

%% C++ function |wall_area| with MATLAB name |clib.ctMatlabInterface.wall_area|
% C++ Signature: double wall_area(int i)

wall_areaDefinition = addFunction(libDef, ...
    "double wall_area(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.wall_area", ...
    "Description", "clib.ctMatlabInterface.wall_area Representation of C++ function wall_area."); % Modify help description values as needed.
defineArgument(wall_areaDefinition, "i", "int32");
defineOutput(wall_areaDefinition, "RetVal", "double");
validate(wall_areaDefinition);

%% C++ function |wall_setArea| with MATLAB name |clib.ctMatlabInterface.wall_setArea|
% C++ Signature: int wall_setArea(int i,double v)

wall_setAreaDefinition = addFunction(libDef, ...
    "int wall_setArea(int i,double v)", ...
    "MATLABName", "clib.ctMatlabInterface.wall_setArea", ...
    "Description", "clib.ctMatlabInterface.wall_setArea Representation of C++ function wall_setArea."); % Modify help description values as needed.
defineArgument(wall_setAreaDefinition, "i", "int32");
defineArgument(wall_setAreaDefinition, "v", "double");
defineOutput(wall_setAreaDefinition, "RetVal", "int32");
validate(wall_setAreaDefinition);

%% C++ function |wall_setThermalResistance| with MATLAB name |clib.ctMatlabInterface.wall_setThermalResistance|
% C++ Signature: int wall_setThermalResistance(int i,double rth)

wall_setThermalResistanceDefinition = addFunction(libDef, ...
    "int wall_setThermalResistance(int i,double rth)", ...
    "MATLABName", "clib.ctMatlabInterface.wall_setThermalResistance", ...
    "Description", "clib.ctMatlabInterface.wall_setThermalResistance Representation of C++ function wall_setThermalResistance."); % Modify help description values as needed.
defineArgument(wall_setThermalResistanceDefinition, "i", "int32");
defineArgument(wall_setThermalResistanceDefinition, "rth", "double");
defineOutput(wall_setThermalResistanceDefinition, "RetVal", "int32");
validate(wall_setThermalResistanceDefinition);

%% C++ function |wall_setHeatTransferCoeff| with MATLAB name |clib.ctMatlabInterface.wall_setHeatTransferCoeff|
% C++ Signature: int wall_setHeatTransferCoeff(int i,double u)

wall_setHeatTransferCoeffDefinition = addFunction(libDef, ...
    "int wall_setHeatTransferCoeff(int i,double u)", ...
    "MATLABName", "clib.ctMatlabInterface.wall_setHeatTransferCoeff", ...
    "Description", "clib.ctMatlabInterface.wall_setHeatTransferCoeff Representation of C++ function wall_setHeatTransferCoeff."); % Modify help description values as needed.
defineArgument(wall_setHeatTransferCoeffDefinition, "i", "int32");
defineArgument(wall_setHeatTransferCoeffDefinition, "u", "double");
defineOutput(wall_setHeatTransferCoeffDefinition, "RetVal", "int32");
validate(wall_setHeatTransferCoeffDefinition);

%% C++ function |wall_setHeatFlux| with MATLAB name |clib.ctMatlabInterface.wall_setHeatFlux|
% C++ Signature: int wall_setHeatFlux(int i,int n)

wall_setHeatFluxDefinition = addFunction(libDef, ...
    "int wall_setHeatFlux(int i,int n)", ...
    "MATLABName", "clib.ctMatlabInterface.wall_setHeatFlux", ...
    "Description", "clib.ctMatlabInterface.wall_setHeatFlux Representation of C++ function wall_setHeatFlux."); % Modify help description values as needed.
defineArgument(wall_setHeatFluxDefinition, "i", "int32");
defineArgument(wall_setHeatFluxDefinition, "n", "int32");
defineOutput(wall_setHeatFluxDefinition, "RetVal", "int32");
validate(wall_setHeatFluxDefinition);

%% C++ function |wall_setExpansionRateCoeff| with MATLAB name |clib.ctMatlabInterface.wall_setExpansionRateCoeff|
% C++ Signature: int wall_setExpansionRateCoeff(int i,double k)

wall_setExpansionRateCoeffDefinition = addFunction(libDef, ...
    "int wall_setExpansionRateCoeff(int i,double k)", ...
    "MATLABName", "clib.ctMatlabInterface.wall_setExpansionRateCoeff", ...
    "Description", "clib.ctMatlabInterface.wall_setExpansionRateCoeff Representation of C++ function wall_setExpansionRateCoeff."); % Modify help description values as needed.
defineArgument(wall_setExpansionRateCoeffDefinition, "i", "int32");
defineArgument(wall_setExpansionRateCoeffDefinition, "k", "double");
defineOutput(wall_setExpansionRateCoeffDefinition, "RetVal", "int32");
validate(wall_setExpansionRateCoeffDefinition);

%% C++ function |wall_setVelocity| with MATLAB name |clib.ctMatlabInterface.wall_setVelocity|
% C++ Signature: int wall_setVelocity(int i,int n)

wall_setVelocityDefinition = addFunction(libDef, ...
    "int wall_setVelocity(int i,int n)", ...
    "MATLABName", "clib.ctMatlabInterface.wall_setVelocity", ...
    "Description", "clib.ctMatlabInterface.wall_setVelocity Representation of C++ function wall_setVelocity."); % Modify help description values as needed.
defineArgument(wall_setVelocityDefinition, "i", "int32");
defineArgument(wall_setVelocityDefinition, "n", "int32");
defineOutput(wall_setVelocityDefinition, "RetVal", "int32");
validate(wall_setVelocityDefinition);

%% C++ function |wall_setEmissivity| with MATLAB name |clib.ctMatlabInterface.wall_setEmissivity|
% C++ Signature: int wall_setEmissivity(int i,double epsilon)

wall_setEmissivityDefinition = addFunction(libDef, ...
    "int wall_setEmissivity(int i,double epsilon)", ...
    "MATLABName", "clib.ctMatlabInterface.wall_setEmissivity", ...
    "Description", "clib.ctMatlabInterface.wall_setEmissivity Representation of C++ function wall_setEmissivity."); % Modify help description values as needed.
defineArgument(wall_setEmissivityDefinition, "i", "int32");
defineArgument(wall_setEmissivityDefinition, "epsilon", "double");
defineOutput(wall_setEmissivityDefinition, "RetVal", "int32");
validate(wall_setEmissivityDefinition);

%% C++ function |reactorsurface_new| with MATLAB name |clib.ctMatlabInterface.reactorsurface_new|
% C++ Signature: int reactorsurface_new(char const * name)

reactorsurface_newDefinition = addFunction(libDef, ...
    "int reactorsurface_new(char const * name)", ...
    "MATLABName", "clib.ctMatlabInterface.reactorsurface_new", ...
    "Description", "clib.ctMatlabInterface.reactorsurface_new Representation of C++ function reactorsurface_new."); % Modify help description values as needed.
defineArgument(reactorsurface_newDefinition, "name", "string", "input", "nullTerminated");
defineOutput(reactorsurface_newDefinition, "RetVal", "int32");
validate(reactorsurface_newDefinition);

%% C++ function |reactorsurface_del| with MATLAB name |clib.ctMatlabInterface.reactorsurface_del|
% C++ Signature: int reactorsurface_del(int i)

reactorsurface_delDefinition = addFunction(libDef, ...
    "int reactorsurface_del(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.reactorsurface_del", ...
    "Description", "clib.ctMatlabInterface.reactorsurface_del Representation of C++ function reactorsurface_del."); % Modify help description values as needed.
defineArgument(reactorsurface_delDefinition, "i", "int32");
defineOutput(reactorsurface_delDefinition, "RetVal", "int32");
validate(reactorsurface_delDefinition);

%% C++ function |reactorsurface_name| with MATLAB name |clib.ctMatlabInterface.reactorsurface_name|
% C++ Signature: int reactorsurface_name(int i,int len,char * nbuf)

reactorsurface_nameDefinition = addFunction(libDef, ...
    "int reactorsurface_name(int i,int len,char * nbuf)", ...
    "MATLABName", "clib.ctMatlabInterface.reactorsurface_name", ...
    "Description", "clib.ctMatlabInterface.reactorsurface_name Representation of C++ function reactorsurface_name."); % Modify help description values as needed.
defineArgument(reactorsurface_nameDefinition, "i", "int32");
defineArgument(reactorsurface_nameDefinition, "len", "int32");
defineArgument(reactorsurface_nameDefinition, "nbuf", "clib.array.ctMatlabInterface.Char", "input", "len"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Char","int8","string", or "char"
defineOutput(reactorsurface_nameDefinition, "RetVal", "int32");
validate(reactorsurface_nameDefinition);

%% C++ function |reactorsurface_setName| with MATLAB name |clib.ctMatlabInterface.reactorsurface_setName|
% C++ Signature: int reactorsurface_setName(int i,char const * name)

reactorsurface_setNameDefinition = addFunction(libDef, ...
    "int reactorsurface_setName(int i,char const * name)", ...
    "MATLABName", "clib.ctMatlabInterface.reactorsurface_setName", ...
    "Description", "clib.ctMatlabInterface.reactorsurface_setName Representation of C++ function reactorsurface_setName."); % Modify help description values as needed.
defineArgument(reactorsurface_setNameDefinition, "i", "int32");
defineArgument(reactorsurface_setNameDefinition, "name", "string", "input", "nullTerminated");
defineOutput(reactorsurface_setNameDefinition, "RetVal", "int32");
validate(reactorsurface_setNameDefinition);

%% C++ function |reactorsurface_install| with MATLAB name |clib.ctMatlabInterface.reactorsurface_install|
% C++ Signature: int reactorsurface_install(int i,int n)

reactorsurface_installDefinition = addFunction(libDef, ...
    "int reactorsurface_install(int i,int n)", ...
    "MATLABName", "clib.ctMatlabInterface.reactorsurface_install", ...
    "Description", "clib.ctMatlabInterface.reactorsurface_install Representation of C++ function reactorsurface_install."); % Modify help description values as needed.
defineArgument(reactorsurface_installDefinition, "i", "int32");
defineArgument(reactorsurface_installDefinition, "n", "int32");
defineOutput(reactorsurface_installDefinition, "RetVal", "int32");
validate(reactorsurface_installDefinition);

%% C++ function |reactorsurface_setkinetics| with MATLAB name |clib.ctMatlabInterface.reactorsurface_setkinetics|
% C++ Signature: int reactorsurface_setkinetics(int i,int n)

reactorsurface_setkineticsDefinition = addFunction(libDef, ...
    "int reactorsurface_setkinetics(int i,int n)", ...
    "MATLABName", "clib.ctMatlabInterface.reactorsurface_setkinetics", ...
    "Description", "clib.ctMatlabInterface.reactorsurface_setkinetics Representation of C++ function reactorsurface_setkinetics."); % Modify help description values as needed.
defineArgument(reactorsurface_setkineticsDefinition, "i", "int32");
defineArgument(reactorsurface_setkineticsDefinition, "n", "int32");
defineOutput(reactorsurface_setkineticsDefinition, "RetVal", "int32");
validate(reactorsurface_setkineticsDefinition);

%% C++ function |reactorsurface_area| with MATLAB name |clib.ctMatlabInterface.reactorsurface_area|
% C++ Signature: double reactorsurface_area(int i)

reactorsurface_areaDefinition = addFunction(libDef, ...
    "double reactorsurface_area(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.reactorsurface_area", ...
    "Description", "clib.ctMatlabInterface.reactorsurface_area Representation of C++ function reactorsurface_area."); % Modify help description values as needed.
defineArgument(reactorsurface_areaDefinition, "i", "int32");
defineOutput(reactorsurface_areaDefinition, "RetVal", "double");
validate(reactorsurface_areaDefinition);

%% C++ function |reactorsurface_setArea| with MATLAB name |clib.ctMatlabInterface.reactorsurface_setArea|
% C++ Signature: int reactorsurface_setArea(int i,double v)

reactorsurface_setAreaDefinition = addFunction(libDef, ...
    "int reactorsurface_setArea(int i,double v)", ...
    "MATLABName", "clib.ctMatlabInterface.reactorsurface_setArea", ...
    "Description", "clib.ctMatlabInterface.reactorsurface_setArea Representation of C++ function reactorsurface_setArea."); % Modify help description values as needed.
defineArgument(reactorsurface_setAreaDefinition, "i", "int32");
defineArgument(reactorsurface_setAreaDefinition, "v", "double");
defineOutput(reactorsurface_setAreaDefinition, "RetVal", "int32");
validate(reactorsurface_setAreaDefinition);

%% C++ function |reactorsurface_addSensitivityReaction| with MATLAB name |clib.ctMatlabInterface.reactorsurface_addSensitivityReaction|
% C++ Signature: int reactorsurface_addSensitivityReaction(int i,int rxn)

reactorsurface_addSensitivityReactionDefinition = addFunction(libDef, ...
    "int reactorsurface_addSensitivityReaction(int i,int rxn)", ...
    "MATLABName", "clib.ctMatlabInterface.reactorsurface_addSensitivityReaction", ...
    "Description", "clib.ctMatlabInterface.reactorsurface_addSensitivityReaction Representation of C++ function reactorsurface_addSensitivityReaction."); % Modify help description values as needed.
defineArgument(reactorsurface_addSensitivityReactionDefinition, "i", "int32");
defineArgument(reactorsurface_addSensitivityReactionDefinition, "rxn", "int32");
defineOutput(reactorsurface_addSensitivityReactionDefinition, "RetVal", "int32");
validate(reactorsurface_addSensitivityReactionDefinition);

%% C++ function |ct_clearReactors| with MATLAB name |clib.ctMatlabInterface.ct_clearReactors|
% C++ Signature: int ct_clearReactors()

ct_clearReactorsDefinition = addFunction(libDef, ...
    "int ct_clearReactors()", ...
    "MATLABName", "clib.ctMatlabInterface.ct_clearReactors", ...
    "Description", "clib.ctMatlabInterface.ct_clearReactors Representation of C++ function ct_clearReactors."); % Modify help description values as needed.
defineOutput(ct_clearReactorsDefinition, "RetVal", "int32");
validate(ct_clearReactorsDefinition);

%% C++ function |rdiag_new| with MATLAB name |clib.ctMatlabInterface.rdiag_new|
% C++ Signature: int rdiag_new()

rdiag_newDefinition = addFunction(libDef, ...
    "int rdiag_new()", ...
    "MATLABName", "clib.ctMatlabInterface.rdiag_new", ...
    "Description", "clib.ctMatlabInterface.rdiag_new Representation of C++ function rdiag_new."); % Modify help description values as needed.
defineOutput(rdiag_newDefinition, "RetVal", "int32");
validate(rdiag_newDefinition);

%% C++ function |rdiag_del| with MATLAB name |clib.ctMatlabInterface.rdiag_del|
% C++ Signature: int rdiag_del(int i)

rdiag_delDefinition = addFunction(libDef, ...
    "int rdiag_del(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.rdiag_del", ...
    "Description", "clib.ctMatlabInterface.rdiag_del Representation of C++ function rdiag_del."); % Modify help description values as needed.
defineArgument(rdiag_delDefinition, "i", "int32");
defineOutput(rdiag_delDefinition, "RetVal", "int32");
validate(rdiag_delDefinition);

%% C++ function |rdiag_detailed| with MATLAB name |clib.ctMatlabInterface.rdiag_detailed|
% C++ Signature: int rdiag_detailed(int i)

rdiag_detailedDefinition = addFunction(libDef, ...
    "int rdiag_detailed(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.rdiag_detailed", ...
    "Description", "clib.ctMatlabInterface.rdiag_detailed Representation of C++ function rdiag_detailed."); % Modify help description values as needed.
defineArgument(rdiag_detailedDefinition, "i", "int32");
defineOutput(rdiag_detailedDefinition, "RetVal", "int32");
validate(rdiag_detailedDefinition);

%% C++ function |rdiag_brief| with MATLAB name |clib.ctMatlabInterface.rdiag_brief|
% C++ Signature: int rdiag_brief(int i)

rdiag_briefDefinition = addFunction(libDef, ...
    "int rdiag_brief(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.rdiag_brief", ...
    "Description", "clib.ctMatlabInterface.rdiag_brief Representation of C++ function rdiag_brief."); % Modify help description values as needed.
defineArgument(rdiag_briefDefinition, "i", "int32");
defineOutput(rdiag_briefDefinition, "RetVal", "int32");
validate(rdiag_briefDefinition);

%% C++ function |rdiag_setThreshold| with MATLAB name |clib.ctMatlabInterface.rdiag_setThreshold|
% C++ Signature: int rdiag_setThreshold(int i,double v)

rdiag_setThresholdDefinition = addFunction(libDef, ...
    "int rdiag_setThreshold(int i,double v)", ...
    "MATLABName", "clib.ctMatlabInterface.rdiag_setThreshold", ...
    "Description", "clib.ctMatlabInterface.rdiag_setThreshold Representation of C++ function rdiag_setThreshold."); % Modify help description values as needed.
defineArgument(rdiag_setThresholdDefinition, "i", "int32");
defineArgument(rdiag_setThresholdDefinition, "v", "double");
defineOutput(rdiag_setThresholdDefinition, "RetVal", "int32");
validate(rdiag_setThresholdDefinition);

%% C++ function |rdiag_setBoldColor| with MATLAB name |clib.ctMatlabInterface.rdiag_setBoldColor|
% C++ Signature: int rdiag_setBoldColor(int i,char const * color)

rdiag_setBoldColorDefinition = addFunction(libDef, ...
    "int rdiag_setBoldColor(int i,char const * color)", ...
    "MATLABName", "clib.ctMatlabInterface.rdiag_setBoldColor", ...
    "Description", "clib.ctMatlabInterface.rdiag_setBoldColor Representation of C++ function rdiag_setBoldColor."); % Modify help description values as needed.
defineArgument(rdiag_setBoldColorDefinition, "i", "int32");
defineArgument(rdiag_setBoldColorDefinition, "color", "string", "input", "nullTerminated");
defineOutput(rdiag_setBoldColorDefinition, "RetVal", "int32");
validate(rdiag_setBoldColorDefinition);

%% C++ function |rdiag_setNormalColor| with MATLAB name |clib.ctMatlabInterface.rdiag_setNormalColor|
% C++ Signature: int rdiag_setNormalColor(int i,char const * color)

rdiag_setNormalColorDefinition = addFunction(libDef, ...
    "int rdiag_setNormalColor(int i,char const * color)", ...
    "MATLABName", "clib.ctMatlabInterface.rdiag_setNormalColor", ...
    "Description", "clib.ctMatlabInterface.rdiag_setNormalColor Representation of C++ function rdiag_setNormalColor."); % Modify help description values as needed.
defineArgument(rdiag_setNormalColorDefinition, "i", "int32");
defineArgument(rdiag_setNormalColorDefinition, "color", "string", "input", "nullTerminated");
defineOutput(rdiag_setNormalColorDefinition, "RetVal", "int32");
validate(rdiag_setNormalColorDefinition);

%% C++ function |rdiag_setDashedColor| with MATLAB name |clib.ctMatlabInterface.rdiag_setDashedColor|
% C++ Signature: int rdiag_setDashedColor(int i,char const * color)

rdiag_setDashedColorDefinition = addFunction(libDef, ...
    "int rdiag_setDashedColor(int i,char const * color)", ...
    "MATLABName", "clib.ctMatlabInterface.rdiag_setDashedColor", ...
    "Description", "clib.ctMatlabInterface.rdiag_setDashedColor Representation of C++ function rdiag_setDashedColor."); % Modify help description values as needed.
defineArgument(rdiag_setDashedColorDefinition, "i", "int32");
defineArgument(rdiag_setDashedColorDefinition, "color", "string", "input", "nullTerminated");
defineOutput(rdiag_setDashedColorDefinition, "RetVal", "int32");
validate(rdiag_setDashedColorDefinition);

%% C++ function |rdiag_setDotOptions| with MATLAB name |clib.ctMatlabInterface.rdiag_setDotOptions|
% C++ Signature: int rdiag_setDotOptions(int i,char const * opt)

rdiag_setDotOptionsDefinition = addFunction(libDef, ...
    "int rdiag_setDotOptions(int i,char const * opt)", ...
    "MATLABName", "clib.ctMatlabInterface.rdiag_setDotOptions", ...
    "Description", "clib.ctMatlabInterface.rdiag_setDotOptions Representation of C++ function rdiag_setDotOptions."); % Modify help description values as needed.
defineArgument(rdiag_setDotOptionsDefinition, "i", "int32");
defineArgument(rdiag_setDotOptionsDefinition, "opt", "string", "input", "nullTerminated");
defineOutput(rdiag_setDotOptionsDefinition, "RetVal", "int32");
validate(rdiag_setDotOptionsDefinition);

%% C++ function |rdiag_setBoldThreshold| with MATLAB name |clib.ctMatlabInterface.rdiag_setBoldThreshold|
% C++ Signature: int rdiag_setBoldThreshold(int i,double v)

rdiag_setBoldThresholdDefinition = addFunction(libDef, ...
    "int rdiag_setBoldThreshold(int i,double v)", ...
    "MATLABName", "clib.ctMatlabInterface.rdiag_setBoldThreshold", ...
    "Description", "clib.ctMatlabInterface.rdiag_setBoldThreshold Representation of C++ function rdiag_setBoldThreshold."); % Modify help description values as needed.
defineArgument(rdiag_setBoldThresholdDefinition, "i", "int32");
defineArgument(rdiag_setBoldThresholdDefinition, "v", "double");
defineOutput(rdiag_setBoldThresholdDefinition, "RetVal", "int32");
validate(rdiag_setBoldThresholdDefinition);

%% C++ function |rdiag_setNormalThreshold| with MATLAB name |clib.ctMatlabInterface.rdiag_setNormalThreshold|
% C++ Signature: int rdiag_setNormalThreshold(int i,double v)

rdiag_setNormalThresholdDefinition = addFunction(libDef, ...
    "int rdiag_setNormalThreshold(int i,double v)", ...
    "MATLABName", "clib.ctMatlabInterface.rdiag_setNormalThreshold", ...
    "Description", "clib.ctMatlabInterface.rdiag_setNormalThreshold Representation of C++ function rdiag_setNormalThreshold."); % Modify help description values as needed.
defineArgument(rdiag_setNormalThresholdDefinition, "i", "int32");
defineArgument(rdiag_setNormalThresholdDefinition, "v", "double");
defineOutput(rdiag_setNormalThresholdDefinition, "RetVal", "int32");
validate(rdiag_setNormalThresholdDefinition);

%% C++ function |rdiag_setLabelThreshold| with MATLAB name |clib.ctMatlabInterface.rdiag_setLabelThreshold|
% C++ Signature: int rdiag_setLabelThreshold(int i,double v)

rdiag_setLabelThresholdDefinition = addFunction(libDef, ...
    "int rdiag_setLabelThreshold(int i,double v)", ...
    "MATLABName", "clib.ctMatlabInterface.rdiag_setLabelThreshold", ...
    "Description", "clib.ctMatlabInterface.rdiag_setLabelThreshold Representation of C++ function rdiag_setLabelThreshold."); % Modify help description values as needed.
defineArgument(rdiag_setLabelThresholdDefinition, "i", "int32");
defineArgument(rdiag_setLabelThresholdDefinition, "v", "double");
defineOutput(rdiag_setLabelThresholdDefinition, "RetVal", "int32");
validate(rdiag_setLabelThresholdDefinition);

%% C++ function |rdiag_setScale| with MATLAB name |clib.ctMatlabInterface.rdiag_setScale|
% C++ Signature: int rdiag_setScale(int i,double v)

rdiag_setScaleDefinition = addFunction(libDef, ...
    "int rdiag_setScale(int i,double v)", ...
    "MATLABName", "clib.ctMatlabInterface.rdiag_setScale", ...
    "Description", "clib.ctMatlabInterface.rdiag_setScale Representation of C++ function rdiag_setScale."); % Modify help description values as needed.
defineArgument(rdiag_setScaleDefinition, "i", "int32");
defineArgument(rdiag_setScaleDefinition, "v", "double");
defineOutput(rdiag_setScaleDefinition, "RetVal", "int32");
validate(rdiag_setScaleDefinition);

%% C++ function |rdiag_setFlowType| with MATLAB name |clib.ctMatlabInterface.rdiag_setFlowType|
% C++ Signature: int rdiag_setFlowType(int i,int iflow)

rdiag_setFlowTypeDefinition = addFunction(libDef, ...
    "int rdiag_setFlowType(int i,int iflow)", ...
    "MATLABName", "clib.ctMatlabInterface.rdiag_setFlowType", ...
    "Description", "clib.ctMatlabInterface.rdiag_setFlowType Representation of C++ function rdiag_setFlowType."); % Modify help description values as needed.
defineArgument(rdiag_setFlowTypeDefinition, "i", "int32");
defineArgument(rdiag_setFlowTypeDefinition, "iflow", "int32");
defineOutput(rdiag_setFlowTypeDefinition, "RetVal", "int32");
validate(rdiag_setFlowTypeDefinition);

%% C++ function |rdiag_setArrowWidth| with MATLAB name |clib.ctMatlabInterface.rdiag_setArrowWidth|
% C++ Signature: int rdiag_setArrowWidth(int i,double v)

rdiag_setArrowWidthDefinition = addFunction(libDef, ...
    "int rdiag_setArrowWidth(int i,double v)", ...
    "MATLABName", "clib.ctMatlabInterface.rdiag_setArrowWidth", ...
    "Description", "clib.ctMatlabInterface.rdiag_setArrowWidth Representation of C++ function rdiag_setArrowWidth."); % Modify help description values as needed.
defineArgument(rdiag_setArrowWidthDefinition, "i", "int32");
defineArgument(rdiag_setArrowWidthDefinition, "v", "double");
defineOutput(rdiag_setArrowWidthDefinition, "RetVal", "int32");
validate(rdiag_setArrowWidthDefinition);

%% C++ function |rdiag_setTitle| with MATLAB name |clib.ctMatlabInterface.rdiag_setTitle|
% C++ Signature: int rdiag_setTitle(int i,char const * title)

rdiag_setTitleDefinition = addFunction(libDef, ...
    "int rdiag_setTitle(int i,char const * title)", ...
    "MATLABName", "clib.ctMatlabInterface.rdiag_setTitle", ...
    "Description", "clib.ctMatlabInterface.rdiag_setTitle Representation of C++ function rdiag_setTitle."); % Modify help description values as needed.
defineArgument(rdiag_setTitleDefinition, "i", "int32");
defineArgument(rdiag_setTitleDefinition, "title", "string", "input", "nullTerminated");
defineOutput(rdiag_setTitleDefinition, "RetVal", "int32");
validate(rdiag_setTitleDefinition);

%% C++ function |rdiag_write| with MATLAB name |clib.ctMatlabInterface.rdiag_write|
% C++ Signature: int rdiag_write(int i,int fmt,char const * fname)

rdiag_writeDefinition = addFunction(libDef, ...
    "int rdiag_write(int i,int fmt,char const * fname)", ...
    "MATLABName", "clib.ctMatlabInterface.rdiag_write", ...
    "Description", "clib.ctMatlabInterface.rdiag_write Representation of C++ function rdiag_write."); % Modify help description values as needed.
defineArgument(rdiag_writeDefinition, "i", "int32");
defineArgument(rdiag_writeDefinition, "fmt", "int32");
defineArgument(rdiag_writeDefinition, "fname", "string", "input", "nullTerminated");
defineOutput(rdiag_writeDefinition, "RetVal", "int32");
validate(rdiag_writeDefinition);

%% C++ function |rdiag_add| with MATLAB name |clib.ctMatlabInterface.rdiag_add|
% C++ Signature: int rdiag_add(int i,int n)

rdiag_addDefinition = addFunction(libDef, ...
    "int rdiag_add(int i,int n)", ...
    "MATLABName", "clib.ctMatlabInterface.rdiag_add", ...
    "Description", "clib.ctMatlabInterface.rdiag_add Representation of C++ function rdiag_add."); % Modify help description values as needed.
defineArgument(rdiag_addDefinition, "i", "int32");
defineArgument(rdiag_addDefinition, "n", "int32");
defineOutput(rdiag_addDefinition, "RetVal", "int32");
validate(rdiag_addDefinition);

%% C++ function |rdiag_findMajor| with MATLAB name |clib.ctMatlabInterface.rdiag_findMajor|
% C++ Signature: int rdiag_findMajor(int i,double threshold,size_t lda,double * a)

rdiag_findMajorDefinition = addFunction(libDef, ...
    "int rdiag_findMajor(int i,double threshold,size_t lda,double * a)", ...
    "MATLABName", "clib.ctMatlabInterface.rdiag_findMajor", ...
    "Description", "clib.ctMatlabInterface.rdiag_findMajor Representation of C++ function rdiag_findMajor."); % Modify help description values as needed.
defineArgument(rdiag_findMajorDefinition, "i", "int32");
defineArgument(rdiag_findMajorDefinition, "threshold", "double");
defineArgument(rdiag_findMajorDefinition, "lda", "uint64");
defineArgument(rdiag_findMajorDefinition, "a", "clib.array.ctMatlabInterface.Double", "input", "lda"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(rdiag_findMajorDefinition, "RetVal", "int32");
validate(rdiag_findMajorDefinition);

%% C++ function |rdiag_setFont| with MATLAB name |clib.ctMatlabInterface.rdiag_setFont|
% C++ Signature: int rdiag_setFont(int i,char const * font)

rdiag_setFontDefinition = addFunction(libDef, ...
    "int rdiag_setFont(int i,char const * font)", ...
    "MATLABName", "clib.ctMatlabInterface.rdiag_setFont", ...
    "Description", "clib.ctMatlabInterface.rdiag_setFont Representation of C++ function rdiag_setFont."); % Modify help description values as needed.
defineArgument(rdiag_setFontDefinition, "i", "int32");
defineArgument(rdiag_setFontDefinition, "font", "string", "input", "nullTerminated");
defineOutput(rdiag_setFontDefinition, "RetVal", "int32");
validate(rdiag_setFontDefinition);

%% C++ function |rdiag_displayOnly| with MATLAB name |clib.ctMatlabInterface.rdiag_displayOnly|
% C++ Signature: int rdiag_displayOnly(int i,int k)

rdiag_displayOnlyDefinition = addFunction(libDef, ...
    "int rdiag_displayOnly(int i,int k)", ...
    "MATLABName", "clib.ctMatlabInterface.rdiag_displayOnly", ...
    "Description", "clib.ctMatlabInterface.rdiag_displayOnly Representation of C++ function rdiag_displayOnly."); % Modify help description values as needed.
defineArgument(rdiag_displayOnlyDefinition, "i", "int32");
defineArgument(rdiag_displayOnlyDefinition, "k", "int32");
defineOutput(rdiag_displayOnlyDefinition, "RetVal", "int32");
validate(rdiag_displayOnlyDefinition);

%% C++ function |rbuild_new| with MATLAB name |clib.ctMatlabInterface.rbuild_new|
% C++ Signature: int rbuild_new()

rbuild_newDefinition = addFunction(libDef, ...
    "int rbuild_new()", ...
    "MATLABName", "clib.ctMatlabInterface.rbuild_new", ...
    "Description", "clib.ctMatlabInterface.rbuild_new Representation of C++ function rbuild_new."); % Modify help description values as needed.
defineOutput(rbuild_newDefinition, "RetVal", "int32");
validate(rbuild_newDefinition);

%% C++ function |rbuild_del| with MATLAB name |clib.ctMatlabInterface.rbuild_del|
% C++ Signature: int rbuild_del(int i)

rbuild_delDefinition = addFunction(libDef, ...
    "int rbuild_del(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.rbuild_del", ...
    "Description", "clib.ctMatlabInterface.rbuild_del Representation of C++ function rbuild_del."); % Modify help description values as needed.
defineArgument(rbuild_delDefinition, "i", "int32");
defineOutput(rbuild_delDefinition, "RetVal", "int32");
validate(rbuild_delDefinition);

%% C++ function |rbuild_init| with MATLAB name |clib.ctMatlabInterface.rbuild_init|
% C++ Signature: int rbuild_init(int i,char const * logfile,int k)

rbuild_initDefinition = addFunction(libDef, ...
    "int rbuild_init(int i,char const * logfile,int k)", ...
    "MATLABName", "clib.ctMatlabInterface.rbuild_init", ...
    "Description", "clib.ctMatlabInterface.rbuild_init Representation of C++ function rbuild_init."); % Modify help description values as needed.
defineArgument(rbuild_initDefinition, "i", "int32");
defineArgument(rbuild_initDefinition, "logfile", "string", "input", "nullTerminated");
defineArgument(rbuild_initDefinition, "k", "int32");
defineOutput(rbuild_initDefinition, "RetVal", "int32");
validate(rbuild_initDefinition);

%% C++ function |rbuild_build| with MATLAB name |clib.ctMatlabInterface.rbuild_build|
% C++ Signature: int rbuild_build(int i,int k,char const * el,char const * dotfile,int idiag,int iquiet)

rbuild_buildDefinition = addFunction(libDef, ...
    "int rbuild_build(int i,int k,char const * el,char const * dotfile,int idiag,int iquiet)", ...
    "MATLABName", "clib.ctMatlabInterface.rbuild_build", ...
    "Description", "clib.ctMatlabInterface.rbuild_build Representation of C++ function rbuild_build."); % Modify help description values as needed.
defineArgument(rbuild_buildDefinition, "i", "int32");
defineArgument(rbuild_buildDefinition, "k", "int32");
defineArgument(rbuild_buildDefinition, "el", "string", "input", "nullTerminated");
defineArgument(rbuild_buildDefinition, "dotfile", "string", "input", "nullTerminated");
defineArgument(rbuild_buildDefinition, "idiag", "int32");
defineArgument(rbuild_buildDefinition, "iquiet", "int32");
defineOutput(rbuild_buildDefinition, "RetVal", "int32");
validate(rbuild_buildDefinition);

%% C++ function |ct_clearReactionPath| with MATLAB name |clib.ctMatlabInterface.ct_clearReactionPath|
% C++ Signature: int ct_clearReactionPath()

ct_clearReactionPathDefinition = addFunction(libDef, ...
    "int ct_clearReactionPath()", ...
    "MATLABName", "clib.ctMatlabInterface.ct_clearReactionPath", ...
    "Description", "clib.ctMatlabInterface.ct_clearReactionPath Representation of C++ function ct_clearReactionPath."); % Modify help description values as needed.
defineOutput(ct_clearReactionPathDefinition, "RetVal", "int32");
validate(ct_clearReactionPathDefinition);

%% C++ function |surf_setCoverages| with MATLAB name |clib.ctMatlabInterface.surf_setCoverages|
% C++ Signature: int surf_setCoverages(int i,double const * c,int norm)

surf_setCoveragesDefinition = addFunction(libDef, ...
    "int surf_setCoverages(int i,double const * c,int norm)", ...
    "MATLABName", "clib.ctMatlabInterface.surf_setCoverages", ...
    "Description", "clib.ctMatlabInterface.surf_setCoverages Representation of C++ function surf_setCoverages."); % Modify help description values as needed.
defineArgument(surf_setCoveragesDefinition, "i", "int32");
defineArgument(surf_setCoveragesDefinition, "c", "clib.array.ctMatlabInterface.Double", "input", "i"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineArgument(surf_setCoveragesDefinition, "norm", "int32");
defineOutput(surf_setCoveragesDefinition, "RetVal", "int32");
validate(surf_setCoveragesDefinition);

%% C++ function |surf_getCoverages| with MATLAB name |clib.ctMatlabInterface.surf_getCoverages|
% C++ Signature: int surf_getCoverages(int i,double * c)

surf_getCoveragesDefinition = addFunction(libDef, ...
    "int surf_getCoverages(int i,double * c)", ...
    "MATLABName", "clib.ctMatlabInterface.surf_getCoverages", ...
    "Description", "clib.ctMatlabInterface.surf_getCoverages Representation of C++ function surf_getCoverages."); % Modify help description values as needed.
defineArgument(surf_getCoveragesDefinition, "i", "int32");
defineArgument(surf_getCoveragesDefinition, "c", "clib.array.ctMatlabInterface.Double", "input", "i"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(surf_getCoveragesDefinition, "RetVal", "int32");
validate(surf_getCoveragesDefinition);

%% C++ function |surf_setConcentrations| with MATLAB name |clib.ctMatlabInterface.surf_setConcentrations|
% C++ Signature: int surf_setConcentrations(int i,double const * c)

surf_setConcentrationsDefinition = addFunction(libDef, ...
    "int surf_setConcentrations(int i,double const * c)", ...
    "MATLABName", "clib.ctMatlabInterface.surf_setConcentrations", ...
    "Description", "clib.ctMatlabInterface.surf_setConcentrations Representation of C++ function surf_setConcentrations."); % Modify help description values as needed.
defineArgument(surf_setConcentrationsDefinition, "i", "int32");
defineArgument(surf_setConcentrationsDefinition, "c", "clib.array.ctMatlabInterface.Double", "input", "i"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(surf_setConcentrationsDefinition, "RetVal", "int32");
validate(surf_setConcentrationsDefinition);

%% C++ function |surf_getConcentrations| with MATLAB name |clib.ctMatlabInterface.surf_getConcentrations|
% C++ Signature: int surf_getConcentrations(int i,double * c)

surf_getConcentrationsDefinition = addFunction(libDef, ...
    "int surf_getConcentrations(int i,double * c)", ...
    "MATLABName", "clib.ctMatlabInterface.surf_getConcentrations", ...
    "Description", "clib.ctMatlabInterface.surf_getConcentrations Representation of C++ function surf_getConcentrations."); % Modify help description values as needed.
defineArgument(surf_getConcentrationsDefinition, "i", "int32");
defineArgument(surf_getConcentrationsDefinition, "c", "clib.array.ctMatlabInterface.Double", "input", "i"); % <MLTYPE> can be "clib.array.ctMatlabInterface.Double", or "double"
defineOutput(surf_getConcentrationsDefinition, "RetVal", "int32");
validate(surf_getConcentrationsDefinition);

%% C++ function |surf_setSiteDensity| with MATLAB name |clib.ctMatlabInterface.surf_setSiteDensity|
% C++ Signature: int surf_setSiteDensity(int i,double s0)

surf_setSiteDensityDefinition = addFunction(libDef, ...
    "int surf_setSiteDensity(int i,double s0)", ...
    "MATLABName", "clib.ctMatlabInterface.surf_setSiteDensity", ...
    "Description", "clib.ctMatlabInterface.surf_setSiteDensity Representation of C++ function surf_setSiteDensity."); % Modify help description values as needed.
defineArgument(surf_setSiteDensityDefinition, "i", "int32");
defineArgument(surf_setSiteDensityDefinition, "s0", "double");
defineOutput(surf_setSiteDensityDefinition, "RetVal", "int32");
validate(surf_setSiteDensityDefinition);

%% C++ function |surf_siteDensity| with MATLAB name |clib.ctMatlabInterface.surf_siteDensity|
% C++ Signature: double surf_siteDensity(int i)

surf_siteDensityDefinition = addFunction(libDef, ...
    "double surf_siteDensity(int i)", ...
    "MATLABName", "clib.ctMatlabInterface.surf_siteDensity", ...
    "Description", "clib.ctMatlabInterface.surf_siteDensity Representation of C++ function surf_siteDensity."); % Modify help description values as needed.
defineArgument(surf_siteDensityDefinition, "i", "int32");
defineOutput(surf_siteDensityDefinition, "RetVal", "double");
validate(surf_siteDensityDefinition);

%% C++ function |surf_setCoveragesByName| with MATLAB name |clib.ctMatlabInterface.surf_setCoveragesByName|
% C++ Signature: int surf_setCoveragesByName(int i,char const * c)

surf_setCoveragesByNameDefinition = addFunction(libDef, ...
    "int surf_setCoveragesByName(int i,char const * c)", ...
    "MATLABName", "clib.ctMatlabInterface.surf_setCoveragesByName", ...
    "Description", "clib.ctMatlabInterface.surf_setCoveragesByName Representation of C++ function surf_setCoveragesByName."); % Modify help description values as needed.
defineArgument(surf_setCoveragesByNameDefinition, "i", "int32");
defineArgument(surf_setCoveragesByNameDefinition, "c", "string", "input", "nullTerminated");
defineOutput(surf_setCoveragesByNameDefinition, "RetVal", "int32");
validate(surf_setCoveragesByNameDefinition);

%% Validate the library definition
validate(libDef);

end
