function n = speciesNames(tp)
% SPECIESNAMES  Get the species names.
% n = speciesNames(tp)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :return:
%     Cell array of strings of all of the species names
%

n = speciesName(tp, 1:nSpecies(tp));
