function nm = speciesName(tp, k)
% SPECIESNAME  Get the name of a species given the index.
% nm = speciesName(tp, k)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     class derived from ThermoPhase)
% :param k:
%     Scalar or array of integer species numbers
% :return:
%     Cell array of strings
%

[m,n] = size(k);
nm = cell(m, n);
for i = 1:m
    for j = 1:n
        nm{i,j} = phase_get(tp.tp_id, 40, k(i,j));
    end
end
