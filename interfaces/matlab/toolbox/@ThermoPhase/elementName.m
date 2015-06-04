function nm = elementName(tp, m)
% ELEMENTNAME  Get the name of an element given its index.
% nm = elementName(tp, m)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :param m:
%     Scalar or vector of integers of element indices
% :return:
%     If m is a scalar integer, the return value will be a string
%     containing the name of the m^th species. If it is an array of
%     integers, the output will be a cell array of
%     the same shape containing the name strings.
%

[mm, nn] = size(m);
nm = cell(mm,nn);
for i = 1:mm
    for j = 1:nn
        nm{i,j} = phase_get(tp.tp_id, 41, m(i,j));
    end
end
