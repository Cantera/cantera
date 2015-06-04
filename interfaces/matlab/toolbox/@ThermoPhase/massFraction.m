function y = massFraction(tp, species)
% MASSFRACTION  Get the mass fraction of a species.
% y = massFraction(tp, species)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :param species:
%     String or cell array of strings of species whose mass
%     fraction is desired
% :return:
%     Scalar or vector double mass fractions
%

y = 0.0;
yarray = massFractions(tp);
if isa(species, 'char')
    k = speciesIndex(tp, species);
    if  k > 0
        y = yarray(k);
    end

elseif isa(species, 'cell')
    n = length(species);
    y = zeros(1, n);
    for j = 1:n
        k = speciesIndex(tp, species{j});
        if k > 0
            y(j) = yarray(k);
        end
    end
end
