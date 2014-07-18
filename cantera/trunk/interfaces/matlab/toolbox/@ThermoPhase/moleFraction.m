function x = moleFraction(tp, species)
% MOLEFRACTION  Get the mole fraction of a species.
% x = moleFraction(tp, species)
% :param tp:
%     Instance of class :mat:func:`ThermoPhase` (or another
%     object that derives from ThermoPhase)
% :param species:
%     String or cell array of strings of species whose mole
%     fraction is desired
% :return:
%     Scalar or vector double mole fractions
%

x = 0.0;
xarray = moleFractions(tp);
if isa(species, 'char')
    k = speciesIndex(tp, species);
    if  k > 0
        x = xarray(k);
    end

elseif isa(species, 'cell')
    n = length(species);
    x = zeros(1, n);
    for j = 1:n
        k = speciesIndex(tp, species{j});
        if k > 0
            x(j) = xarray(k);
        end
    end
end
