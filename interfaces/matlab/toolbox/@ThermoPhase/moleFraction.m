function x = moleFraction(tp, species)

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
