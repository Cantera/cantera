function x = moleFraction(s, species)

x = 0.0;
xarray = moleFractions(s);
if isa(species,'char')
    k = speciesIndex(s, species);
    if  k > 0
        x = xarray(k);
    end

elseif isa(species,'cell')
    n = length(species);
    x = zeros(1, n);
    for j = 1:n
        k = speciesIndex(s, species{j});
        if k > 0
            x(j) = xarray(k);
        end
    end
end
