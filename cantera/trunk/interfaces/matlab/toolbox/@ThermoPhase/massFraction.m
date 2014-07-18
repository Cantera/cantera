function y = massFraction(tp, species)

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
