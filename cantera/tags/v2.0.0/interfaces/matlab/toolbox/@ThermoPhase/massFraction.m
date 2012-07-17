function y = massFraction(s, species)

y = 0.0;
yarray = massFractions(s);
if isa(species,'char')
    k = speciesIndex(s, species);
    if  k > 0
        y = yarray(k);
    end

elseif isa(species,'cell')
    n = length(species);
    y = zeros(1, n);
    for j = 1:n
        k = speciesIndex(s, species{j});
        if k > 0
            y(j) = yarray(k);
        end
    end
end
