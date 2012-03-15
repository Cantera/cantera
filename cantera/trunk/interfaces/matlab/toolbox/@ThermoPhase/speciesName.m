function nm = speciesName(a, k)
% SPECIESNAME - Name of species k.
%    If k is a scalar integer, the return value will be a string
%    containing the name of the kth species. If it is an array of
%    integers, the output will be a cell array of
%    the same shape containing the name strings.
%
[m, n] = size(k);
nm = cell(m,n);
for i = 1:m
    for j = 1:n
        nm{i,j} = phase_get(a.tp_id, 40, k(i,j));
    end
end
