function nm = elementName(a, m)
% ELEMENTNAME -  Name of element with index m.
%
%    If m is a scalar integer, the return value will be a string
%    containing the name of the m^th species. If it is an array of
%    integers, the output will be a cell array of
%    the same shape containing the name strings.
%
[m, n] = size(k);
nm = {};
for i = 1:m
  for j = 1:n
    nm{i,j} = phase_get(a.tp_id, 41, k(i,j));
  end
end


