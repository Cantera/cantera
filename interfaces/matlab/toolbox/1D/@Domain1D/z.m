function zz = z(d, n)
% GRID - 
%   
if nargin == 1
  for i = 1:nPoints(d)
    zz(i) = domain_methods(d.dom_id, 19, i);
  end
else
  m = length(n);
  for i = 1:m
    zz(i) = domain_methods(d.dom_id, 19, n(i));
  end
end


