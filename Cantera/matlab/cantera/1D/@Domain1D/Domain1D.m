function d = Domain1D(a, b, c, d, e)
% DOMAIN1D - Create a new one-dimensional domain. 
%
d.dom_id = -1;

if nargin == 1
  d.dom_id = domain_methods(0, a);
elseif nargin == 2
  if a == 1
    if isa(b,'Solution')
      d.dom_id = domain_methods(0, 1, thermo_hndl(b), kinetics_hndl(b), ...
				trans_hndl(b));
    else
      error('Wrong argument type. Expecting instance of class Solution.')
    end
  end
end
if d.dom_id < 0
  error(geterr);
end
d = class(d, 'Domain1D');
