function d = setTolerances(d, rtol, atol, typ)
% SETTOLERANCES - 
%   
ityp = 0;
if nargin == 4 
  switch typ
   case 'ts'
    itype = -1;
   case 'time'
    itype = -1;
   case 'ss'
    itype = 1;
   case 'steady'
    itype = 1;
  end
end
domain_methods(d.dom_id, 52, rtol, atol, ityp);
