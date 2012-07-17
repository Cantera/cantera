function m = OutletRes(id)
% OUTLET - Return a Domain1D instance representing an outlet reservoir.
%
m = Domain1D(-2);
if nargin == 0
    setID(m,'outletres');
else
    setID(m,id);
end
