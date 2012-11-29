function m = Outlet(id)
% OUTLET - Return a Domain1D instance representing an outlet.
%
m = Domain1D(5);
if nargin == 0
    setID(m,'outlet');
else
    setID(m,id);
end
