function m = Surface(id)
% SURFACE - Return a Domain1D instance representing a non-reacting
% surface.
m = Domain1D(3);
if nargin == 0
  setID(m,'surface');
else
  setID(m,id);
end
