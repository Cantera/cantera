function m = SymmPlane(id)
% SYMMPLANE - Return a Domain1D instance representing a symmetry plane.
%
m = Domain1D(4);
if nargin == 0
    setID(m,'symmetry_plane');
else
    setID(m,id);
end
