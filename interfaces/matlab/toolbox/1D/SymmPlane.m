function m = SymmPlane(id)
% SYMMPLANE  Create a symmetry plane domain.
% m = SymmPlane(id)
% :param id:
%     String ID of the symmetry plane.
% :return:
%     Instance of class :mat:func:`Domain1D` representing a symmetry
%     plane.
%

m = Domain1D(4);
if nargin == 0
    setID(m, 'symmetry_plane');
else
    setID(m, id);
end
