function m = Surface(id, surface_mech)
% SURFACE  Create a surface domain.
% m = Surface(id, surface_mech)
% :param id:
%     String ID of surface
% :param surface_mech:
%     Instance of class :mat:func:`Interface` defining
%     the surface reaction mechanism to be used. Optional.
% :return:
%     Instance of class :mat:func:`Domain1D` representing a
%     non-reacting or reacting surface.
%

if nargin < 2
    m = Domain1D(3);
    if nargin == 0
        setID(m, 'surface');
    elseif nargin == 1
        setID(m, id);
    end
else
    m = Domain1D(6, surface_mech);
    setID(m, id);
end
