function m = Surface(id, surface_mech)
% SURFACE - Return a Domain1D instance representing a non-reacting
% or reacting surface.
if nargin < 2
    m = Domain1D(3);
    if nargin == 0
        setID(m,'surface');
    elseif nargin == 1
        setID(m,id);
    end
else
    m = Domain1D(6, surface_mech);
    setID(m,id);
end
