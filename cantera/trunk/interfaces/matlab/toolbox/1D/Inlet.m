function m = Inlet(id)
% INLET - Return a Domain1D instance representing an inlet.
%
%   Note that an inlet can only be a terminal domain - it must be
%   either the leftmost or rightmost domain in a stack.
m = Domain1D(2);
if nargin == 0
    setID(m,'inlet');
else
    setID(m,id);
end
