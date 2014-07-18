function m = Inlet(id)
% INLET  Create an inlet domain.
% m = Inlet(id)
% Note that an inlet can only be a terminal domain - it must be
% either the leftmost or rightmost domain in a stack.
%
% :param id:
%     String name of the inlet.
% :return:
%     Instance of class :mat:func:`Domain1D` representing an inlet.
%

m = Domain1D(2);
if nargin == 0
    setID(m, 'inlet');
else
    setID(m, id);
end
