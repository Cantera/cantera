function m = AxisymmetricFlow(gas, id)
% AXISYMMETRICFLOW - Axisymmetric flow.
%
%    Return a Domain1D instance representing an axisymmetric flow.
%
m = Domain1D(1, gas);
if nargin == 1
    setID(m,'flow');
else
    setID(m,id);
end
