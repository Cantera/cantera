function m = AxisymmetricFlow(gas, id)
% AXISYMMETRICFLOW  Create an axisymmetric flow domain.
% m = AxisymmetricFlow(gas, id)
% :param gas:
%     Instance of class :mat:func:`Solution`
% :param id:
%     String, ID of the flow
% :return:
%     Domain1D instance representing an axisymmetric flow.
%

m = Domain1D(1, gas);
if nargin == 1
    setID(m, 'flow');
else
    setID(m, id);
end
