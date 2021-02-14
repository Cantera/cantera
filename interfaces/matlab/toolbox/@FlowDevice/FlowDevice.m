function x = FlowDevice(typ)
% FLOWDEVICE  FlowDevice class constructor.
% x = FlowDevice(typ)
% Base class for devices that allow flow between reactors.
% :mat:func:`FlowDevice` objects are assumed to be adiabatic,
% non-reactive, and have negligible internal volume, so that they are
% internally always in steady-state even if the upstream and downstream
% reactors are not. The fluid enthalpy, chemical composition, and mass
% flow rate are constant across a :mat:func:`FlowDevice`, and the
% pressure difference equals the difference in pressure between the
% upstream and downstream reactors.
%
% See also: :mat:func:`MassFlowController`, :mat:func:`Valve`
%
% :param typ:
%     Type of :mat:func:`FlowDevice` to be created. ``typ='MassFlowController'``
%     for :mat:func:`MassFlowController`,  ``typ='PressureController'`` for
%     :mat:func:`PressureController` and ``typ='Valve'`` for
%     :mat:func:`Valve`
% :return:
%     Instance of class :mat:func:`FlowDevice`
%

if nargin == 0
    typ = 'MassFlowController';
end

x.type = char(typ);
x.index = flowdevicemethods(0, x.type);
if x.index < 0
    error(geterr);
end
x.upstream = -1;
x.downstream = -1;
x = class(x, 'FlowDevice');
