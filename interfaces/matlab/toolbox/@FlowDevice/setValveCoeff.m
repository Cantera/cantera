function setValveCoeff(f, k)
% SETVALVECOEFF  Set the valve coefficient :math:`K`.
% setValveCoeff(f, k)
% The mass flow rate [kg/s] is computed from the expression
%
% .. math:: \dot{m} = K(P_{upstream} - P_{downstream})
%
% as long as this produces a positive value.  If this expression is
% negative, zero is returned.
%
% See also: :mat:func:`Valve`
%
% :param f:
%     Instance of class :mat:func:`Valve`
% :param k:
%     Value of the valve coefficient. Units: kg/Pa-s
%

if ~strcmp(f.type, 'Valve')
    error('Valve coefficient can only be set for valves')
end
ok = flowdevicemethods(4, f.index, k);
if ok < 0
    error(geterr);
end
