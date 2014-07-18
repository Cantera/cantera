function v = viscosity(a)
% VISCOSITY  Get the dynamic viscosity.
% v = viscosity(a)
% :param a:
%     Instance of class :mat:func:`Transport` (or another
%     object derived from Transport)
%     for which the viscosity is desired.
% :return:
%     Dynamic viscosity. Units: Pa*s
%

v = trans_get(a.id, 1);
if v == -1.0
    error(geterr);
elseif v < 0.0
    error('exception raised');
end
