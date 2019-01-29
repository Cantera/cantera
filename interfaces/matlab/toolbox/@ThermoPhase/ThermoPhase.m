function t = ThermoPhase(r)
% THERMOPHASE  ThermoPhase class constructor.
% t = ThermoPhase(r)
% :param r:
%     An instance of class :mat:func:`XML_Node`.
% :return:
%     Instance of class :mat:func:`ThermoPhase`
%

if nargin ~= 1
    error('ThermoPhase expects 1 input argument.');
end

t.owner = 1;
hr = xml_hndl(r);
t.tp_id = thermo_get(hr, 0);
if t.tp_id < 0
    error(geterr);
end
t = class(t, 'ThermoPhase');
