function t = ThermoPhase(src, id)
% THERMOPHASE  ThermoPhase class constructor.
% t = ThermoPhase(src, id)
% :param src:
%     Input string of YAML file name.
% :param id:
%     ID of the phase to import as specified in the input file. (optional)
% :return:
%     Instance of class :mat:func:`ThermoPhase`
%

if nargin > 2
    error('ThermoPhase expects 1 or 2 input arguments.');
end

if nargin == 1
    id = '-';
end

t.owner = 1;
t.tp_id = thermo_get(0, 0, src, id);
if t.tp_id < 0
    error(geterr);
end
t = class(t, 'ThermoPhase');
