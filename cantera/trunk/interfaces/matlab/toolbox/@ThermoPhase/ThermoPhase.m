function t = ThermoPhase(r)
% THERMOPHASE  ThermoPhase class constructor.
% t = ThermoPhase(r)
% :param r:
%     If ``r`` is an instance of class :mat:func:`ThermoPhase`,
%     a copy of the instance is returned. Otherwise, ``r`` must
%     be an instance of class :mat:func:`XML_Node`.
% :return:
%     Instance of class :mat:func:`ThermoPhase`
%

if nargin == 1
    if isa(r, 'ThermoPhase')
        % create a copy
        t = r;
        return
    elseif isa(r, 'XML_Node')
        t.owner = 1;
        hr = xml_hndl(r);
        t.tp_id = thermo_get(hr, 0);
        if t.tp_id < 0
            error(geterr);
        end
    else
        t.owner = 0;
        t.tp_id = r;
    end
    t = class(t, 'ThermoPhase');
else
    error('ThermoPhase expects 1 input argument.');
end
