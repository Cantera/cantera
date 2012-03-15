function t = ThermoPhase(r)
%THERMOPHASE Cantera ThermoPhase class constructor
%
if nargin == 1
    if isa(r,'ThermoPhase')
        % create a copy
        t = r;
        return
    elseif isa(r,'XML_Node')
        t.owner = 1;
        hr = hndl(r);
        t.tp_id = thermo_get(hr,0);
        if t.tp_id < 0
            error(geterr);
        end
    else
        t.owner = 0;
        t.tp_id = r;
    end
    t = class(t,'ThermoPhase');
else
    error('wrong number of arguments');
end
