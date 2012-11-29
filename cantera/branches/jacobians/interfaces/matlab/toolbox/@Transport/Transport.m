function tr = Transport(xml_phase, th, model, loglevel)
%TRANSPORT  Transport class constructor.
%
%  k = TRANSPORT(model, p, loglevel) creates a transport
%  manager for phase object p.
%
%  The 'model' parameter is a string that specifies the transport
%  model. The phase object must have already been created.
%
tr.id = 0;
if nargin == 4
    tr.th = th;
    if strcmp(model, 'default')
        try
            node = child(xml_phase,'transport');
            tr.model = attrib(node,'model');
        catch
            tr.model = 'None';
        end
    else
        tr.model = model;
    end
    tr.id = trans_get(hndl(th), -1, tr.model, loglevel) ;
    tr = class(tr,'Transport');
elseif isa(model,'Transport')
    tr = model;
else
    error('syntax error');
end
