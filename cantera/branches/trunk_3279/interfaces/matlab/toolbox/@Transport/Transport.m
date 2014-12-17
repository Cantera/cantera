function tr = Transport(r, th, model, loglevel)
% TRANSPORT  Transport class constructor.
% tr = Transport(r, th, model, loglevel)
% Create a new instance of class :mat:func:`Transport`. One, three,
% or four arguments may be supplied. If one argument is given,
% it must be an instance of class :mat:func:`Transport`, and
% a copy will be returned. If three or four arguments are given,
% the first two must be an instance of class :mat:func:`XML_Node`
% and an instance of class :mat:func:`ThermoPhase` respectively.
% The third argument is the type of modeling desired, specified
% by the string ``'default'``, ``'Mix'`` or ``'Multi'``.
% ``'default'`` uses the default transport specified in the
% :mat:func:`XML_Node`. The fourth argument is
% the logging level desired.
%
% :param r:
%     Either an instance of class :mat:func:`Transport` or an
%     instance of class :mat:func:`XML_Node`
% :param th:
%     Instance of class :mat:func:`ThermoPhase`
% :param model:
%     String indicating the transport model to use. Possible values
%     are ``'default'``, ``'Mix'``, and ``'Multi'``
% :param loglevel:
%     Level of diagnostic logging. Default if not specified is 4.
% :return:
%     Instance of class :mat:func:`Transport`
%

tr.id = 0;
if nargin == 1
    if isa(r, 'Transport')
        tr = r;
    else
        error(['Unless the first argument is an instance of class ' ...
              'Transport, there must be more than one argument'])
    end
else
    if nargin == 3
        loglevel = 4;
    end
    if ~isa(r, 'XML_Node')
        error(['The first argument must either be an instance of class ' ...
              'Transport or XML_Node'])
    elseif ~isa(th, 'ThermoPhase')
        error('The second argument must be an instance of class ThermoPhase')
    else
        tr.th = th;
        if strcmp(model, 'default')
            try
                node = child(r, 'transport');
                tr.model = attrib(node, 'model');
            catch
                tr.model = 'None';
            end
        else
            tr.model = model;
        end
        tr.id = trans_get(thermo_hndl(th), -1, tr.model, loglevel) ;
        tr = class(tr, 'Transport');
    end
end
