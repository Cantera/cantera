function tr = Transport(th, model, loglevel)
% TRANSPORT  Transport class constructor.
% tr = Transport(r, th, model, loglevel)
% Create a new instance of class :mat:func:`Transport`. One to three arguments
% may be supplied. The first must be an instance of class
% :mat:func:`ThermoPhase`. The second (optional) argument is the type of
% model desired, specified by the string ``'default'``, ``'Mix'`` or
% ``'Multi'``. ``'default'`` uses the default transport specified in the
% :mat:func:`XML_Node`. The third argument is the logging level desired.
%
% :param th:
%     Instance of class :mat:func:`ThermoPhase`
% :param model:
%     String indicating the transport model to use. Possible values
%     are ``'default'``, ``'None'``, ``'Mix'``, and ``'Multi'``.
%     Optional.
% :param loglevel:
%     Level of diagnostic logging. Default if not specified is 4.
% :return:
%     Instance of class :mat:func:`Transport`
%

tr.id = 0;
if nargin == 2
    model = 'default';
end

if nargin < 3
    loglevel = 4;
end

if ~isa(th, 'ThermoPhase')
    error('The first argument must be an instance of class ThermoPhase')
else
    tr.th = th;
    if strcmp(model, 'default')
        tr.id = trans_get(thermo_hndl(th), -2, loglevel);
    else
        tr.id = trans_get(thermo_hndl(th), -1, model, loglevel);
    end
    tr = class(tr, 'Transport');
end
