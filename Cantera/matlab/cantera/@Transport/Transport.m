function tr = Transport(model, th, loglevel)
%TRANSPORT  Transport class constructor.
%
%  k = TRANSPORT(model, p, loglevel) creates a transport
%  manager for phase object p. 
%
%  The 'model' parameter is a string that specifies the transport
%  model. The phase object must have already been created.
%
tr.id = 0;
if nargin == 3
  tr.th = th;
  tr.model = model;
  tr.id = trans_get(hndl(th), -1, model, loglevel) ;
  tr = class(tr,'Transport');
elseif isa(model,'Transport')
  tr = model;
else
  error('syntax error')
end

