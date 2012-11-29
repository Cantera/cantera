function x = FlowDevice(typ)
%
if nargin == 0
    typ = 1;
end
x.index = flowdevicemethods(0,typ);
if x.index < 0
    error(geterr);
end
x.type = typ;
x.upstream = -1;
x.downstream = -1;
x = class(x,'FlowDevice');
