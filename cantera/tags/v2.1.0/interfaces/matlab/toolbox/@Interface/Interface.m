function s = Interface(src, id, p1, p2, p3, p4)
% Interface - class Interface constructor.
%
doc = XML_Node('doc',src);
node = findByID(doc,id);
t = ThermoPhase(node);
if nargin == 2
    k = Kinetics(node,t);
elseif nargin == 3
    k = Kinetics(node,t,p1);
elseif nargin == 4
    k = Kinetics(node,t,p1,p2);
elseif nargin == 5
    k = Kinetics(node,t,p1,p2,p3);
elseif nargin == 6
    k = Kinetics(node,t,p1,p2,p3,p4);
end

s.kin = k;
s.th = t;
s = class(s,'Interface',t,k);
