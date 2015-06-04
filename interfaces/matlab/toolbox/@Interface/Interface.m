function s = Interface(src, id, p1, p2, p3, p4)
% INTERFACE  Interface class constructor.
% s = Interface(src, id, p1, p2, p3, p4)
% See :ref:`sec-interfaces`.
%
% See also: :mat:func:`importEdge`, :mat:func:`importInterface`
%
% :param src:
%     CTI or CTML file containing the interface or edge phase.
% :param id:
%     Name of the interface or edge phase in the CTI or CTML file.
% :param p1:
%     Adjoining phase to the interface.
% :param p2:
%     Adjoining phase to the interface.
% :param p3:
%     Adjoining phase to the interface.
% :param p4:
%     Adjoining phase to the interface.
% :return:
%     Instance of class :mat:func:`Interface`
%

doc = XML_Node('doc', src);
node = findByID(doc, id);
t = ThermoPhase(node);
if nargin == 2
    k = Kinetics(node, t);
elseif nargin == 3
    k = Kinetics(node, t, p1);
elseif nargin == 4
    k = Kinetics(node, t, p1, p2);
elseif nargin == 5
    k = Kinetics(node, t, p1, p2, p3);
elseif nargin == 6
    k = Kinetics(node, t, p1, p2, p3, p4);
end

s.kin = k;
s.th = t;
s = class(s,'Interface', t, k);
